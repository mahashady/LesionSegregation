from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Basic helpers
# -----------------------------

def require_columns(df: pd.DataFrame, cols: List[str], name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{name} missing required columns: {missing}. Present: {list(df.columns)}")


def normalize_prob_vector(x: np.ndarray) -> np.ndarray:
    s = float(np.sum(x))
    if s <= 0:
        raise ValueError("Cannot normalize: sum <= 0")
    return x / s

def vprint(verbose: bool, *args, **kwargs):
    if verbose:
        print(*args, **kwargs)

# -----------------------------
# Likelihood handling
# -----------------------------

def build_likelihood_long(
    tumor_df: pd.DataFrame,
    likelihood_df: pd.DataFrame,
    d_values: list[int] | None = None,
    min_prob: float = 1e-300,
) -> pd.DataFrame:
    """
    Expand observed tumor data across D states and attach P(C_t | M_t, D).

    Parameters
    ----------
    tumor_df : pd.DataFrame
        Must contain columns: sample, M, C, 
    likelihood_df : pd.DataFrame
        Must contain columns: M, C, LAD, prob
    d_values : list[int] or None
        If None, inferred from likelihood_df.
    min_prob : float
        Probability floor to avoid exact zeros in downstream EM.

    Returns
    -------
    pd.DataFrame
        Columns: sample, M, C, LAD, prob
    """
    require_columns(tumor_df, ["sample", "M", "C"], "tumor_df")
    require_columns(likelihood_df, ["M", "C", "LAD", "prob"], "likelihood_df")

    tumor_df = tumor_df.copy()
    likelihood_df = likelihood_df.copy()

    if d_values is None:
        d_values = sorted(likelihood_df["LAD"].unique().tolist())

    expanded = (
        tumor_df.assign(_tmp=1)
        .merge(pd.DataFrame({"LAD": d_values, "_tmp": 1}), on="_tmp")
        .drop(columns="_tmp")
    )

    merged = expanded.merge(
        likelihood_df[["M", "C", "LAD", "prob"]],
        on=["M", "C", "LAD"],
        how="left",
    )
    if merged["prob"].isna().any():
        bad = merged.loc[merged["prob"].isna(), ["sample", "M", "C", "LAD"]]
        raise ValueError(
            "Some observed (M, C, LAD) combinations are missing from likelihood_df.\n"
            f"Examples:\n{bad.head(10)}"
        )

    merged["prob"] = merged["prob"].clip(lower=min_prob)
    return merged


def build_likelihood_matrix(
    tumor_df: pd.DataFrame,
    likelihood_df: pd.DataFrame,
    d_values: list[int] | None = None,
    min_prob: float = 1e-300,
) -> tuple[np.ndarray, list[int], pd.DataFrame]:
    """
    Build tumor-by-LAD likelihood matrix L where

        L[t, d] = P(C_t | M_t, D=d)
    Rows correspond to samples, columns to possible LAD values.

    Returns
    -------
    L : np.ndarray
        Shape (T, K)
    d_values : list[int]
        Ordered LAD states
    long_df : pd.DataFrame
        Long-form merged table
    """
    long_df = build_likelihood_long(
        tumor_df=tumor_df,
        likelihood_df=likelihood_df,
        d_values=d_values,
        min_prob=min_prob,
    )
    if d_values is None:
        d_values = sorted(long_df["LAD"].unique().tolist())
    else:
        d_values = list(d_values)

    pivot = long_df.pivot_table(index="sample", columns="LAD", values="prob", aggfunc="first")
    pivot = pivot.reindex(index=tumor_df["sample"].tolist(), columns=d_values)

    if pivot.isna().any().any():
        raise RuntimeError("NaNs present in likelihood matrix after pivot/reindex.")

    L = pivot.to_numpy(dtype=float)

    return L, d_values, long_df

# -----------------------------
# Multistate EM
# -----------------------------

@dataclass
class EMResult:
    pi_hat: np.ndarray
    responsibilities: np.ndarray
    loglik_trace: List[float]
    converged: bool
    n_iter: int
    
def em_fit_mixture_from_likelihoods(
    L: np.ndarray,
    max_iter: int = 1000,
    tol: float = 1e-10,
    pi_init: np.ndarray | None = None,
) -> EMResult:
    """
    Fit mixture proportions pi from known component likelihoods by EM.

    Model:
        p(data_t) = sum_D pi_D * L[t, D]

    Parameters
    ----------
    L : np.ndarray
        Shape (T, K), where L[t, k] = P(observation_t | D_k)
    max_iter : int
        Maximum number of EM iterations
    tol : float
        Convergence threshold on absolute change in log-likelihood
    pi_init : np.ndarray or None
        Optional initial pi vector of shape (K,)

    Returns
    -------
    object of class EMResult
        - pi_hat: np.ndarray, shape (K,)
        - responsibilities: np.ndarray, shape (T, K)
        - loglik_trace: list[float]
        - converged: bool
        - n_iter: int
    """
    L = np.asarray(L, dtype=float)

    if L.ndim != 2:
        raise ValueError("L must be a 2D array.")

    T, K = L.shape

    if np.any(L < 0):
        raise ValueError("L must be nonnegative.")

    if np.any(L.sum(axis=1) <= 0):
        raise ValueError("Each tumor must have positive likelihood under at least one D.")

    if pi_init is None:
        pi = np.full(K, 1.0 / K)
    else:
        pi = np.asarray(pi_init, dtype=float)
        if pi.shape != (K,):
            raise ValueError(f"pi_init must have shape ({K},)")
        if np.any(pi < 0):
            raise ValueError("pi_init must be nonnegative.")
        if pi.sum() <= 0:
            raise ValueError("pi_init must sum to a positive value.")
        pi = pi / pi.sum()

    loglik_trace: list[float] = []
    converged = False

    for _ in range(max_iter):
        numerator = L * pi[None, :]
        denom = numerator.sum(axis=1, keepdims=True)
        responsibilities = numerator / denom
        pi_new = responsibilities.mean(axis=0)

        loglik = float(np.sum(np.log(denom[:, 0])))
        loglik_trace.append(loglik)
        if len(loglik_trace) > 1:
            if abs(loglik_trace[-1] - loglik_trace[-2]) < tol:
                pi = pi_new
                converged = True
                break

        pi = pi_new

    numerator = L * pi[None, :]
    denom = numerator.sum(axis=1, keepdims=True)
    responsibilities = numerator / denom

    return EMResult(
        pi_hat=pi,
        responsibilities=responsibilities,
        loglik_trace=loglik_trace,
        converged=converged,
        n_iter=len(loglik_trace),
    )


def summarize_posteriors(
    tumor_df: pd.DataFrame,
    responsibilities: np.ndarray,
    d_values: list[int],
) -> pd.DataFrame:
    """
    Convert posterior responsibility matrix to a tumor-level summary table.

    Returns
    -------
    pd.DataFrame
        Columns: sample, M, C, post_D1 (post_D2, ...), D_map, post_max
    """
    post_df = tumor_df[["sample", "M", "C"]].copy().reset_index(drop=True)

    if responsibilities.shape[0] != len(post_df):
        raise ValueError("Number of rows in responsibilities does not match tumor_df.")

    for j, d in enumerate(d_values):
        post_df[f"post_LAD{d}"] = responsibilities[:, j]

    map_idx = responsibilities.argmax(axis=1)
    post_df["LAD_map"] = [d_values[i] for i in map_idx]
    post_df["post_max"] = responsibilities.max(axis=1)

    return post_df


def summarize_pi(
    pi_hat: np.ndarray,
    d_values: list[int],
) -> pd.DataFrame:
    """
    Return the estimated LAD mixture proportions in a tidy table.
    """
    return pd.DataFrame({
        "LAD": d_values,
        "pi_hat": np.asarray(pi_hat, dtype=float),
    })

def bootstrap_pi(
    tumor_df: pd.DataFrame,
    likelihood_df: pd.DataFrame,
    d_values: list[int] | None = None,
    n_boot: int = 500,
    ci: float = 0.95,
    random_seed: int = 1,
    max_iter: int = 1000,
    tol: float = 1e-10,
    min_prob: float = 1e-300,
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Bootstrap confidence intervals for pi by resampling tumors with replacement.

    Parameters
    ----------
    tumor_df : pd.DataFrame
        Observed tumor summaries
    likelihood_df : pd.DataFrame
        Lookup table with P(C | M, D)
    d_values : list[int] or None
        LAD states
    n_boot : int
        Number of bootstrap samples
    ci : float
        Confidence level, e.g. 0.95
    random_seed : int
        RNG seed
    max_iter : int
        Max EM iterations
    tol : float
        EM convergence threshold
    min_prob : float
        Probability floor when merging

    Returns
    -------
    summary_df : pd.DataFrame
        Columns:
        - D
        - pi_hat
        - boot_se
        - ci_low
        - ci_high
        - n_boot_success
    boot_mat : np.ndarray
        Bootstrap estimates, shape (n_successful_boot, K)
    """
    rng = np.random.default_rng(random_seed)

    fit0 = run_joint_inference(
        tumor_df=tumor_df,
        likelihood_df=likelihood_df,
        d_values=d_values,
        max_iter=max_iter,
        tol=tol,
        min_prob=min_prob,
        bootstrap=False,
    )

    pi_hat = fit0["pi_df"]["pi_hat"].to_numpy()
    d_values = fit0["pi_df"]["LAD"].tolist()

    tumor_df = tumor_df.reset_index(drop=True).copy()
    T = len(tumor_df)

    boot_estimates = []

    for b in range(n_boot):
        idx = rng.integers(0, T, size=T)
        boot_df = tumor_df.iloc[idx].copy().reset_index(drop=True)

        # Make IDs unique within the bootstrap sample
        boot_df["sample"] = [f"boot{b}_row{i}" for i in range(len(boot_df))]

        try:
            fit_b = run_joint_inference(
                tumor_df=boot_df,
                likelihood_df=likelihood_df,
                d_values=d_values,
                max_iter=max_iter,
                tol=tol,
                min_prob=min_prob,
                bootstrap=False,
            )
            boot_estimates.append(fit_b["pi_df"]["pi_hat"].to_numpy())
        except Exception:
            continue

    if len(boot_estimates) == 0:
        raise RuntimeError("All bootstrap fits failed.")

    boot_mat = np.vstack(boot_estimates)

    alpha = 1 - ci
    q_low = 100 * (alpha / 2)
    q_high = 100 * (1 - alpha / 2)

    ci_low = np.percentile(boot_mat, q_low, axis=0)
    ci_high = np.percentile(boot_mat, q_high, axis=0)
    boot_se = boot_mat.std(axis=0, ddof=1)

    summary_df = pd.DataFrame({
        "LAD": d_values,
        "pi_hat": pi_hat,
        "boot_se": boot_se,
        "ci_early": ci_low,
        "ci_late": ci_high,
        "n_boot_success": boot_mat.shape[0],
    })

    return summary_df, boot_mat

def run_joint_inference(
    tumor_df: pd.DataFrame,
    likelihood_df: pd.DataFrame,
    d_values: list[int] | None = None,
    max_iter: int = 1000,
    tol: float = 1e-10,
    min_prob: float = 1e-300,
    bootstrap: bool = False,
    n_boot: int = 500,
    ci: float = 0.95,
    random_seed: int = 1,
) -> dict:
    """
    Main function to run joint LAD inference.

    Parameters
    ----------
    tumor_df : pd.DataFrame
        Must contain: sample, M (number of multi-allelic sites), C (number of chroms)
    likelihood_df : pd.DataFrame
        Must contain: M, C, LAD, prob
    d_values : list[int] or None
        LAD states to use. If None, inferred from likelihood_df.
    max_iter : int
        Max EM iterations
    tol : float
        EM convergence threshold
    min_prob : float
        Probability floor when building likelihood matrix
    bootstrap : bool
        Whether to compute bootstrap confidence intervals for pi
    n_boot : int
        Number of bootstrap samples if bootstrap=True
    ci : float
        Confidence level for bootstrap intervals
    random_seed : int
        RNG seed for bootstrap

    Returns
    -------
    dict
        - likelihood_long : long-form table with P(C_t | M_t, LAD)
        - likelihood_matrix : np.ndarray, shape (T, K)
        - d_values : list[int]
        - em_result : raw EM result dict
        - pi_df : pd.DataFrame with estimated pi
        - posterior_df : pd.DataFrame with per-tumor posterior probabilities
        - bootstrap_pi_df : pd.DataFrame or None
        - bootstrap_mat : np.ndarray or None
    """
    tumor_df = tumor_df.copy()

    if tumor_df["sample"].duplicated().any():
        raise ValueError("tumor_df contains duplicated sample values.")


    L, d_values, likelihood_long = build_likelihood_matrix(
        tumor_df=tumor_df,
        likelihood_df=likelihood_df,
        d_values=d_values,
        min_prob=min_prob,
    )

    em_result = em_fit_mixture_from_likelihoods(
        L=L,
        max_iter=max_iter,
        tol=tol,
    )

    pi_df = summarize_pi(
        pi_hat=em_result.pi_hat,
        d_values=d_values,
    )

    posterior_df = summarize_posteriors(
        tumor_df=tumor_df,
        responsibilities=em_result.responsibilities,
        d_values=d_values,
    )

    bootstrap_pi_df = None
    bootstrap_mat = None

    if bootstrap:
        bootstrap_pi_df, bootstrap_mat = bootstrap_pi(
            tumor_df=tumor_df,
            likelihood_df=likelihood_df,
            d_values=d_values,
            n_boot=n_boot,
            ci=ci,
            random_seed=random_seed,
            max_iter=max_iter,
            tol=tol,
            min_prob=min_prob,
        )

    return {
        "likelihood_long": likelihood_long,
        "likelihood_matrix": L,
        "d_values": d_values,
        "em_result": em_result,
        "pi_df": pi_df,
        "posterior_df": posterior_df,
        "bootstrap_pi_df": bootstrap_pi_df,
        "bootstrap_mat": bootstrap_mat,
    }

# -----------------------------
# Binary inference
# -----------------------------

@dataclass
class BinaryFit:
    q_hat: float
    post_low: np.ndarray
    post_high: np.ndarray
    loglik: float

def golden_section_maximize(
    f,
    a: float = 0.0,
    b: float = 1.0,
    tol: float = 1e-10,
    max_iter: int = 5000,
) -> float:
    gr = (np.sqrt(5.0) + 1.0) / 2.0
    c = b - (b - a) / gr
    d = a + (b - a) / gr
    fc = f(c)
    fd = f(d)

    it = 0
    while (b - a) > tol and it < max_iter:
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - (b - a) / gr
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + (b - a) / gr
            fd = f(d)
        it += 1

    return (a + b) / 2.0


def build_group_likelihoods(
    tumor_df: pd.DataFrame,
    likelihood_df: pd.DataFrame,
    low_lads: List[int] = [1, 2, 3, 4],
    min_prob: float = 1e-300,
) -> Tuple[np.ndarray, np.ndarray, List[int]]:
    L, d_values, _ = build_likelihood_matrix(
        tumor_df=tumor_df,
        likelihood_df=likelihood_df,
        d_values=None,
        min_prob=min_prob,
    )

    d_values = list(d_values)
    low_set = set(low_lads)
    high_lads = [d for d in d_values if d not in low_set]
    if len(high_lads) == 0:
        raise ValueError("Late LAD set is empty. Check low_lads or LAD states present in likelihood_df.")

    lad_to_j = {d: j for j, d in enumerate(d_values)}
    low_js = [lad_to_j[d] for d in low_lads if d in lad_to_j]
    if len(low_js) == 0:
        raise ValueError(f"None of low_lads={low_lads} found in LAD states: {d_values}")
    high_js = [lad_to_j[d] for d in high_lads]

    L_low = np.mean(L[:, low_js], axis=1)
    L_high = np.mean(L[:, high_js], axis=1)

    if np.any(L_low <= 0) or np.any(L_high <= 0):
        raise RuntimeError("Non-positive group likelihood encountered; check min_prob flooring.")

    return L_low, L_high, d_values

def fit_binary_mixture(L_low: np.ndarray, L_high: np.ndarray, tol: float = 1e-10) -> BinaryFit:
    L_low = np.asarray(L_low, dtype=float)
    L_high = np.asarray(L_high, dtype=float)
    if L_low.shape != L_high.shape:
        raise ValueError("L_early and L_late must have the same shape.")
    if np.any(L_low <= 0) or np.any(L_high <= 0):
        raise ValueError("Likelihoods must be >0. Check min_prob flooring.")

    def loglik(q: float) -> float:
        mix = q * L_low + (1.0 - q) * L_high
        return float(np.sum(np.log(mix)))

    q_hat = golden_section_maximize(loglik, a=0.0, b=1.0, tol=tol)

    denom = q_hat * L_low + (1.0 - q_hat) * L_high
    post_low = (q_hat * L_low) / denom
    post_high = 1.0 - post_low

    return BinaryFit(q_hat=q_hat, post_low=post_low, post_high=post_high, loglik=loglik(q_hat))

def bootstrap_q(
    tumor_df: pd.DataFrame,
    likelihood_df: pd.DataFrame,
    low_lads: List[int],
    n_boot: int = 2000,
    ci: float = 0.95,
    random_seed: int = 1,
    min_prob: float = 1e-300,
    tol: float = 1e-10,
) -> pd.DataFrame:
    rng = np.random.default_rng(random_seed)
    tumor_df = tumor_df.reset_index(drop=True).copy()
    T = len(tumor_df)

    q_vals = np.empty(n_boot, dtype=float)

    for b in range(n_boot):
        idx = rng.integers(0, T, size=T)
        boot_df = tumor_df.iloc[idx].copy().reset_index(drop=True)
        boot_df["sample"] = [f"boot{b}_row{i}" for i in range(T)]

        L_low, L_high, _ = build_group_likelihoods(
            tumor_df=boot_df,
            likelihood_df=likelihood_df,
            low_lads=low_lads,
            min_prob=min_prob,
        )
        fit = fit_binary_mixture(L_low, L_high, tol=tol)
        q_vals[b] = fit.q_hat

    alpha = 1.0 - ci
    ci_low = float(np.quantile(q_vals, alpha / 2.0))
    ci_high = float(np.quantile(q_vals, 1.0 - alpha / 2.0))

    return pd.DataFrame({
        "q_hat_bootstrap_mean": [float(np.mean(q_vals))],
        "q_hat_bootstrap_sd": [float(np.std(q_vals, ddof=1))],
        "ci_early": [ci_low],
        "ci_late": [ci_high],
        "n_boot": [int(n_boot)],
        "ci_level": [float(ci)],
    })


def permutation_test_q(
    tumor_df: pd.DataFrame,
    likelihood_df: pd.DataFrame,
    low_lads: List[int],
    n_perm: int = 5000,
    random_seed: int = 1,
    min_prob: float = 1e-300,
    tol: float = 1e-10,
) -> pd.DataFrame:
    rng = np.random.default_rng(random_seed)
    tumor_df = tumor_df.reset_index(drop=True).copy()

    L_low_obs, L_high_obs, _ = build_group_likelihoods(
        tumor_df=tumor_df,
        likelihood_df=likelihood_df,
        low_lads=low_lads,
        min_prob=min_prob,
    )
    fit_obs = fit_binary_mixture(L_low_obs, L_high_obs, tol=tol)
    q_obs = fit_obs.q_hat

    q_perm = np.empty(n_perm, dtype=float)

    m_vals = tumor_df["M"].astype(int).to_numpy()
    groups: Dict[int, np.ndarray] = {m: np.where(m_vals == m)[0] for m in np.unique(m_vals)}

    for b in range(n_perm):
        perm_df = tumor_df.copy()
        for m, idx in groups.items():
            perm_df.loc[idx, "C"] = rng.permutation(perm_df.loc[idx, "C"].to_numpy())

        L_low, L_high, _ = build_group_likelihoods(
            tumor_df=perm_df,
            likelihood_df=likelihood_df,
            low_lads=low_lads,
            min_prob=min_prob,
        )
        fit_b = fit_binary_mixture(L_low, L_high, tol=tol)
        q_perm[b] = fit_b.q_hat

    p_one_sided = float(np.mean(q_perm >= q_obs))

    return pd.DataFrame({
        "q_obs": [float(q_obs)],
        "p_value_one_sided": [p_one_sided],
        "n_perm": [int(n_perm)],
        "q_perm_mean": [float(np.mean(q_perm))],
        "q_perm_sd": [float(np.std(q_perm, ddof=1))],
    })

def run_binary_inference(
    tumor_df: pd.DataFrame,
    likelihood_df: pd.DataFrame,
    low_lads: List[int],
    min_prob: float = 1e-300,
    tol: float = 1e-10,
    do_bootstrap: bool = False,
    n_boot: int = 200,
    ci: float = 0.95,
    do_permutation: bool = False,
    n_perm: int = 1000,
    random_seed: int = 72,
    verbose: bool = True,
) -> Dict[str, pd.DataFrame]:

    vprint(verbose, f"[Binary] Building likelihoods...")
    L_low, L_high, d_values = build_group_likelihoods(
        tumor_df=tumor_df,
        likelihood_df=likelihood_df,
        low_lads=low_lads,
        min_prob=min_prob,
    )
    
    vprint(verbose, f"[Binary] Fitting mixture...")
    fit = fit_binary_mixture(L_low, L_high, tol=tol)

    posterior = tumor_df.loc[:, ["sample", "M", "C"]].copy().reset_index(drop=True)
    posterior["post_early_LAD"] = fit.post_low
    posterior["post_late_LAD"] = fit.post_high
    posterior["class_map"] = np.where(posterior["post_early_LAD"] >= 0.5, "early LAD", "late LAD")
    posterior["post_max"] = posterior[["post_early_LAD", "post_late_LAD"]].max(axis=1)

    out: Dict[str, pd.DataFrame] = {}
    out["binary_q_hat"] = pd.DataFrame({
        "q_hat": [fit.q_hat],
        "loglik": [fit.loglik],
        "low_lads": [",".join(map(str, low_lads))],
        "high_lads": [",".join(map(str, [d for d in d_values if d not in set(low_lads)]))],
    })
    out["binary_posterior"] = posterior

    if do_bootstrap:
        vprint(verbose, f"[Binary] Running bootstrap (n={n_boot})...")
        out["binary_q_CI"] = bootstrap_q(
            tumor_df=tumor_df,
            likelihood_df=likelihood_df,
            low_lads=low_lads,
            n_boot=n_boot,
            ci=ci,
            random_seed=random_seed,
            min_prob=min_prob,
            tol=tol,
        )

    if do_permutation:
        vprint(verbose, f"[Binary] Running permutation test (n={n_perm})...")
        out["binary_permutation"] = permutation_test_q(
            tumor_df=tumor_df,
            likelihood_df=likelihood_df,
            low_lads=low_lads,
            n_perm=n_perm,
            random_seed=random_seed,
            min_prob=min_prob,
            tol=tol,
        )

    return out

# -----------------------------
# Save / plot helpers
# -----------------------------

def save_joint_inference_results(result: dict, outdir: str | Path, prefix: str) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    result["pi_df"].to_csv(outdir / f"{prefix}.pi.tsv", sep="\t", index=False)
    result["posterior_df"].to_csv(outdir / f"{prefix}.posterior.tsv", sep="\t", index=False)
    result["likelihood_long"].to_csv(outdir / f"{prefix}.likelihood_long.tsv", sep="\t", index=False)

    if result["bootstrap_pi_df"] is not None:
        result["bootstrap_pi_df"].to_csv(
            outdir / f"{prefix}.pi_bootstrap.tsv",
            sep="\t",
            index=False
        )

def save_binary_inference_results(result: dict, outdir: str | Path, prefix: str) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    result["binary_q_hat"].to_csv(outdir / f"{prefix}.binary_q_hat.tsv", sep="\t", index=False)
    result["binary_posterior"].to_csv(outdir / f"{prefix}.binary_posterior.tsv", sep="\t", index=False)

    if "binary_q_CI" in result:
        result["binary_q_CI"].to_csv(outdir / f"{prefix}.binary_q_CI.tsv", sep="\t", index=False)

    if "binary_permutation" in result:
        result["binary_permutation"].to_csv(outdir / f"{prefix}.binary_permutation.tsv", sep="\t", index=False)

def plot_pi_with_ci(result: dict, outfile: str | Path | None = None):
    pi_df = result["pi_df"].copy()
    boot_df = result["bootstrap_pi_df"]

    fig, ax = plt.subplots(figsize=(5, 4))

    if boot_df is not None:
        plot_df = pi_df.merge(boot_df[["LAD", "ci_early", "ci_late"]], on="LAD", how="left")
        yerr = [
            plot_df["pi_hat"] - plot_df["ci_early"],
            plot_df["ci_late"] - plot_df["pi_hat"],
        ]
        ax.bar(plot_df["LAD"].astype(str), plot_df["pi_hat"])
        ax.errorbar(
            x=range(len(plot_df)),
            y=plot_df["pi_hat"],
            yerr=yerr,
            fmt="none",
            capsize=4,
        )
    else:
        plot_df = pi_df
        ax.bar(plot_df["LAD"].astype(str), plot_df["pi_hat"])

    ax.set_xlabel("LAD")
    ax.set_ylabel("Estimated proportion")
    ax.set_title("Joint inferred LAD distribution")
    ax.set_ylim(0, max(0.05, plot_df["pi_hat"].max() * 1.25))
    fig.tight_layout()

    if outfile is not None:
        fig.savefig(outfile, dpi=300, bbox_inches="tight")

    return fig

def plot_binary_q_point_estimate(
    result: dict,
    outfile: str | Path | None = None,
):
    """
    Point estimate plot for inferred Early LAD proportion q_hat.
    """

    q_hat = float(result["binary_q_hat"]["q_hat"].iloc[0])
    ci_df = result.get("binary_q_CI", None)

    fig, ax = plt.subplots(figsize=(3.2, 4.2))

    # Point estimate
    ax.errorbar(
        x=[0],
        y=[q_hat],
        yerr=None,
        fmt="o",
        color="black",
        markersize=6,
        capsize=0,
        zorder=3,
    )

    # Confidence interval
    if ci_df is not None:
        ci_low = float(ci_df["ci_early"].iloc[0])
        ci_high = float(ci_df["ci_late"].iloc[0])

        yerr = [[q_hat - ci_low], [ci_high - q_hat]]

        ax.errorbar(
            x=[0],
            y=[q_hat],
            yerr=yerr,
            fmt="o",
            color="black",
            ecolor="black",
            elinewidth=1.5,
            capsize=5,
            capthick=1.5,
            markersize=6,
            zorder=4,
        )

    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(0, 1.1)

    ax.set_xticks([0])
    ax.set_xticklabels(
        ["Early LAD (LAD 1-4)\nproportion"],
        fontsize=11,
    )

    ax.set_ylabel("Estimated proportion", fontsize=12)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.tick_params(axis="y", labelsize=11)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(axis="y", linestyle="--", alpha=0.25)

    fig.tight_layout()

    if outfile is not None:
        fig.savefig(outfile, dpi=300, bbox_inches="tight")

    return fig

def plot_binary_q_with_ci(result: dict, outfile: str | Path | None = None):

    q_hat = float(result["binary_q_hat"]["q_hat"].iloc[0])
    ci_df = result.get("binary_q_CI", None)

    fig, ax = plt.subplots(figsize=(3.5, 4))

    x = [0, 1]
    heights = [q_hat, 1.0 - q_hat]

    ax.bar(
        x,
        heights,
        width=0.5,
        color=["#4C72B0", "#DD8452"],
        edgecolor="black",
        linewidth=0.8,
    )

    # Confidence interval for Early group
    if ci_df is not None:

        ci_low = float(ci_df["ci_early"].iloc[0])
        ci_high = float(ci_df["ci_late"].iloc[0])

        yerr = [[q_hat - ci_low], [ci_high - q_hat]]

        ax.errorbar(
            x=[0],
            y=[q_hat],
            yerr=yerr,
            fmt="none",
            color="black",
            linewidth=1.2,
            capsize=4,
            capthick=1.2,
            zorder=10,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(
        ["Early\n(LAD 1–4)", "Late\n(LAD 5+)"],
        fontsize=11,
    )

    ax.set_ylabel(
        "Estimated proportion",
        fontsize=12,
    )

    ax.set_ylim(0, 1.1)
    ax.set_yticks(np.arange(0, 1.01, 0.2))

    ax.tick_params(
        axis="y",
        labelsize=11,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.3,
    )

    fig.tight_layout()

    if outfile is not None:
        fig.savefig(
            outfile,
            dpi=300,
            bbox_inches="tight",
        )

    return fig

def plot_posterior_confidence(result: dict, outfile: str | Path | None = None):
    posterior_df = result["posterior_df"].copy()

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.hist(posterior_df["post_max"], bins=20)
    ax.set_xlabel("Max posterior probability per tumor")
    ax.set_ylabel("Number of tumors")
    ax.set_title("Confidence of tumor LAD assignments")
    fig.tight_layout()

    if outfile is not None:
        fig.savefig(outfile, dpi=300, bbox_inches="tight")

    return fig





