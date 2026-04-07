"""
Joint LAD inference across tumors using a mixture model.

Model
-----
For each tumor t with observed (M_t, C_t),

    P(C_t | M_t) = sum_D pi_D * P(C_t | M_t, D)

where:
- D is LAD state
- pi_D are the unknown population mixture proportions
- P(C | M, D) comes from likelihood_df

This module:
1. builds the tumor-by-D likelihood matrix from likelihood_df
2. fits pi by EM
3. computes per-tumor posterior probabilities P(D | C_t, M_t)
4. optionally computes bootstrap confidence intervals for pi
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import json
from pathlib import Path
import matplotlib.pyplot as plt

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
    req_t = {"sample", "M", "C"}
    req_l = {"M", "C", "LAD", "prob"}

    if not req_t.issubset(tumor_df.columns):
        raise ValueError(f"tumor_df missing columns: {req_t - set(tumor_df.columns)}")
    if not req_l.issubset(likelihood_df.columns):
        raise ValueError(f"likelihood_df missing columns: {req_l - set(likelihood_df.columns)}")

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
    print(merged.head())
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
    print(long_df.head(n=10))
    if d_values is None:
        d_values = sorted(long_df["LAD"].unique().tolist())

    L = (
        long_df.pivot(index="sample", columns="LAD", values="prob")
        .reindex(columns=d_values)
        .to_numpy()
    )

    return L, d_values, long_df


def em_fit_mixture_from_likelihoods(
    L: np.ndarray,
    max_iter: int = 1000,
    tol: float = 1e-10,
    pi_init: np.ndarray | None = None,
) -> dict:
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
    dict
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

    for i in range(max_iter):
        print("Iter ", i)
        print(pi)
        numerator = L * pi[None, :]
        denom = numerator.sum(axis=1, keepdims=True)
        responsibilities = numerator / denom
        pi_new = responsibilities.mean(axis=0)
        print(pi_new)
        print("_____________")

        loglik = float(np.sum(np.log(denom[:, 0])))
        loglik_trace.append(loglik)
        if len(loglik_trace) > 1:
            print(abs(loglik_trace[-1] - loglik_trace[-2]))
            if abs(loglik_trace[-1] - loglik_trace[-2]) < tol:
                pi = pi_new
                converged = True
                break

        pi = pi_new

    numerator = L * pi[None, :]
    denom = numerator.sum(axis=1, keepdims=True)
    responsibilities = numerator / denom

    return {
        "pi_hat": pi,
        "responsibilities": responsibilities,
        "loglik_trace": loglik_trace,
        "converged": converged,
        "n_iter": len(loglik_trace),
    }


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
        "ci_low": ci_low,
        "ci_high": ci_high,
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
        pi_hat=em_result["pi_hat"],
        d_values=d_values,
    )

    posterior_df = summarize_posteriors(
        tumor_df=tumor_df,
        responsibilities=em_result["responsibilities"],
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

def plot_pi_with_ci(result: dict, outfile: str | Path | None = None):
    pi_df = result["pi_df"].copy()
    boot_df = result["bootstrap_pi_df"]

    fig, ax = plt.subplots(figsize=(5, 4))

    if boot_df is not None:
        plot_df = pi_df.merge(boot_df[["LAD", "ci_low", "ci_high"]], on="LAD", how="left")
        yerr = [
            plot_df["pi_hat"] - plot_df["ci_low"],
            plot_df["ci_high"] - plot_df["pi_hat"],
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


def main():
    tumors_summary_file = "/home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/multi_stat_by_chr.enriched.platinum.SBS17b_prop_lt_0.1.csv"
    likelihood_file = "/home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/LAD_inference/likelihood_matrix_for_LAD_inference.csv"
    outdir = Path("/home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/LAD_inference")
    prefix = "joint_LAD_inference.platinum.SBS17b_prop_lt_0.1"
    tumor_df = pd.read_csv(tumors_summary_file, sep="\t")
    tumor_df = tumor_df[tumor_df["M"] != 1]
    print(tumor_df.head())

    likelihood_df = pd.read_csv(likelihood_file, sep="\t")
    print(likelihood_df.head())

    result = run_joint_inference(
        tumor_df=tumor_df,
        likelihood_df=likelihood_df,
        bootstrap=False,
        n_boot=1000,
        random_seed=42,
    )

    save_joint_inference_results(result, outdir=outdir, prefix=prefix)

    plots_dir = outdir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    fig1 = plot_pi_with_ci(result, outfile=plots_dir / f"{prefix}.pi_plot.png")
    plt.close(fig1)

    fig2 = plot_posterior_confidence(result, outfile=plots_dir / f"{prefix}.posterior_confidence.png")
    plt.close(fig2)

    print(result["pi_df"])
    print(result["posterior_df"].head())

    if result["bootstrap_pi_df"] is not None:
        print(result["bootstrap_pi_df"])


  

if __name__=="__main__":
    main()    