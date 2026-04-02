"""
inference.py

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
        Must contain columns:
        - tumor_id
        - M
        - C
    likelihood_df : pd.DataFrame
        Must contain columns:
        - M
        - C
        - D
        - prob
    d_values : list[int] or None
        If None, inferred from likelihood_df.
    min_prob : float
        Probability floor to avoid exact zeros in downstream EM.

    Returns
    -------
    pd.DataFrame
        Columns:
        - tumor_id
        - M
        - C
        - D
        - prob
    """
    req_t = {"tumor_id", "M", "C"}
    req_l = {"M", "C", "D", "prob"}

    if not req_t.issubset(tumor_df.columns):
        raise ValueError(f"tumor_df missing columns: {req_t - set(tumor_df.columns)}")
    if not req_l.issubset(likelihood_df.columns):
        raise ValueError(f"likelihood_df missing columns: {req_l - set(likelihood_df.columns)}")

    tumor_df = tumor_df.copy()
    likelihood_df = likelihood_df.copy()

    if d_values is None:
        d_values = sorted(likelihood_df["D"].unique().tolist())

    expanded = (
        tumor_df.assign(_tmp=1)
        .merge(pd.DataFrame({"D": d_values, "_tmp": 1}), on="_tmp")
        .drop(columns="_tmp")
    )

    merged = expanded.merge(
        likelihood_df[["M", "C", "D", "prob"]],
        on=["M", "C", "D"],
        how="left",
    )

    if merged["prob"].isna().any():
        bad = merged.loc[merged["prob"].isna(), ["tumor_id", "M", "C", "D"]]
        raise ValueError(
            "Some observed (M, C, D) combinations are missing from likelihood_df.\n"
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
    Build tumor-by-D likelihood matrix L where

        L[t, d] = P(C_t | M_t, D=d)

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
        d_values = sorted(long_df["D"].unique().tolist())

    L = (
        long_df.pivot(index="tumor_id", columns="D", values="prob")
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
        Columns:
        - tumor_id
        - M
        - C
        - post_D1, post_D2, ...
        - D_map
        - post_max
    """
    post_df = tumor_df[["tumor_id", "M", "C"]].copy().reset_index(drop=True)

    if responsibilities.shape[0] != len(post_df):
        raise ValueError("Number of rows in responsibilities does not match tumor_df.")

    for j, d in enumerate(d_values):
        post_df[f"post_D{d}"] = responsibilities[:, j]

    map_idx = responsibilities.argmax(axis=1)
    post_df["D_map"] = [d_values[i] for i in map_idx]
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
        "D": d_values,
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
    d_values = fit0["pi_df"]["D"].tolist()

    tumor_df = tumor_df.reset_index(drop=True).copy()
    T = len(tumor_df)

    boot_estimates = []

    for b in range(n_boot):
        idx = rng.integers(0, T, size=T)
        boot_df = tumor_df.iloc[idx].copy().reset_index(drop=True)

        # Make IDs unique within the bootstrap sample
        boot_df["tumor_id"] = [f"boot{b}_row{i}" for i in range(len(boot_df))]

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
        "D": d_values,
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
        Must contain:
        - tumor_id
        - M
        - C
    likelihood_df : pd.DataFrame
        Must contain:
        - M
        - C
        - D
        - prob
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
        - likelihood_long : long-form table with P(C_t | M_t, D)
        - likelihood_matrix : np.ndarray, shape (T, K)
        - d_values : list[int]
        - em_result : raw EM result dict
        - pi_df : pd.DataFrame with estimated pi
        - posterior_df : pd.DataFrame with per-tumor posterior probabilities
        - bootstrap_pi_df : pd.DataFrame or None
        - bootstrap_mat : np.ndarray or None
    """
    tumor_df = tumor_df.copy()

    if tumor_df["tumor_id"].duplicated().any():
        raise ValueError("tumor_df contains duplicated tumor_id values.")

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

def main(args):
    json_samples = ""
    df_multistat_by_chr = "/home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/ALL_Hartwig_bi_multi.by_chrom.txt"
    output_multi_summary = "/home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/multi_stat_by_chr.enriched.platinum.SBS17_less10.csv"
    with open(json_samples, "r") as f:
        data = json.load(f)
    list_of_samples = data["enriched"]

    tumor_df = build_tumor_summary(
        input_file: df_multistat_by_chr,
        output_file: output_multi_summary,
        samples_of_interest: list_of_samples,
        include_quatro: bool = True
    )

    result = run_joint_inference(
        tumor_df=tumor_df,
        likelihood_df=likelihood_df,
        bootstrap=True,
        n_boot=1000,
        random_seed=42,
    )

    print(result["pi_df"])
    print(result["posterior_df"].head())

    if result["bootstrap_pi_df"] is not None:
        print(result["bootstrap_pi_df"])
    

  

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)    