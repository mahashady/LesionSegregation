from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils.LAD_inference_utils import run_binary_inference


def sample_C_from_lookup(
    M: int,
    LAD: int,
    likelihood_df: pd.DataFrame,
    rng: np.random.Generator,
) -> int:
    """
    Sample one C from lookup table P(C | M, LAD).
    """
    sub = likelihood_df[(likelihood_df["M"] == M) & (likelihood_df["LAD"] == LAD)]
    if sub.empty:
        raise ValueError(f"No lookup rows found for M={M}, LAD={LAD}")

    probs = sub["prob"].to_numpy(dtype=float)
    C_vals = sub["C"].to_numpy(dtype=int)

    probs = probs / probs.sum()
    return int(rng.choice(C_vals, p=probs))


def simulate_dataset_from_lookup(
    observed_tumor_df: pd.DataFrame,
    likelihood_df: pd.DataFrame,
    q_early: float,
    early_lads: Sequence[int] = (1, 2, 3, 4),
    late_lads: Sequence[int] = (5, 6, 7, 8, 9, 10),
    rng: np.random.Generator | None = None,
    within_group_probs_early: Sequence[float] | None = None,
    within_group_probs_late: Sequence[float] | None = None,
) -> pd.DataFrame:
    """
    Simulate one cohort with:
      - same number of tumors as observed_tumor_df
      - same exact M values as observed_tumor_df
      - early/late LAD ratio = q_early / (1-q_early)

    Returns
    -------
    pd.DataFrame
        Columns: sample, M, C, LAD_true, group_true
    """
    if rng is None:
        rng = np.random.default_rng()

    df = observed_tumor_df.copy().reset_index(drop=True)

    if not {"sample", "M"}.issubset(df.columns):
        raise ValueError("observed_tumor_df must contain columns: sample, M")

    early_lads = list(early_lads)
    late_lads = list(late_lads)

    if within_group_probs_early is None:
        within_group_probs_early = np.repeat(1 / len(early_lads), len(early_lads))
    else:
        within_group_probs_early = np.asarray(within_group_probs_early, dtype=float)
        within_group_probs_early = within_group_probs_early / within_group_probs_early.sum()

    if within_group_probs_late is None:
        within_group_probs_late = np.repeat(1 / len(late_lads), len(late_lads))
    else:
        within_group_probs_late = np.asarray(within_group_probs_late, dtype=float)
        within_group_probs_late = within_group_probs_late / within_group_probs_late.sum()

    n = len(df)
    group_true = rng.choice(["early", "late"], size=n, p=[q_early, 1 - q_early])

    LAD_true = []
    C_sim = []

    for i in range(n):
        M_i = int(df.loc[i, "M"])

        if group_true[i] == "early":
            lad_i = int(rng.choice(early_lads, p=within_group_probs_early))
        else:
            lad_i = int(rng.choice(late_lads, p=within_group_probs_late))

        c_i = sample_C_from_lookup(M=M_i, LAD=lad_i, likelihood_df=likelihood_df, rng=rng)

        LAD_true.append(lad_i)
        C_sim.append(c_i)

    out = pd.DataFrame({
        "sample": [f"sim_{i+1}" for i in range(n)],
        "M": df["M"].astype(int).to_numpy(),
        "C": np.asarray(C_sim, dtype=int),
        "LAD_true": np.asarray(LAD_true, dtype=int),
        "group_true": group_true,
    })
    return out


def run_binary_calibration(
    observed_tumor_df: pd.DataFrame,
    likelihood_df: pd.DataFrame,
    true_q_values: Iterable[float],
    n_reps: int = 20,
    early_lads: Sequence[int] = (1, 2, 3, 4),
    late_lads: Sequence[int] = (5, 6, 7, 8, 9, 10),
    inference_low_lads: Sequence[int] | None = None,
    random_seed: int = 72,
    min_prob: float = 1e-300,
    tol: float = 1e-10,
    bootstrap: bool = False,
    permutation: bool = False,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    For each true early ratio q, simulate n_reps cohorts and run binary inference.

    Returns
    -------
    pd.DataFrame
        One row per simulated cohort with true and inferred ratios.
    """
    if inference_low_lads is None:
        inference_low_lads = list(early_lads)

    rng = np.random.default_rng(random_seed)
    rows = []

    true_q_values = list(true_q_values)
    total_runs = len(true_q_values) * n_reps
    run_counter = 0

    if verbose:
        print("Starting binary calibration")
        print(f"  Number of true q values: {len(true_q_values)}")
        print(f"  Replicates per q: {n_reps}")
        print(f"  Total simulation runs: {total_runs}")
        print(f"  Early LADs: {list(early_lads)}")
        print(f"  Late LADs: {list(late_lads)}")

    for q_true in true_q_values:
        if verbose:
            print(f"\n=== True early ratio q = {q_true:.3f} ===")

        for rep in range(n_reps):
            run_counter += 1
            if verbose:
                print(f"[{run_counter}/{total_runs}] replicate {rep + 1}/{n_reps}")

            sim_seed = int(rng.integers(0, 2**32 - 1))
            sim_rng = np.random.default_rng(sim_seed)

            sim_df = simulate_dataset_from_lookup(
                observed_tumor_df=observed_tumor_df,
                likelihood_df=likelihood_df,
                q_early=q_true,
                early_lads=early_lads,
                late_lads=late_lads,
                rng=sim_rng,
            )

            fit = run_binary_inference(
                tumor_df=sim_df[["sample", "M", "C"]],
                likelihood_df=likelihood_df,
                low_lads=list(inference_low_lads),
                min_prob=min_prob,
                tol=tol,
                do_bootstrap=bootstrap,
                n_boot=200,
                ci=0.95,
                do_permutation=permutation,
                n_perm=1000,
                random_seed=sim_seed,
            )

            q_hat = float(fit["binary_q_hat"]["q_hat"].iloc[0])

            if verbose:
                print(
                    f"    q_hat={q_hat:.3f} | "
                    f"n_early_true={(sim_df['group_true'] == 'early').sum()} | "
                    f"n_late_true={(sim_df['group_true'] == 'late').sum()}"
                )

            rows.append({
                "q_true": float(q_true),
                "rep": rep + 1,
                "q_hat": q_hat,
                "n_tumors": len(sim_df),
                "mean_M": sim_df["M"].mean(),
                "mean_C": sim_df["C"].mean(),
                "n_early_true": int((sim_df["group_true"] == "early").sum()),
                "n_late_true": int((sim_df["group_true"] == "late").sum()),
            })

    if verbose:
        print("\nBinary calibration finished.")

    return pd.DataFrame(rows)


def plot_binary_calibration(
    calib_df: pd.DataFrame,
    outfile: str | Path | None = None,
    observed_q_hat: float | None = None,
):
    """
    Boxplot:
      x = true early ratio
      y = inferred early ratio
    """
    q_levels = sorted(calib_df["q_true"].unique())
    data = [calib_df.loc[calib_df["q_true"] == q, "q_hat"].to_numpy() for q in q_levels]

    fig, ax = plt.subplots(figsize=(7, 5))

    positions = q_levels
    ax.boxplot(data, positions=positions, widths=0.05)

    # overlay points
    for pos, q in zip(positions, q_levels):
        y = calib_df.loc[calib_df["q_true"] == q, "q_hat"].to_numpy()
        x = np.random.normal(loc=pos, scale=0.01, size=len(y))
        ax.plot(x, y, "o", alpha=0.8,markersize=2,color="#5B4B8A")

    # diagonal reference in q-space
    # map true q values to their box positions
    ax.plot(positions, q_levels, linestyle="--", color="grey")

    if observed_q_hat is not None:
        ax.axhline(
            y=observed_q_hat,
            linestyle="--",
            linewidth=2,
            label=f"Observed q̂ = {observed_q_hat:.2f}",
        )
        ax.legend(loc="lower right")

    ax.set_xticks(positions)
    ax.set_xticklabels([str(q) for q in q_levels])
    ax.set_xlabel("True early LAD proportion")
    ax.set_ylabel("Inferred early LAD proportion")
    ax.set_title("Binary LAD calibration")
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlim(-0.1, 1.1)

    fig.tight_layout()

    if outfile is not None:
        fig.savefig(outfile, dpi=300, bbox_inches="tight")

    return fig

def plot_binary_calibration_publication(
    calib_df: pd.DataFrame,
    outfile: str | Path | None = None,
    observed_q_hat: float | None = None,
    observed_ci: tuple[float, float] | None = None,    
):
    """
    Publication-style calibration plot:
      x = true early LAD proportion
      y = inferred early LAD proportion

    Shows:
      - individual simulation replicates as light points
      - median q_hat per true q
      - 95% empirical interval
      - diagonal y = x reference line
      - optional observed q_hat as horizontal line
    """

    plot_df = calib_df.copy()
    q_levels = np.array(sorted(plot_df["q_true"].unique()), dtype=float)

    summary = (
        plot_df
        .groupby("q_true")["q_hat"]
        .agg(
            median="median",
            q025=lambda x: np.quantile(x, 0.025),
            q975=lambda x: np.quantile(x, 0.975),
        )
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=(5.5, 4.2))

    rng = np.random.default_rng(1)

    # Individual replicates
    for q in q_levels:
        y = plot_df.loc[plot_df["q_true"] == q, "q_hat"].to_numpy()
        x = rng.normal(loc=q, scale=0.008, size=len(y))

        ax.scatter(
            x,
            y,
            s=12,
            color="black",
            alpha=0.35,
            linewidths=0,
            zorder=2,
        )

    # 95% empirical intervals
    yerr = np.vstack([
        summary["median"] - summary["q025"],
        summary["q975"] - summary["median"],
    ])

    ax.errorbar(
        summary["q_true"],
        summary["median"],
        yerr=yerr,
        fmt="o",
        color="black",
        ecolor="black",
        elinewidth=1.4,
        capsize=4,
        capthick=1.4,
        markersize=5,
        zorder=4,
        label="Median ± 95% interval",
    )

    # Calibration diagonal
    ax.plot(
        [0, 1],
        [0, 1],
        linestyle="--",
        color="red",
        linewidth=1.2,
        zorder=1,
        label="Perfect calibration",
    )

    # Observed estimate
    if observed_q_hat is not None:
        if observed_ci is not None:
            ci_low, ci_high = observed_ci

            ax.errorbar(
                observed_q_hat,
                observed_q_hat,
                xerr=[
                    [observed_q_hat - ci_low],
                    [ci_high - observed_q_hat],
                ],
                yerr=[
                    [observed_q_hat - ci_low],
                    [ci_high - observed_q_hat],
                ],
                fmt="o",
                markersize=9,
                color="#D55E00",
                ecolor="#D55E00",
                markeredgecolor="black",
                markeredgewidth=0.8,
                elinewidth=1.6,
                capsize=4,
                capthick=1.6,
                zorder=6,
                label=f"Observed q̂ = {observed_q_hat:.2f} 95% CI",
            )
        else:
            ax.scatter(
                observed_q_hat,
                observed_q_hat,
                s=120,
                color="#D55E00",
                edgecolor="black",
                linewidth=0.8,
                zorder=6,
                label=f"Observed q̂ = {observed_q_hat:.2f}",
            )

    ax.set_xlabel("Simulated early LAD proportion", fontsize=12)
    ax.set_ylabel("Inferred early LAD proportion", fontsize=12)

    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)

    ax.set_xticks(np.arange(0, 1.01, 0.25))
    ax.set_yticks(np.arange(0, 1.01, 0.25))

    ax.tick_params(axis="both", labelsize=11)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(
        axis="both",
        linestyle="--",
        alpha=0.25,
        linewidth=0.8,
    )

    ax.legend(
        frameon=False,
        fontsize=9,
        loc="lower right",
    )

    fig.tight_layout()

    if outfile is not None:
        fig.savefig(outfile, dpi=300, bbox_inches="tight")

    return fig