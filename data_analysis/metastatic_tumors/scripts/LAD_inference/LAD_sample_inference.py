"""
sample_specific.py

Sample-specific LAD inference for tumors with high number of multiallelic sites.

Implements:
- simulate_divisions (same logic as original R)
- per-sample simulation of C under different LAD states
- extraction of P(C_obs | D) likelihoods
- plotting
- main function to run the workflow
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


# ------------------------------------------------------------
# Core simulation (same logic as R)
# ------------------------------------------------------------

def simulate_divisions(input_vector, N, rng=None):
    if rng is None:
        rng = np.random.default_rng()

    vec = np.array(input_vector, dtype=int).copy()

    if N == 0:
        return vec

    for _ in range(N):
        for i in range(len(vec)):
            current_value = vec[i]

            if current_value == 2:
                vec[i] = rng.choice([2, 1, 0], p=[1/4, 1/2, 1/4])
            elif current_value == 1:
                vec[i] = rng.choice([1, 0], p=[1/2, 1/2])
            else:
                vec[i] = 0

    # preserve original behavior
    if vec.sum() == 0:
        vec[0] = 1

    return vec


# ------------------------------------------------------------
# Simulation for ONE sample
# ------------------------------------------------------------

def permute_mult_sample_table(
    data: pd.DataFrame,
    sample_id: str,
    n_iter: int = 1000,
    rng: np.random.Generator | None = None,
) -> pd.DataFrame:

    if rng is None:
        rng = np.random.default_rng()

    iv = np.repeat(2, 22)

    local_data = data.loc[data["sample"] == sample_id].copy()
    if local_data.empty:
        raise ValueError(f"Sample {sample_id} not found.")

    dist_rows = []

    # observed values
    C_obs = int((local_data["MultA"] > 0.1).sum())
    M_obs = int(local_data["MultA"].sum())

    for mrca in range(5):  # 0..4
        D = mrca + 1

        for _ in range(n_iter):

            result_vector = simulate_divisions(iv, mrca, rng=rng)

            denom = result_vector.sum()
            if denom == 0:
                continue

            probabilities = (
                local_data["BiA"].to_numpy(dtype=float)
                / local_data["Total_mutations"].to_numpy(dtype=float)
                * (result_vector / denom)
            )

            if probabilities.sum() <= 0:
                continue

            probabilities = probabilities / probabilities.sum()

            selected_chr = rng.choice(
                np.arange(22),
                size=M_obs,
                replace=True,
                p=probabilities
            )

            C_sim = int(np.unique(selected_chr).size)

            dist_rows.append({
                "C_sim": C_sim,
                "C_obs": C_obs,
                "D": D,
                "M": M_obs,
            })

    return pd.DataFrame(dist_rows)


# ------------------------------------------------------------
# Convert simulation → likelihoods
# ------------------------------------------------------------

def sample_specific_likelihoods(dist: pd.DataFrame) -> pd.DataFrame:
    """
    Convert simulation table into likelihood P(C_obs | D)
    """

    if dist.empty:
        raise ValueError("Empty simulation table")

    C_obs = int(dist["C_obs"].iloc[0])

    # empirical distribution
    tab = (
        dist.groupby("D")["C_sim"]
        .value_counts(normalize=True)
        .rename("prob")
        .reset_index()
    )

    tab = tab.rename(columns={"C_sim": "C"})

    # keep only observed C
    out = tab[tab["C"] == C_obs].copy()

    # ensure all D present
    all_D = pd.DataFrame({"D": sorted(dist["D"].unique())})
    out = all_D.merge(out, on="D", how="left")

    out["prob"] = out["prob"].fillna(0.0)
    out["C_obs"] = C_obs

    return out[["D", "C_obs", "prob"]]


# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------

def plot_sample(dist: pd.DataFrame, sample_id: str, savepath=None):
    states = sorted(dist["D"].unique(), reverse=True)
    C_obs = dist["C_obs"].iloc[0]

    fig, axes = plt.subplots(
        nrows=len(states),
        figsize=(5, 7),
        sharex=True
    )

    if len(states) == 1:
        axes = [axes]

    max_c = int(dist["C_sim"].max())
    bins = np.arange(0.5, max_c + 1.5, 1)

    for ax, D in zip(axes, states):
        sub = dist.loc[dist["D"] == D, "C_sim"]

        ax.hist(sub, bins=bins, density=True, alpha=0.7)
        ax.axvline(C_obs, linewidth=1.5)
        ax.set_ylabel(f"LAD {D}")

    axes[-1].set_xlabel("N chroms with multiallelic sites")
    fig.suptitle(sample_id)
    fig.tight_layout()

    if savepath is not None:
        fig.savefig(savepath, dpi=300, bbox_inches="tight")

    return fig


# ------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------

def run_sample_specific(
    data: pd.DataFrame,
    sample_id: str,
    n_iter: int = 2000,
    outdir: str | Path | None = None,
    seed: int = 1,
):
    """
    Main function to run sample-specific LAD inference.

    Returns
    -------
    dict with:
      - dist: simulation table
      - likelihoods: P(C_obs | D)
      - fig: matplotlib figure
    """

    rng = np.random.default_rng(seed)

    dist = permute_mult_sample_table(
        data=data,
        sample_id=sample_id,
        n_iter=n_iter,
        rng=rng,
    )

    likelihoods = sample_specific_likelihoods(dist)

    savepath = None
    if outdir is not None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        savepath = outdir / f"{sample_id}_sample_specific.png"

    fig = plot_sample(dist, sample_id, savepath)

    return {
        "dist": dist,
        "likelihoods": likelihoods,
        "fig": fig,
    }
def main(args):
    CHR_OE = load_chr_oe(
    "../results/ALL_Hartwig_bi_multi.txt",
    "../results/ALL_Hartwig_bi_multi.by_chrom.txt",
    )
    result = run_sample_specific(
        data=CHR_OE,
        sample_id="CPCT02020393",
        n_iter=2000,
        outdir="plots/",
    )

    print(result["likelihoods"])    
    

  

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)        