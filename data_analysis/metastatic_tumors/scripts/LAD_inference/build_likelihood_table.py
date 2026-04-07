"""
Build likelihood lookup table for LAD inference.

This script simulates P(C | M, D) and saves it for reuse.
C - number of chromosomes with multiallelic sites
M - number of multiallelic sites
D - LAD
"""

from pathlib import Path
import pandas as pd
import json

from utils.LAD_simul_utils import build_likelihood_df


def main() -> None:
    out_file = Path(
        "/home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/LAD_inference/likelihood_matrix_for_LAD_inference.csv"
    )

    out_file.parent.mkdir(parents=True, exist_ok=True)

    print("Building likelihood table...")
    print("M range: 2–15")
    print("LAD range: 1–7")
    print("Iterations per (M,LAD): 20000")

    likelihood_df = build_likelihood_df(
        M_values=range(2, 16),   # number of multiallelic sites
        D_values=[1, 2, 3, 4, 5, 6, 7],  # LAD states
        n_iter=20000,             # simulations per cell
        seed=42,                  # reproducibility
        pseudocount=0.0,
        complete_grid=True,
    )

    # sanity check
    check = likelihood_df.groupby(["M", "LAD"])["prob"].sum()
    if not (check.round(6) == 1).all():
        raise ValueError("Probabilities do not sum to 1!")

    likelihood_df.to_csv(out_file, sep="\t", index=False)

    # save metadata
    metadata = {
        "M_values": list(range(2, 16)),
        "D_values": [1, 2, 3, 4, 5],
        "n_iter": 20000,
        "seed": 42,
    }

    meta_file = out_file.with_suffix(".json")
    with open(meta_file, "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"Saved {len(likelihood_df)} rows to {out_file}")
    print(f"Metadata saved to {meta_file}")


if __name__ == "__main__":
    main()