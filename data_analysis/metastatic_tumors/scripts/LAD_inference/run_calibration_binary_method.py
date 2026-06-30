#!/usr/bin/env python3
"""
Binary LAD calibration simulation.

Simulates cohorts with known early/late LAD ratio and checks
whether binary inference recovers it.

Outputs:
  - TSV with results
  - boxplot PNG
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils.LAD_inference_utils import run_binary_inference
from utils.calibration_utils import *

def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--tumor_file", required=True)
    p.add_argument("--likelihood_file", required=True)
    p.add_argument("--outdir", required=True)

    p.add_argument("--q_values", default="0.1,0.25,0.5,0.75,0.9")
    p.add_argument("--n_reps", type=int, default=20)

    p.add_argument("--early_lads", default="1,2,3,4")
    p.add_argument("--late_lads", default="5,6,7,8,9,10")

    p.add_argument("--seed", type=int, default=72)

    return p.parse_args()


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("======================================")
    print("Binary LAD calibration simulation")
    print("======================================")

    print("\nLoading data...")
    tumor_df = pd.read_csv(args.tumor_file, sep="\t")
    tumor_df = tumor_df[tumor_df["M"] > 1]
    likelihood_df = pd.read_csv(args.likelihood_file, sep="\t")

    q_values = [float(x) for x in args.q_values.split(",")]
    early_lads = [int(x) for x in args.early_lads.split(",")]
    late_lads = [int(x) for x in args.late_lads.split(",")]

    print(f"  Tumors loaded: {len(tumor_df)}")
    print(f"  q_values: {q_values}")
    print(f"  reps per q: {args.n_reps}")
    print(f"  early LADs: {early_lads}")
    print(f"  late LADs: {late_lads}")

    print("\nStarting calibration...\n")

    results = run_binary_calibration(
        observed_tumor_df=tumor_df,
        likelihood_df=likelihood_df,
        true_q_values=q_values,
        n_reps=args.n_reps,
        early_lads=early_lads,
        late_lads=late_lads,
        random_seed=args.seed,
    )

    print("\nCalibration finished.")
    print(tumor_df.head())
    observed_fit = run_binary_inference(
        tumor_df=tumor_df[["sample", "M", "C"]],
        likelihood_df=likelihood_df,
        low_lads=early_lads,
        do_bootstrap=True,
    )
    observed_q_hat = float(observed_fit["binary_q_hat"]["q_hat"].iloc[0])
    observed_ci_df = observed_fit.get("binary_q_CI", None)
    print(observed_ci_df)
    if observed_ci_df is not None:
        observed_ci_low = float(observed_ci_df["ci_early"].iloc[0])
        observed_ci_high = float(observed_ci_df["ci_late"].iloc[0])

    low_lads_str = "-".join(map(str, early_lads))
    late_lads_str = "-".join(map(str, late_lads))
    prefix = f"binary_calibration.lowLAD_{low_lads_str}.lateLAD_{late_lads_str}"

    out_tsv = outdir / f"{prefix}.tsv"
    print(f"Saving results saved to {out_tsv}")
    results.to_csv(out_tsv, sep="\t", index=False)

    plot_file = outdir / f"{prefix}.png"
    print(f"Generating plot saved to {plot_file}")
    plot_binary_calibration(results, outfile=plot_file, observed_q_hat=observed_q_hat)

    plot_file = outdir / f"{prefix}.for_paper.png"
    print(f"Generating plot saved to {plot_file}")
    plot_binary_calibration_publication(results, outfile=plot_file, observed_q_hat=observed_q_hat, observed_ci=(observed_ci_low, observed_ci_high))

    print("\nDone.")


if __name__ == "__main__":
    main()