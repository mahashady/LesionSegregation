#!/usr/bin/env python3
"""
Compute LAD distribution without EM.

This script:
- computes per-tumor posterior P(LAD | data)
- averages across tumors
- outputs final LAD distribution
- optionally plots it
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def compute_parametric_posterior(
    tumor_df: pd.DataFrame,
    likelihood_df: pd.DataFrame,
    d_values: list[int] | None = None,
) -> pd.DataFrame:
    """
    Compute LAD distribution by averaging per-tumor normalized likelihoods.

    Returns
    -------
    pd.DataFrame with columns:
        LAD, prob
    """

    # build L matrix 
    long_df = likelihood_df.merge(
        tumor_df[["sample", "M", "C"]],
        on=["M", "C"],
        how="inner",
    )

    if d_values is None:
        d_values = sorted(long_df["LAD"].unique())

    L = (
        long_df.pivot(index="sample", columns="LAD", values="prob")
        .reindex(columns=d_values)
        .fillna(0.0)
        .to_numpy()
    )

    # normalize per tumor (row-wise)
    row_sums = L.sum(axis=1, keepdims=True)
    if np.any(row_sums == 0):
        raise ValueError("Some tumors have zero likelihood across all LADs")

    posterior = L / row_sums

    # average across tumors
    pi_hat = posterior.mean(axis=0)

    return pd.DataFrame({
        "LAD": d_values,
        "prob": pi_hat,
    })

def plot_lad_distribution(df, outfile=None):
    fig, ax = plt.subplots(figsize=(5, 4))

    ax.bar(df["LAD"], df["prob"])

    ax.set_xlabel("LAD")
    ax.set_ylabel("Probability")
    ax.set_title("Posterior LAD distribution")

    ax.set_ylim(0, df["prob"].max() * 1.2)

    fig.tight_layout()

    if outfile:
        fig.savefig(outfile, dpi=300)

    return fig

def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--tumor_file", required=True)
    p.add_argument("--likelihood_file", required=True)
    p.add_argument("--outdir", required=True)

    p.add_argument("--prefix", default="lad_no_em")

    return p.parse_args()


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("\nLoading data...")
    tumor_df = pd.read_csv(args.tumor_file, sep="\t")
    likelihood_df = pd.read_csv(args.likelihood_file, sep="\t")

    print(f"  Tumors: {len(tumor_df)}")

    print("\nComputing parametric posterior ...")

    posterior_df = compute_parametric_posterior(
        tumor_df=tumor_df,
        likelihood_df=likelihood_df,
    )

    print("\nResult:")
    print(posterior_df)

    # save table
    out_tsv = outdir / f"{args.prefix}.tsv"
    posterior_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nSaved table to {out_tsv}")

    # plot
    out_plot = outdir / f"{args.prefix}.png"
    plot_lad_distribution(posterior_df, out_plot)

    print("\nDone.")


if __name__ == "__main__":
    main()