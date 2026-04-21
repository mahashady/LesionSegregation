#!/usr/bin/env python3
"""
Build likelihood lookup table for LAD inference.

This script simulates P(C | M, LAD) and saves it for reuse.

Outputs:
    - likelihood look up table (TSV)
    - metadata JSON

Modes:
    - no_recombination
    - recombination
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from utils.LAD_simul_utils import build_likelihood_df


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build likelihood lookup table for LAD inference")

    p.add_argument(
        "--outdir",
        default="/home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/LAD_inference",
        help="Output directory",
    )
    p.add_argument(
        "--prefix",
        default="likelihood_matrix_for_LAD_inference",
        help="Base output prefix",
    )

    p.add_argument(
        "--model",
        choices=["no_recombination", "recombination"],
        default="no_recombination",
        help="Simulation model",
    )

    p.add_argument(
        "--M_min",
        type=int,
        default=2,
        help="Minimum M value",
    )
    p.add_argument(
        "--M_max",
        type=int,
        default=15,
        help="Maximum M value (inclusive)",
    )

    p.add_argument(
        "--LAD_values",
        default="1,2,3,4,5",
        help="Comma-separated LAD states",
    )

    p.add_argument(
        "--n_iter",
        type=int,
        default=1000,
        help="Simulations per (M, LAD)",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed",
    )
    p.add_argument(
        "--pseudocount",
        type=float,
        default=0.0,
        help="Optional pseudocount smoothing",
    )
    p.add_argument(
        "--no_complete_grid",
        action="store_true",
        help="Disable complete C grid output",
    )

    # recombination-specific
    p.add_argument(
        "--recomb_lambda",
        type=float,
        default=2.0,
        help="Poisson mean number of recombinations per chromosome per division (used only in recombination mode)",
    )
    p.add_argument(
        "--win_size",
        type=int,
        default=10_000_000,
        help="Bin size for chromosome segmentation in recombination mode",
    )

    return p.parse_args()


def make_output_prefix(
    base_prefix: str,
    model: str,
    M_min: int,
    M_max: int,
    lad_values: list[int],
    n_iter: int,
    recomb_lambda: float | None = None,
    win_size: int | None = None,
) -> str:
    """
    Build informative output prefix.
    """
    lad_label = f"LAD{min(lad_values)}-{max(lad_values)}"
    m_label = f"M{M_min}-{M_max}"
    iter_label = f"n{n_iter}"

    if model == "recombination":
        return (
            f"{base_prefix}."
            f"{model}."
            f"{m_label}."
            f"{lad_label}."
            f"lambda{recomb_lambda}."
            f"win{win_size}."
            f"{iter_label}"
        )
    else:
        return (
            f"{base_prefix}."
            f"{model}."
            f"{m_label}."
            f"{lad_label}."
            f"{iter_label}"
        )


def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    M_values = range(args.M_min, args.M_max + 1)
    LAD_values = [int(x) for x in args.LAD_values.split(",") if x.strip()]
    complete_grid = not args.no_complete_grid

    out_prefix = make_output_prefix(
        base_prefix=args.prefix,
        model=args.model,
        M_min=args.M_min,
        M_max=args.M_max,
        lad_values=LAD_values,
        n_iter=args.n_iter,
        recomb_lambda=args.recomb_lambda if args.model == "recombination" else None,
        win_size=args.win_size if args.model == "recombination" else None,
    )

    out_file = outdir / f"{out_prefix}.tsv"
    meta_file = outdir / f"{out_prefix}.json"

    print("Building likelihood table...")
    print(f"Model: {args.model}")
    print(f"M range: {args.M_min}–{args.M_max}")
    print(f"LAD values: {LAD_values}")
    print(f"Iterations per (M,LAD): {args.n_iter}")
    if args.model == "recombination":
        print(f"Recombination lambda: {args.recomb_lambda}")
        print(f"Window size: {args.win_size}")
    print(f"Output: {out_file}")

    likelihood_df = build_likelihood_df(
        M_values=M_values,
        D_values=LAD_values,
        n_iter=args.n_iter,
        seed=args.seed,
        pseudocount=args.pseudocount,
        complete_grid=complete_grid,
        model=args.model,
        recomb_lambda=args.recomb_lambda if args.model == "recombination" else None,
        win_size=args.win_size,
    )

    # sanity check
    check = likelihood_df.groupby(["M", "LAD"])["prob"].sum()
    if not (check.round(6) == 1).all():
        raise ValueError("Probabilities do not sum to 1 within each (M, LAD) block.")

    likelihood_df.to_csv(out_file, sep="\t", index=False)

    metadata = {
        "model": args.model,
        "M_values": list(M_values),
        "LAD_values": LAD_values,
        "n_iter": args.n_iter,
        "seed": args.seed,
        "pseudocount": args.pseudocount,
        "complete_grid": complete_grid,
        "recomb_lambda": args.recomb_lambda if args.model == "recombination" else None,
        "win_size": args.win_size if args.model == "recombination" else None,
        "output_file": str(out_file),
    }

    with open(meta_file, "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"Saved {len(likelihood_df)} rows to {out_file}")
    print(f"Metadata saved to {meta_file}")


if __name__ == "__main__":
    main()