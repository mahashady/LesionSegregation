from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

from utils.LAD_inference_utils import (
    run_binary_inference,
    run_joint_inference,
    save_binary_inference_results,
    save_joint_inference_results,
    plot_binary_q_with_ci,
    plot_binary_q_point_estimate,
    plot_pi_with_ci,
    plot_posterior_confidence,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="LAD inference runner")
    p.add_argument("--tumor_tsv", required=True, help="TSV with columns: sample, M, C")
    p.add_argument("--likelihood_tsv", required=True, help="TSV with columns: M, C, LAD, prob")
    p.add_argument("--output_dir", required=True, help="Output directory")
    p.add_argument("--prefix", required=True, help="Prefix for output files")

    p.add_argument("--mode", choices=["binary", "multistate"], default="binary")
    p.add_argument("--low_lads", default="1,2,3,4", help="Comma-separated LADs in the low group")
    p.add_argument("--min_prob", type=float, default=1e-300)
    p.add_argument("--tol", type=float, default=1e-10)
    p.add_argument("--max_iter", type=int, default=1000)

    p.add_argument("--filter_M_not_equal", type=int, default=None, help="Drop tumors with this M")
    p.add_argument("--bootstrap", action="store_true")
    p.add_argument("--n_boot", type=int, default=1000)
    p.add_argument("--ci", type=float, default=0.95)

    p.add_argument("--permutation", action="store_true")
    p.add_argument("--n_perm", type=int, default=5000)
    p.add_argument("--seed", type=int, default=42)

    return p.parse_args()


def run_lad_inference(
    tumor_df: pd.DataFrame,
    likelihood_df: pd.DataFrame,
    mode: str = "binary",
    low_lads: tuple[int, ...] = (1, 2),
    min_prob: float = 1e-300,
    tol: float = 1e-10,
    bootstrap: bool = False,
    n_boot: int = 1000,
    ci: float = 0.95,
    permutation: bool = False,
    n_perm: int = 5000,
    random_seed: int = 72,
    max_iter: int = 1000,
):
    if mode == "binary":
        return run_binary_inference(
            tumor_df=tumor_df,
            likelihood_df=likelihood_df,
            low_lads=list(low_lads),
            min_prob=min_prob,
            tol=tol,
            do_bootstrap=bootstrap,
            n_boot=n_boot,
            ci=ci,
            do_permutation=permutation,
            n_perm=n_perm,
            random_seed=random_seed,
        )

    if mode == "multistate":
        return run_joint_inference(
            tumor_df=tumor_df,
            likelihood_df=likelihood_df,
            d_values=None,
            max_iter=max_iter,
            tol=tol,
            min_prob=min_prob,
        )

    raise ValueError("mode must be 'binary' or 'multistate'")


def main() -> None:
    args = parse_args()

    tumor_df = pd.read_csv(args.tumor_tsv, sep="\t")
    likelihood_df = pd.read_csv(args.likelihood_tsv, sep="\t")

    if args.filter_M_not_equal is not None:
        tumor_df = tumor_df[tumor_df["M"] != args.filter_M_not_equal].copy()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    plots_dir = outdir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    low_lads = tuple(int(x) for x in args.low_lads.split(",") if x.strip())

    result = run_lad_inference(
        tumor_df=tumor_df,
        likelihood_df=likelihood_df,
        mode=args.mode,
        low_lads=low_lads,
        min_prob=args.min_prob,
        tol=args.tol,
        bootstrap=args.bootstrap,
        n_boot=args.n_boot,
        ci=args.ci,
        permutation=args.permutation,
        n_perm=args.n_perm,
        random_seed=args.seed,
        max_iter=args.max_iter,
    )

    if args.mode == "binary":
        save_binary_inference_results(result, outdir=outdir, prefix=args.prefix)

        fig = plot_binary_q_with_ci(result, outfile=plots_dir / f"{args.prefix}.binary_q_plot.png")
        plt.close(fig)

        fig = plot_binary_q_point_estimate(result,outfile=plots_dir / f"{args.prefix}.binary_q_point_estimate.png")
        plt.close(fig)

        print(result["binary_q_hat"])
        print(result["binary_posterior"].head())

        if "binary_q_CI" in result:
            print(result["binary_q_CI"])
        if "binary_permutation" in result:
            print(result["binary_permutation"])

    else:
        save_joint_inference_results(result, outdir=outdir, prefix=args.prefix)

        fig1 = plot_pi_with_ci(result, outfile=plots_dir / f"{args.prefix}.pi_plot.png")
        plt.close(fig1)

        fig2 = plot_posterior_confidence(result, outfile=plots_dir / f"{args.prefix}.posterior_confidence.png")
        plt.close(fig2)

        print(result["pi_df"])
        print(result["posterior_df"].head())


if __name__ == "__main__":
    main()