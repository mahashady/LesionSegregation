from __future__ import annotations

import numpy as np
import pandas as pd
import json

def build_tumor_summary(
    input_file: str,
    output_file: str | None = None,
    samples_of_interest: list[str] | None = None,
    include_quatro: bool = True,
    verbose: bool = True,
) -> pd.DataFrame:

    df = pd.read_csv(input_file, sep="\t")
    samples_of_interest[:] = [s + "T" for s in samples_of_interest]
    print(df.shape)
    print(samples_of_interest)
    # filter samples if needed
    if samples_of_interest is not None:
        df = df[df["sample"].isin(samples_of_interest)].copy()
    print(df.shape)
    # select multiallelic columns
    tri_cols = [c for c in df.columns if "triallelic" in c.lower() or "triallilic" in c.lower()]
    quad_cols = [c for c in df.columns if "quatro" in c.lower() or "qiuatro" in c.lower()]

    if not tri_cols:
        raise ValueError("No triallelic columns found.")

    multi_cols = tri_cols + quad_cols if include_quatro else tri_cols

    # compute per chromosome multiallelic
    df = df.copy()
    df["multi_per_chr"] = df[multi_cols].sum(axis=1)

    # aggregate per tumor
    tumor_df = (
        df.groupby("sample", as_index=False)
        .agg(
            M=("multi_per_chr", "sum"),
            C=("multi_per_chr", lambda x: int((x > 0).sum())),
        )
    )

    if output_file is not None:
        tumor_df.to_csv(output_file, sep="\t", index=False)

    return tumor_df

def main(args):
    json_samples = "../../results/samples_lists.Platinum.json"
    df_multistat_by_chr = "../../results/ALL_Hartwig_bi_multi.by_chrom.txt"
    output_multi_summary = "../../results/multi_stat_by_chr.enriched.platinum.SBS17_less10.csv"
    with open(json_samples, "r") as f:
        data = json.load(f)
    list_of_samples = data["enriched"]

    tumor_df = build_tumor_summary(
        input_file=df_multistat_by_chr,
        output_file=output_multi_summary,
        samples_of_interest= list_of_samples,
        include_quatro=True
    )

 

if __name__ == "__main__":
    main(None)