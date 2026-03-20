#$ -S /usr/bin/python
import argparse
from glob import glob
import string
import gzip
import re
import os
import subprocess
import statistics
import json
import numpy as np
import random
import pandas as pd
from hmmlearn import hmm
from scipy.stats.stats import pearsonr   
from collections import Counter
from pathlib import Path

from hmm_models import run_HMM

def condense2bed(df, state_col):
    """
    Collapse consecutive rows with the same state into intervals (BED-style).
    
    Parameters
    ----------
    df : pd.DataFrame with info to condense
        Must have columns 'chr' and 'pos' and a state column.
    state_col : str
        Name of the column containing states.
    """

    df = df[~df[state_col].isna()]
    # Ensure data is sorted by chromosome and position
    df = df.sort_values(["chr", "pos"]).reset_index(drop=True)
    
    # Detect where state or chromosome changes
    state_change = (df[state_col] != df[state_col].shift(1)) | (df['chr'] != df['chr'].shift(1))
    
    # Assign a unique group id to each contiguous block
    group_id = state_change.cumsum()
    
    # Aggregate start and end positions per group
    result_df = df.groupby(['chr', group_id, state_col], as_index=False).agg(
        start_pos=('pos', 'first'),
        end_pos=('pos', 'last')
    )[['chr', 'start_pos', 'end_pos', state_col]]
    
    # Rename state column to 'state' for BED-style output
    #result_df = result_df.rename(columns={state_col: 'state'})
    
    return result_df


def compute_sites_and_states_props(states_intervals, df_all_sorted, state_col="state"):
    """
    Compute counts and proportions of each state for intervals and sites.

    Parameters
    ----------
    states_intervals : pd.DataFrame
        DataFrame containing intervals with assigned HMM states.
        Must have a column specified by `state_col` that contains the state assignment per interval.

    df_all_sorted : pd.DataFrame
        DataFrame containing per-site HMM state assignments.
        Must have a column specified by `state_col`.

    state_col : str, default "state"
        Name of the column in both DataFrames that contains the state assignments.
    """
    if state_col == "HMM_multi_state":
        all_states = ["M1","M2","M3"]
    elif state_col == "HMM_as_state":
        all_states = ["A1","A2","A3","A4","A5"]
    else:
        all_states = sorted(set(states_intervals[state_col]).union(df_all_sorted[state_col]))
  
    # Intervals
    intervals_counts = Counter(states_intervals[state_col])
    intervals_counts_full = {s: intervals_counts.get(s, 0) for s in all_states}
    print(intervals_counts_full)
    total_intervals = sum(intervals_counts_full.values())
    intervals_props = {k: v / total_intervals for k, v in intervals_counts_full.items()}
    print(intervals_props)

    # Sites
    sites_counts = Counter(df_all_sorted[state_col])
    sites_counts_full = {s: sites_counts.get(s, 0) for s in all_states}
    total_sites = sum(sites_counts_full.values())
    sites_props = {k: v / total_sites for k, v in sites_counts_full.items()}
    
    # Convert to Series or DataFrame for convenience
    intervals_props = pd.Series(intervals_props).sort_index()
    sites_props = pd.Series(sites_props).sort_index()
    
    return intervals_counts,intervals_props,sites_counts, sites_props

def summarize_hmm(intervals_df, sites_df, state_col, state_names, emissions,unique_categories, target_category):
    """
    Returns:
      row_values: list of values in the order matching generated columns
      columns: list of column names
    """
    # ---- find category index safely ----
    if target_category not in unique_categories:
        raise ValueError(f"{target_category} not found in categories {unique_categories}")

    cat_idx = list(unique_categories).index(target_category)

    # ---- compute counts/props ----
    (seg_counts, seg_props, site_counts, site_props) = compute_sites_and_states_props(
        intervals_df, sites_df, state_col=state_col
    )

    # ---- emissions for chosen category ----
    em_vals = [state_em[cat_idx] for state_em in emissions[:len(state_names)]]
    em_cols = [f"emission_{state_col}_{target_category}_{s}" for s in state_names]

    seg_count_vals = [seg_counts.get(s, 0) for s in state_names]
    seg_count_cols = [f"segments_{state_col}_{s}" for s in state_names]

    seg_prop_vals = [seg_props.get(s, 0) for s in state_names]
    seg_prop_cols = [f"segments_{state_col}_{s}_prop" for s in state_names]

    site_count_vals = [site_counts.get(s, 0) for s in state_names]
    site_count_cols = [f"sites_{state_col}_{s}" for s in state_names]

    site_prop_vals = [site_props.get(s, 0) for s in state_names]
    site_prop_cols = [f"sites_{state_col}_{s}_prop" for s in state_names]

    values = em_vals + seg_count_vals + seg_prop_vals + site_count_vals + site_prop_vals
    cols = em_cols + seg_count_cols + seg_prop_cols + site_count_cols + site_prop_cols

    return values, cols


def main():
    in_files = glob("../../data/encoded_data_all/*.hmm")
    out_multi_dir = Path("../results/HMM_multi_result")
    out_as_dir = Path("../results/HMM_as_result")
    out_multi_dir.mkdir(parents=True, exist_ok=True)
    out_as_dir.mkdir(parents=True, exist_ok=True)

    hmm_summary_rows = []
    summary_cols = None

    multi_state_names = ["M1", "M2", "M3"]
    as_state_names = ["A1", "A2", "A3", "A4", "A5"]
    unique_categories_multi = ["B", "M"]
    target_category_multi = "M"
    unique_categories_as = ["T>N", "A>N"]
    target_category_as = "T>N"
    for file_name in in_files:
        file_path = Path(file_name)
        sample = file_path.stem
        print(f"Processing {sample}")

        df = pd.read_csv(file_path, sep=",", low_memory=False)
        df = df.sort_values(["chr", "pos"], ascending=[True, True]).reset_index(drop=True)
        
        # ----- HMM_multi -----
        _, emissions_multi, _, Z_multi, _ = run_HMM(df[["code_multi"]].to_numpy(), "HMM_multi")

        df["HMM_multi_state"] = Z_multi
        multi_intervals = condense2bed(df,"HMM_multi_state")

        # ----- HMM_as (only rows where code_as == A>N or T>N) -----
        mask_as = df["code_as"].ne("other")
        df_as = df.loc[mask_as, ["chr", "pos", "code_as"]].copy()

        _, emissions_as, _, Z_as, _ = run_HMM(df_as[["code_as"]].to_numpy(), "HMM_as")
        df["HMM_as_state"] = pd.Series([pd.NA] * len(df), dtype="object")
        df.loc[mask_as, "HMM_as_state"] = Z_as
        # Important: build intervals using the filtered df to avoid NAs because of not A>N/T>N mutations
        df_as_for_bed = df.loc[mask_as, ["chr", "pos", "HMM_as_state"]].copy()
        as_intervals = condense2bed(df_as_for_bed, "HMM_as_state")

        # ----- Write per-sample outputs -----
        df.to_csv(out_multi_dir / f"{sample}.tsv", index=False, sep="\t")
        multi_intervals.to_csv(out_multi_dir / f"{sample}_intervals_HMMmulti.bed", index=False, sep="\t")
        as_intervals.to_csv(out_as_dir / f"{sample}_intervals_HMMas.bed", index=False, sep="\t")

        # ----- Summarize -----
        multi_vals, multi_cols = summarize_hmm(
            multi_intervals, df, "HMM_multi_state", multi_state_names, emissions_multi, unique_categories_multi, target_category_multi
        )
        as_vals, as_cols = summarize_hmm(
            as_intervals, df, "HMM_as_state", as_state_names, emissions_as, unique_categories_as, target_category_as
        )
        row = [sample] + multi_vals + as_vals
        cols = ["Sample"] + multi_cols + as_cols

        if summary_cols is None:
            summary_cols = cols

        hmm_summary_rows.append(row)

    df_hmm_summary = pd.DataFrame(hmm_summary_rows, columns=summary_cols)
    print(df_hmm_summary.head())

    outfile_summary = Path("../results/hmm_summary.tsv")
    df_hmm_summary.to_csv(outfile_summary, index=False, sep="\t")

 
if __name__=="__main__":
    main()