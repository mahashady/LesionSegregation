#$ -S /usr/bin/python

import argparse
import glob
import string
import gzip
import re
import os
import subprocess
import statistics
import json
import random
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.special import gammaln

def multinomial_loglik(counts, probs):
    counts = np.asarray(counts, dtype=float)
    probs = np.asarray(probs, dtype=float)

    eps = 1e-12
    probs = np.clip(probs, eps, 1.0)

    return np.sum(counts * np.log(probs))

def infer_lad_row(counts, expected_dict, m1_prop, threshold=0.05):
    # Override rule
    if m1_prop < threshold:
        return 1, None
    if m1_prop > (1-threshold):
        return "late", None
    loglik = {
        mrca: multinomial_loglik(counts, probs)
        for mrca, probs in expected_dict.items()
    }
    best_mrca = max(loglik, key=loglik.get)
    return best_mrca, loglik


def infer_lad_euclid_row(props, expected_dict, m1_prop, threshold=0.05):
    # Override rule
    if m1_prop < threshold:
        return 1
    if m1_prop > (1-threshold):
        return "late"
    dists = {
        mrca: np.linalg.norm(props - probs)
        for mrca, probs in expected_dict.items()
    }
    return min(dists, key=dists.get)

def no_semicolon_and_not_empty(s: pd.Series) -> pd.Series:
        s = s.fillna("")
        return (~s.str.contains(";", regex=False)) & (s != "")

def main(args):

    expected = {
        "1": [0, 0, 1],
        "2": [1/4, 2/4, 1/4],
        "3": [9/16, 6/16, 1/16],
        "4": [49/64, 14/64, 1/64],
        "5": [225/256, 30/256, 1/256]
    }

    # ----------------------------
    # 1) Load input tables
    # ----------------------------
    # HMM summary 
    df_hmm_summary = pd.read_csv("../results/hmm_summary.tsv", sep="\t")
    count_cols_segm = ["segments_HMM_multi_state_M1", "segments_HMM_multi_state_M2", "segments_HMM_multi_state_M3"]
    prop_cols_segm = ["segments_HMM_multi_state_M1_prop", "segments_HMM_multi_state_M2_prop", "segments_HMM_multi_state_M3_prop"]
    count_cols_sites = ["sites_HMM_multi_state_M1", "sites_HMM_multi_state_M2", "sites_HMM_multi_state_M3"]
    prop_cols_sites = ["sites_HMM_multi_state_M1_prop", "sites_HMM_multi_state_M2_prop", "sites_HMM_multi_state_M3_prop"]
    # Ploidy summary
    df_ploidy_summary = pd.read_csv("../results/Ploidy_summary.tsv", sep=r"\s+", header=0)

    # Mice lines summary: C3H + CAST
    c3h_df = pd.read_csv(
        "../../data/C3H_summary_with_sign_driver_positions.txt",
        sep=",",
        header=0,
    )
    c3h_df["driver_positions_c3h_mutId"] = c3h_df["driver_positions"]
    c3h_df = c3h_df.rename(columns={c3h_df.columns[0]: "Sample"})
    c3h_df["mice_line"] = "C3H"

    cast_df = pd.read_csv(
        "../../data/CAST_summary_with_sign_driver_positions.txt",
        sep=",",
        header=0,
    )
    cast_df = cast_df.rename(columns={cast_df.columns[0]: "Sample"})
    cast_df["mice_line"] = "CAST"

    df_mice_lines_summary = pd.concat([c3h_df, cast_df], ignore_index=True)

    # ----------------------------
    # 2) Merge tables
    # ----------------------------
    print(df_hmm_summary.head())
    print(df_ploidy_summary.head())

    full_table = df_hmm_summary.merge(df_ploidy_summary, on="Sample", how="inner")
    full_table = full_table.merge(df_mice_lines_summary, on="Sample", how="inner")

    # ----------------------------
    # 3) Ploidy_by_HMM_as and mixtures
    # ----------------------------
    full_table["Ploidy_by_HMM_as"] = "one_tumour"

    as_props = full_table[[f"sites_HMM_as_state_A{i}_prop" for i in range(1, 6)]]
    # tumors with >=4 states having prop >0.05 are mixtures
    full_table.loc[(as_props > 0.05).sum(axis=1) >= 4, "Ploidy_by_HMM_as"] = "mixture_tumours"
    print(len(full_table[full_table["Ploidy_by_HMM_as"] == "mixture_tumours"]))
    full_table["tumour_mixture"] = 0
    mask_mixtures = ((full_table["Mixture_tumours"] == "mixture_tumours") & (full_table["subsets_fit"] == "Good")) | (
        full_table["Ploidy_by_HMM_as"] == "mixture_tumours"
        )
    full_table.loc[mask_mixtures, "tumour_mixture"] = 1

    # ----------------------------
    # 4) Symmetry classification
    # ----------------------------
    full_table["Symmetry"] = "Asymmetric"


    # If the most extreme asymmetry states are absent OR
    # most extreme asymmetry states have low emissions OR
    # opposite asymmetries in clonal and subclonal mutations and good fit of vaf in clonal and subclonal
    sym_mask = (
        ((full_table["sites_HMM_as_state_A1_prop"] < 0.05) & (full_table["sites_HMM_as_state_A5_prop"] < 0.05))
        | (full_table["emission_HMM_as_state_T>N_A1"] < 0.8)
        | (full_table["emission_HMM_as_state_T>N_A5"] > 0.2)
        | ((full_table["Mixture_tumours"] == "Symmetric") & (full_table["subsets_fit"] == "Good"))
        )
    print(sym_mask)
    full_table.loc[sym_mask, "Symmetry"] = "Symmetric"

    # ----------------------------
    # 5) Export mixtures subset
    # ----------------------------
    df_mixtures = full_table.loc[full_table["tumour_mixture"] == 1].copy()
    df_mixtures.to_csv("../results/Summary_mixtures.txt", index=False)
    print("Mixed tumors info has been written to ../results/Summary_mixtures.txt")


    # ----------------------------
    # 6) Divisions subset and LAD assignment
    # ----------------------------
    df_no_mixtures = full_table.loc[full_table["tumour_mixture"] != 1].copy()

    df_no_mixtures["LAD_ML_segm"] = df_no_mixtures.apply(
        lambda row: infer_lad_row(
            row[count_cols_segm].values, 
            expected, 
            row["segments_HMM_multi_state_M1_prop"])[0],
        axis=1
    )   

    df_no_mixtures["LAD_EUCLID_segm"] = df_no_mixtures.apply(
        lambda row: infer_lad_euclid_row(
            row[prop_cols_segm].values.astype(float),
            expected,
            row["segments_HMM_multi_state_M1_prop"]),
        axis=1
    )
    df_no_mixtures["LAD_ML_sites"] = df_no_mixtures.apply(
        lambda row: infer_lad_row(
            row[count_cols_sites].values, 
            expected, 
            row["sites_HMM_multi_state_M1_prop"])[0],
        axis=1
    )   

    df_no_mixtures["LAD_EUCLID_sites"] = df_no_mixtures.apply(
        lambda row: infer_lad_euclid_row(
            row[prop_cols_sites].values.astype(float),
            expected,
            row["sites_HMM_multi_state_M1_prop"]),
        axis=1
    )

    lad_cols = [
        "LAD_ML_segm",
        "LAD_EUCLID_segm",
        "LAD_ML_sites",
        "LAD_EUCLID_sites",
    ]

    df_no_mixtures.loc[
        df_no_mixtures["Symmetry"] == "Symmetric",
        lad_cols
    ] = 0
    # Save full summary
    df_no_mixtures.to_csv("../results/Summary_LAD_with_symmetrical_no_mixtures.txt", index=False)
    print("Full summary no mixtures has been written to ../results/Summary_LAD_with_symmetrical_no_mixtures.txt")

    # ----------------------------
    # 7) Drivers info
    # ----------------------------

    #divisions_one_driver = divisions.loc[no_semicolon_and_not_empty(divisions["drivers"])].copy()
    df_summary_one_knownDriver = df_no_mixtures.loc[no_semicolon_and_not_empty(df_no_mixtures["knownDrivers"])].copy()

    lad_by_line_all = (
        df_no_mixtures.groupby(["mice_line", "LAD_ML_segm"]).size().reset_index(name="all_samples")
        )
    lad_by_line_one_driver = (
        df_summary_one_knownDriver.groupby(["mice_line", "LAD_ML_segm"]).size().reset_index(name="one_driver")
        )
    lad_by_line = lad_by_line_all.merge(
        lad_by_line_one_driver, on=["mice_line", "LAD_ML_segm"], how="left"
        )
    lad_by_line["one_driver"] = lad_by_line["one_driver"].fillna(0).astype(int)

    lad_by_line.to_csv("../results/LAD_ML_segm_by_mice_line_agg.txt", index=False)
    print(lad_by_line)
    print("LAD by mice line has been written to ../results/LAD_ML_segm_by_mice_line_agg.txt")

    # By driver (one_driver)
    df_summary_one_knownDriver["driver_with_poisiton"] = (
        df_summary_one_knownDriver["drivers"].astype(str) + "_" + df_summary_one_knownDriver["driver_positions_c3h_mutId"].astype(str)
        )

    lad_by_line_by_driver_one_driver = (
        df_summary_one_knownDriver.groupby(["mice_line", "LAD_ML_segm", "driver_with_poisiton"])
        .size()
        .reset_index(name="Freq")
        )
    lad_by_line_by_driver_one_driver.to_csv("../results/LAD_ML_segm_by_mice_line_by_driver_agg.txt", index=False)
    print("LAD by mice line by driver has been written to ../results/LAD_ML_segm_by_mice_line_by_driver_agg.txt")

    # Many drivers summaries
    lad_many_drivers = df_no_mixtures.loc[df_no_mixtures["drivers"].fillna("").str.contains(";", regex=False)].copy()
    lad_many_drivers["comments"] = ""

    lad_many_drivers.loc[
        (lad_many_drivers["LAD_ML_segm"] == "1") & (lad_many_drivers["segments_HMM_multi_state_M3_prop"] == 1),
        "comments"
        ] = "4_cells_survived"

    lad_many_drivers.loc[
        (lad_many_drivers["LAD_ML_segm"] == "1")
        & (lad_many_drivers["segments_HMM_multi_state_M2_prop"] > 0.05)
        & (lad_many_drivers["segments_HMM_multi_state_M3_prop"] > 0.05),
        "comments"
    ] = "3_cells_survived"

    lad_many_drivers.loc[
        (lad_many_drivers["LAD_ML_segm"] == "1") & (lad_many_drivers["segments_HMM_multi_state_M2_prop"] == 1),
        "comments"
    ] = "not_all_survived+different_clone_size???"

    lad_many_drivers.loc[
        (~lad_many_drivers["LAD_ML_segm"].isin(["0", "1"]))
        & (lad_many_drivers["segments_HMM_multi_state_M3_prop"] == 0)
        & (lad_many_drivers["segments_HMM_multi_state_M1_prop"] > 0.05),
        "comments"
    ] = "LAD=1 or LAD=2not_all_survived+different_clone_size??"

    lad_many_drivers_out = lad_many_drivers[["Sample", "mice_line", "drivers", "LAD_ML_segm", "comments"]].copy()

    df_mixtures_many_drivers = df_mixtures.loc[df_mixtures["drivers"].fillna("").str.contains(";", regex=False)].copy()
    df_mixtures_many_drivers["comments"] = np.where(
        df_mixtures_many_drivers["Ploidy_by_HMM_as"] == "Tetra",
        "as_by_HMM_as",
        "as_by_clon-subclon_HMM_as",
    )
    df_mixtures_many_drivers["division"] = "Mixture"
    df_mixtures_many_drivers_out = df_mixtures_many_drivers[["Sample", "mice_line", "drivers", "division", "comments"]].copy()

    df_mixtures_many_drivers = pd.concat([lad_many_drivers_out, df_mixtures_many_drivers_out], ignore_index=True)
    df_mixtures_many_drivers.to_csv("../results/Summary_many_drivers.txt", index=False)
    print("Summary_many_drivers has been written to ../results/Summary_many_drivers.txt")




if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)

