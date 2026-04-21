#!/usr/bin/env python3

import argparse
import csv
import json
import os
import itertools
import pandas as pd
from collections import defaultdict
from utils.expected_spectra_utils import *

GENOME_CONTEXT_DEFAULT = "../data/genome_counts_tribases.json"
GENOME_SIGNATURES_DEFAULT = "../data/Pan_full.processes.tsv"
SAMPLE_SUFFIX_DEFAULT = "_bi_spectrum.txt"

SIGNATURES_BY_TREATMENT = {
    "Platinum": [
        "25_1",
        "37_1",
        "21_SBS31_0.953955_1",
        "14_1",
        "31_SBS17b_0.968799_1",
    ],
    "Alkylating": [
        "cyclophosphamide_557117b73fe2",
        "38_SBS2_0.996907_1",
        "20_SBS13_0.948838_1",
    ],
}

SAMPLE_COHORTS = {
    "alkylating_enriched": {
        "sample_list": "../results/samples_lists/samples_lists.Alkylating.json",
        "setting_key": "enriched",
        "samples_dir": "../results/bi-allelic_sites/spectra/bi_spectrum_by_sample_enriched.chemo_alkyl_immuno/",
        "output": "../results/multi-allelic_sites/spectra/expected/enriched.Alkylating.expected_triallelic.csv",
    },
    "alkylating_NONenriched": {
        "sample_list": "../results/samples_lists/samples_lists.Alkylating.json",
        "setting_key": "NONenriched",
        "samples_dir": "../results/bi-allelic_sites/spectra/bi_spectrum_by_sample_NONenriched.chemo_alkyl_immuno/",
        "output": "../results/multi-allelic_sites/spectra/expected/NONenriched.Alkylating.expected_triallelic.csv",
    },
    "platinum_enriched_low_SBS17b": {
        "sample_list": "../results/samples_lists/samples_lists.Platinum.sigfilt_SBS17b_prop_lt_0.1.json",
        "setting_key": "enriched",
        "samples_dir": "../results/bi-allelic_sites/spectra/bi_spectrum_by_sample_enriched.chemo_alkyl_immuno/",
        "output": "../results/multi-allelic_sites/spectra/expected/enriched.Platinum.sigfilt_SBS17b_prop_lt_0.1.expected_triallelic.csv",
    },
    "platinum_NONenriched_low_SBS17b": {
        "sample_list": "../results/samples_lists/samples_lists.Platinum.sigfilt_SBS17b_prop_lt_0.1.json",
        "setting_key": "NONenriched",
        "samples_dir": "../results/bi-allelic_sites/spectra/bi_spectrum_by_sample_NONenriched.chemo_alkyl_immuno/",
        "output": "../results/multi-allelic_sites/spectra/expected/NONenriched.Platinum.sigfilt_SBS17b_prop_lt_0.1.expected_triallelic.csv",
    },
    "platinum_enriched_high_SBS17b": {
        "sample_list": "../results/samples_lists/samples_lists.Platinum.sigfilt_SBS17b_prop_gt_0.1.json",
        "setting_key": "enriched",
        "samples_dir": "../results/bi-allelic_sites/spectra/bi_spectrum_by_sample_enriched.chemo_alkyl_immuno/",
        "output": "../results/multi-allelic_sites/spectra/expected/enriched.Platinum.sigfilt_SBS17b_prop_gt_0.1.expected_triallelic.csv",
    },
    "platinum_NONenriched_high_SBS17b": {
        "sample_list": "../results/samples_lists/samples_lists.Platinum.sigfilt_SBS17b_prop_gt_0.1.json",
        "setting_key": "NONenriched",
        "samples_dir": "../results/bi-allelic_sites/spectra/bi_spectrum_by_sample_NONenriched.chemo_alkyl_immuno/",
        "output": "../results/multi-allelic_sites/spectra/expected/NONenriched.Platinum.sigfilt_SBS17b_prop_gt_0.1.expected_triallelic.csv",
    },
}

def run_one_sample_cohort(cohort_name, genome_context, sample_suffix):
    cfg = SAMPLE_COHORTS[cohort_name]
    sample_lists = load_json(cfg["sample_list"])
    samples = sample_lists[cfg["setting_key"]]

    total_expected = defaultdict(float)

    for sample in samples:
        filepath = os.path.join(cfg["samples_dir"], f"{sample}{sample_suffix}")
        mutation_rates = load_sample_mutation_rates(filepath)
        sample_expected = compute_expected_triallelic(mutation_rates, genome_context)
        #print(sample_expected)

        for label, count in sample_expected.items():
            total_expected[label] += count

    ensure_parent_dir(cfg["output"])
    with open(cfg["output"], "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["context", "expected_count"])
        for label in sorted(total_expected):
            writer.writerow([label, total_expected[label]])

    print(f"[OK] {cohort_name} -> {cfg['output']}")

def run_samples(args):
    genome_context = load_json(args.genome_context)

    if args.run_all:
        cohort_names = list(SAMPLE_COHORTS.keys())
    else:
        cohort_names = args.cohort

    for cohort_name in cohort_names:
        if cohort_name not in SAMPLE_COHORTS:
            raise ValueError(f"Unknown cohort: {cohort_name}")
        run_one_sample_cohort(
            cohort_name=cohort_name,
            genome_context=genome_context,
            sample_suffix=args.sample_suffix,
        )


def run_signatures(args):
    genome_context = load_json(args.genome_context)

    df_signatures = pd.read_csv(args.signatures_file, sep="\t")
    df_experimental = pd.read_csv(args.experimental_file, sep="\t")

    df = df_signatures.merge(
        df_experimental,
        left_on="Unnamed: 0",
        right_on="MutationType"
    )

    df["Unnamed: 0"] = df["Unnamed: 0"].astype(str)
    df["context"] = df["Unnamed: 0"].apply(parse_signature_context)

    signatures_of_interest = SIGNATURES_BY_TREATMENT[args.treatment]

    results = {}
    for signature in signatures_of_interest:
        mutation_rates = dict(zip(df["context"], df[signature]))
        results[signature] = compute_expected_triallelic(mutation_rates, genome_context)

    out_df = pd.DataFrame(results).fillna(0.0)
    out_df.index.name = "context"

    output = "../results/multi-allelic_sites/spectra/expected/Signatures." + args.treatment + ".expected_triallelic.csv"
    ensure_parent_dir(output)
    out_df.to_csv(output)

    print(f"[OK] signatures {args.treatment} -> {output}")


def build_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_samples = subparsers.add_parser("samples")
    p_samples.add_argument(
        "--genome-context",
        default=GENOME_CONTEXT_DEFAULT,
    )
    p_samples.add_argument(
        "--sample-suffix",
        default=SAMPLE_SUFFIX_DEFAULT,
    )
    p_samples.add_argument(
        "--run-all",
        action="store_true",
        help="Run all predefined sample cohorts.",
    )
    p_samples.add_argument(
        "--cohort",
        nargs="+",
        choices=list(SAMPLE_COHORTS.keys()),
        help="Run one or more predefined sample cohorts.",
    )
    p_samples.set_defaults(func=run_samples)

    p_signatures = subparsers.add_parser("signatures")
    p_signatures.add_argument(
        "--treatment",
        required=True,
        choices=list(SIGNATURES_BY_TREATMENT.keys()),
    )
    p_signatures.add_argument(
        "--genome-context",
        default=GENOME_CONTEXT_DEFAULT,
    )
    p_signatures.add_argument(
        "--signatures-file",
        default=GENOME_SIGNATURES_DEFAULT,
    )
    p_signatures.add_argument(
        "--experimental-file",
        default="../data/human_sbs96_unfiltered_v1_0.txt",
    )

    p_signatures.set_defaults(func=run_signatures)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "samples" and not args.run_all and not args.cohort:
        parser.error("For 'samples', provide either --run-all or --cohort.")

    args.func(args)


if __name__ == "__main__":
    main()