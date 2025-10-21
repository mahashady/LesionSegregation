#$ -S /usr/bin/python
import argparse
import glob
import string
import gzip
import re
import os
import subprocess
import pandas as pd
import numpy as np
import random
import statistics
import json
import csv
from collections import defaultdict
import itertools


setting = "NONenriched"

# Load genome composition from a JSON file
with open("/workspace/projects/bladder_ts/results/cord_blood_tws/genome_counts_tribases.json", "r") as f:
    genome_context = json.load(f)  

# Load lists of samples
# with open("/workspace/projects/lesion_segregation/metastatic_tumors/results/samples_lists.Platinum.SBS_more10.json", "r") as fl:
with open("/workspace/projects/lesion_segregation/metastatic_tumors/results/samples_lists.Alkylating.json", "r") as fl:
    lists_of_samples = json.load(fl)  

def load_sample_mutation_rates(filepath):
    mutation_rates = {}
    with open(filepath, "r") as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            if row[0] == "context":
                continue  # Skip header or malformed rows
            context, rate = row[1].strip(), float(row[6])
            mutation_rates[context] = rate
    return mutation_rates

def group_by_trinuc(mutation_rates):
    grouped = defaultdict(list)
    for label, rate in mutation_rates.items():
        trinuc = label[0:3]
        grouped[trinuc].append((label, rate))
    return grouped

def compute_expected_triallelic(mutation_rates, genome_context):
    grouped = group_by_trinuc(mutation_rates)
    print(grouped)
    expected = defaultdict(float)

    for trinuc, mutations in grouped.items():
        context_count = genome_context.get(trinuc)
        if context_count is None:
            continue

        for (label1, rate1), (label2, rate2) in itertools.combinations(mutations, 2):
            joint_prob = rate1 * rate2
            expected_count = context_count * joint_prob
            expected[label1] += expected_count
            expected[label2] += expected_count

    return expected

def main(args):
    samples_dir = "/workspace/projects/lesion_segregation/metastatic_tumors/results/bi_spectrum_by_sample_" + setting + ".chemo_alkyl_immuno/"
    total_expected = defaultdict(float)
    list_of_samples = lists_of_samples[setting]
    for sample in list_of_samples:
        filepath = samples_dir + sample + "_bi_spectrum.txt"
        mutation_rates = load_sample_mutation_rates(filepath)
        sample_expected = compute_expected_triallelic(mutation_rates, genome_context)
        for label, count in sample_expected.items():
            total_expected[label] += count

    # with open("/workspace/projects/lesion_segregation/metastatic_tumors/results/expected_triallelic_spectrum/" + setting + ".Platinum.SBS_more10.expected_triallelic.csv", "w", newline="") as f:
    with open("/workspace/projects/lesion_segregation/metastatic_tumors/results/expected_triallelic_spectrum/" + setting + ".Alkylating.expected_triallelic.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["context", "expected_count"])
        for label in sorted(total_expected.keys()):
            writer.writerow([label, total_expected[label]])

    print("✅ Cohort-level triallelic spectrum saved to 'cohort_expected_triallelic.csv'")



if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)