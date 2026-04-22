#$ -S /usr/bin/python
import argparse
import pandas as pd
import json
import csv
from collections import defaultdict
import itertools

# Setting could be enriched or NONenriched
# Treatment could be Platinum.SBS_more10 or Alkylating
setting = "NONenriched"
treatment = "Alkylating"
# Load genome composition from a JSON file
with open("../data/genome_counts_tribases.json", "r") as f:
    genome_context = json.load(f)  

# Load lists of samples
with open("../results/samples_lists." + treatment + ".json", "r") as fl:
    lists_of_samples = json.load(fl)  





def main(args):
    samples_dir = "../results/bi_spectrum_by_sample_" + setting + ".chemo_alkyl_immuno/"
    total_expected = defaultdict(float)
    list_of_samples = lists_of_samples[setting]
    for sample in list_of_samples:
        filepath = samples_dir + sample + "_bi_spectrum.txt"
        mutation_rates = load_sample_mutation_rates(filepath)
        sample_expected = compute_expected_triallelic(mutation_rates, genome_context)
        for label, count in sample_expected.items():
            total_expected[label] += count

    with open("../results/expected_triallelic_spectrum/" + setting + "." + treatment + ".expected_triallelic.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["context", "expected_count"])
        for label in sorted(total_expected.keys()):
            writer.writerow([label, total_expected[label]])

    print("Cohort-level triallelic spectrum saved to 'cohort_expected_triallelic.csv'")



if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)