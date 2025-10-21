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
from compute_expected_triallelic_spectrum import compute_expected_triallelic


# Load genome composition from a JSON file
with open("../data/genome_counts_tribases.json", "r") as f:
    genome_context = json.load(f)  

def main(args):
    # treatmnet could be Platinum or Alkylating
    treatment = "Platinum"
    df_signatures = pd.read_csv("../data//Pan_full.processes.tsv", sep="\t")
    df_signatures_experimental = pd.read_csv("../data/human_sbs96_unfiltered_v1_0.txt", sep="\t")
    print(df_signatures.head())
    print(df_signatures_experimental.head())
    df_signatures = df_signatures.merge(
        df_signatures_experimental,
        left_on="Unnamed: 0",
        right_on="MutationType"
    )
    print(df_signatures.head())
    df_signatures["Unnamed: 0"] = df_signatures["Unnamed: 0"].astype(str)
    df_signatures["context"] = df_signatures["Unnamed: 0"].apply(
    lambda x: x[0] + x[2] + x[6] + ">" + x[4]
    )
    if treatment == "Platinum":
        signatures_of_interest = ["25_1", "37_1", "21_SBS31_0.953955_1", "14_1","31_SBS17b_0.968799_1"]
    elif treatment == "Alkylating":
        signatures_of_interest = ["cyclophosphamide_557117b73fe2", "38_SBS2_0.996907_1", "20_SBS13_0.948838_1"]
    else:
        print("Bad treatment type")
    df_signatures = df_signatures[["context"] + signatures_of_interest]
    print(df_signatures.head())

    signatures_expected_list = []
    for signature in signatures_of_interest:
        mutation_rates = dict(zip(df_signatures["context"], df_signatures[signature]))
        signature_expected = compute_expected_triallelic(mutation_rates, genome_context)
        signatures_expected_list.append(signature_expected)
        print(signature_expected)
    print(len(signatures_expected_list))

    df = pd.DataFrame(signatures_expected_list).T
    df.columns = signatures_of_interest
    print(df.head())

    outfile = "../results/expected_triallelic_spectrum/signatures.expected_triallelic." + treatment  + ".csv"
    df.to_csv(outfile, index=True)

    print("Expected triallelic spectrum of platinum signatures saved to 'signatures.expected_triallelic.csv'")


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)