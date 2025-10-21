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


def main(args):
    with open("/workspace/projects/lesion_segregation/metastatic_tumors/results/samples_lists.Platinum.SBS_less10.json", "r") as fl:
# with open("/workspace/projects/lesion_segregation/metastatic_tumors/results/samples_lists.Alkylating.json", "r") as fl:
        list_of_samples = json.load(fl)[setting]
    df_muts_stat = pd.read_table("/workspace/projects/lesion_segregation/metastatic_tumors/results/all_"+ setting  + ".chemo_alkyl_immuno.txt")
    print(df_muts_stat.shape)
    print(sorted(df_muts_stat["n_multi"]))

    df_muts_stat = df_muts_stat[df_muts_stat["patientIdentifier"].isin(list_of_samples)]
    print(sorted(df_muts_stat["n_multi"]))
    print(np.median(df_muts_stat["n_multi"]))

    df_signatures = pd.read_table("/workspace/projects/mutfootprints/mutfootprints/data/HMF/20190502/signature_extraction/results/SignatureAnalyzer/snvs/exposures/Pan_full/Pan_full.exposures.tsv",index_col=0)
    print(df_signatures.head)
    print(df_signatures.shape)
    list_of_samples = [a + 'T' for a in list_of_samples]
    df_signatures = df_signatures[df_signatures.columns.intersection(list_of_samples)]
    print(df_signatures.shape)
    print(df_signatures.index.tolist())
    print(np.mean(df_signatures[df_signatures.index == "21_SBS31_0.953955_1"]))

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)