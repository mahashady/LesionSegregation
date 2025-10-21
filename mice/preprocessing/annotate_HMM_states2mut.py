#$ -S /usr/bin/python
import argparse
import glob
import string
import gzip
import re
import os
import subprocess
import pandas as pd
import statistics
""" 
"""
def merge_HMMas(sample_name):
    df_HMM_input = pd.read_csv("../LAD/results/HMM/input_HMM/" + sample_name + ".hmm", sep=" ", header=None)
    df_HMM_input.columns = ['Chr', 'Pos', 'Change', 'Context', 'Nref', 'Nalt', 'Sample','Multi_class','As_class']
    df_HMM_multi = pd.read_csv("../LAD/results/HMM/output_HMM/" + sample_name + ".hmm.PloidyHMM", header=None)
    df_HMM_multi.columns = ['HMM_multi_state']
    print(df_HMM_input.head())
    print(df_HMM_multi.head())
    df_HMM_input['HMM_multi_state'] = df_HMM_multi["HMM_multi_state"]
    df_HMM_input["mut_id"] = df_HMM_input["Chr"].astype("str") + ":" + df_HMM_input["Pos"].astype("str")
    print(df_HMM_input.head())
    print(df_HMM_multi.head())
    return df_HMM_input


def main(args):
    df_MRCA = pd.read_csv("../LAD/results/Summary_divisions_with_symmetrical_no_mixtures.txt", sep=",")
    print(df_MRCA.head())
    samples = df_MRCA["sample"]
    for sample_name in samples:
        print(sample_name)
        input_file = "../data/mutations_vs_genes/" + sample_name + ".with_gene_annot.nodMat"
        if os.path.isfile(input_file):
            df_HMM_multi = merge_HMMas(sample_name)
            path2out ="../data/mutations_vs_genes_vs_HMM_multi_state/"    
            df_input = pd.read_csv(input_file, sep=",")
            print(df_input.head())
            df_input["mut_id"] = df_input["chr"].astype("str") + ":" + df_input["pos"].astype("str")
            print(df_input.head())
            df_result = pd.merge(df_input, df_HMM_multi, on="mut_id")
            print(df_result.head())
            df_result.to_csv(path2out + sample_name + ".gene.HMMmulti_state.nodMat", index=False, sep= ",") 

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)