#$ -S /usr/bin/python
import argparse
import glob
import gzip
import os
import subprocess
import pandas as pd
from bgreference import hg19
from utils.utils import read_vcf_to_mut_dict



""" 
"""
nucleotides = ['A','T','G','C']
chromosomes = [str(i) for i in range(1,23)]


def make_empty_dict():
    dict = {}
    for nucleotide in nucleotides:
        dict[nucleotide] = {}
        for chromosome in chromosomes:
            dict[nucleotide][chromosome] = {}
    return dict

def make_empty_dict_chroms():
    dict = {}
    for chromosome in chromosomes:
        dict[chromosome] = {}
    return dict

def write_single_alt_sites(outfile, sample, mut_dict, chromosomes, nucleotides):
    for chromosome in chromosomes:
        for nucleotide in nucleotides:
            for pos, mutations in sorted(mut_dict[nucleotide][chromosome].items()):
                if len(mutations) == 1:
                    context = hg19(chromosome, pos - 1, 3)
                    outfile.write(
                        f"{sample}\t{chromosome}\t{nucleotide}\t{pos}\t{context}\t{mutations[0]}\n"
                    )


def main(args):
    subset = "enriched"

    #create output folder
    output_folder = "../results/bi-allelic_sites/bi_sites_by_sample_" + subset + ".chemo_alkyl_immuno/"
    os.makedirs(output_folder, exist_ok=True)

    i = 0
    all_subset_file = open("../results/enrichment/all_" + subset + ".chemo_alkyl_immuno.txt")
    df_subset = pd.read_csv(all_subset_file, sep = "\t")
    multi_samples_list = get_multi_sample_list(df_subset, subset)
    for sample in multi_samples_list:
        print(sample)
        # PATH CONTAINING SOMATIC VCFs 
        file_name = glob.glob("/data/bbg/datasets/hartwig/20230914/somatic/" + sample + "T/purple/" + sample + "*.purple.somatic.vcf.gz")[0]
        outfile = open("../results/bi-allelic_sites/bi_sites_by_sample_" + subset + ".chemo_alkyl_immuno/" + sample + "_bi_sites.txt", "a")
        mut_dict = read_vcf_to_mut_dict(file_name, sample_name = sample+"T")
        print("Sample dictionary created")
        write_single_alt_sites(outfile, sample, mut_dict, chromosomes, nucleotides)
    

  

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)
