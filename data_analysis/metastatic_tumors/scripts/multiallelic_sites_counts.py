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
from collections import Counter


""" 
"""
nucleotides = ['A','T','G','C']
chromosomes = [str(i) for i in range(1,23)]

def make_empty_dict():
    dict = {}
    for nucleotide in nucleotides:
        dict[nucleotide] = {}
        for chromosome in chromosomes:
            dict[nucleotide][chromosome] = []
    return dict


def main(args):
    i = 0
    # HERE the PATH CONTAINING SOMATIC VCFs SHOULD BE UPDATED
    files = glob.glob(".../hartwig/20230914/somatic/*/purple/*.purple.somatic.vcf.gz")
    outfile_full = open("../results/ALL_Hartwig_bi_multi.by_chrom.txt", "a")
    outfile_full.write("sample\tchr\tA>N\tT>N\tG>N\tC>N\tA>N_triallelic\tT>N_triallilic\tG>N_triallelic\tC>N_triallelic\tA>N_cuatro\tT>N_cuatro\tG>N_cuatro\tC>N_cuatro\n")
    outfile_short = open("../results/ALL_Hartwig_bi_multi.txt", "a")
    outfile_short.write("sample\tn_biallelic\tn_multi\n")
    for file in files:
        sample = file.split("/")[6]
        i +=1
        print(i)
        print(sample)
        dict = make_empty_dict()
        f=gzip.open(file,'rb')
        file_content=f.readlines()
        for line in file_content:
            line = line.decode("utf-8")
            if line[0] != "#":
                data = line.strip("\r\n").split("\t")
                chrom, pos, ref, alt = data[0], data[1], data[3], data[4]
                if chrom in chromosomes and ref in nucleotides and alt in nucleotides:
                    dict[ref][chrom].append(pos)
        
        full_result = [0,0]
        for chromosome in chromosomes:
            result = {'A':[],'T':[],'G':[],'C':[]}
            for nucleotide in nucleotides:
                n_bi_multi = Counter(Counter(dict[nucleotide][chromosome]).values())
                n_biallelic = n_bi_multi[1]
                n_trillelic = n_bi_multi[2] 
                n_cuatro = n_bi_multi[3]
                result[nucleotide] = [n_biallelic, n_trillelic, n_quatro]
            outfile_full.write(sample + "\t" + chromosome + "\t" + str(result['A'][0])+ "\t" + str(result['T'][0])+ "\t" + str(result['G'][0])+ "\t" + str(result['C'][0])+ "\t" + str(result['A'][1])+ "\t" + str(result['T'][1])+ "\t" + str(result['G'][1])+ "\t" + str(result['C'][1])+ "\t" + str(result['A'][2])+ "\t" + str(result['T'][2])+ "\t" + str(result['G'][2])+ "\t" + str(result['C'][2]) + "\n")
            all_bi = result['A'][0] + result['T'][0] + result['G'][0] + result['C'][0]
            all_multi = result['A'][1] + result['T'][1] + result['G'][1] + result['C'][1] + result['A'][2] + result['T'][2] + result['G'][2] + result['C'][2]
            full_result[0] += all_bi
            full_result[1] += all_multi
        outfile_short.write(sample + "\t" + str(full_result[0]) + "\t" + str(full_result[1]) + "\n")
        


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)


			
















