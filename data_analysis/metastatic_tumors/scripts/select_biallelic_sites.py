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
from bgreference import hg19
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
            dict[nucleotide][chromosome] = {}
    return dict

def make_empty_dict_chroms():
    dict = {}
    for chromosome in chromosomes:
        dict[chromosome] = {}
    return dict

def main(args):
    subset = "NONenriched"
    i = 0
    multi_samples_file = open("../results/all_" + subset + ".chemo_alkyl_immuno.txt")
    multi_samples_list = []
    for line in multi_samples_file:
        multi_sample = line.strip("\r\n").split("\t")[0].strip('\"')
        n_multiallelic = line.strip("\r\n").split("\t")[3]
        if multi_sample not in multi_samples_list and multi_sample!="patientIdentifier":
            if int(n_multiallelic) > 0:
                multi_samples_list.append(multi_sample)
    
    for sample in multi_samples_list:
        print(sample)
        # HERE the PATH CONTAINING SOMATIC VCFs SHOULD BE UPDATED
        file_name = glob.glob("./hartwig/20230914/somatic/"+sample+"*/purple/" + sample + "*.purple.somatic.vcf.gz")[0]
        if os.path.isfile("../results/bi_sites_by_sample_" + subset + ".chemo_alkyl_immuno/" + sample + "_bi_sites.txt") == False: 
            print("Creating file")
            outfile = open("../results/bi_sites_by_sample_" + subset + ".chemo_alkyl_immuno/" + sample + "_bi_sites.txt", "a")
            i +=1
            #print(i)
            #print(sample)
            dict = make_empty_dict()
            f=gzip.open(file_name,'rb')
            file_content=f.readlines()
            for line in file_content:
                line = line.decode("utf-8")
                if line[0] != "#":
                    data = line.strip("\r\n").split("\t")
                    chrom, pos, ref, alt = data[0], int(data[1]), data[3], data[4]
                    tumour_ad, tumour_af = data[10].split(":")[1], data[10].split(":")[2]
                    mutation = alt + ":" + tumour_ad + ":" + tumour_af
                    if chrom in chromosomes and ref in nucleotides and alt in nucleotides:
                        if pos not in dict[ref][chrom]:
                            dict[ref][chrom][pos]=[]
                        dict[ref][chrom][pos].append(mutation)
            
            print("Dict created")
            for chromosome in chromosomes:
                result = {'A':[],'T':[],'G':[],'C':[]}
                for nucleotide in nucleotides:
                    for pos in dict[nucleotide][chromosome]:
                        if len(dict[nucleotide][chromosome][pos]) == 1:
                            context = hg19(chromosome, int(pos)-1, 3)
                            outfile.write(sample + "\t" + chromosome + "\t" + nucleotide + "\t" + str(pos)+ "\t" + context + "\t" + dict[nucleotide][chromosome][pos][0] + "\n")
    

  

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)
