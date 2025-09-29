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

def load_germline(file):
    germline_dict = make_empty_dict_chroms()
    if os.path.isfile(file):
        f=gzip.open(file,'rb')
        file_content=f.readlines()
        for line in file_content:
            line = line.decode("utf-8")
            if line[0] != "#":
                data = line.strip("\r\n").split("\t")
                chrom, pos, ref, alt, gt_n = data[0], int(data[1]), data[3], data[4], data[9]                
                if chrom in chromosomes and gt_n.split(":")[0] == "0/1":
                    germline_dict[chrom][pos] = ref + ":" + alt + ":" + gt_n
    else:
        print("Germline file is absent")
    return(germline_dict)



def main(args):
    subset = "NONenriched"
    summary_outfile = open("../results/summary_multi_sites_with_germline_" + subset + "_chemo.alkyl.immuno.txt", "a")
    i = 0
    n_multi = 0
    n_multi_with_germline = 0
    multi_samples_file = open("../results/all_" + subset + ".chemo_alkyl_immuno.txt")
    multi_samples_list = []
    for line in multi_samples_file:
        multi_sample = line.strip("\r\n").split(",")[0].strip('\"')
        n_multiallelic = line.strip("\r\n").split(",")[3]
        if multi_sample not in multi_samples_list and multi_sample!="patientIdentifier":
            if int(n_multiallelic) > 0:
                multi_samples_list.append(multi_sample)
    
    for sample in multi_samples_list:
        print(sample)
        file_name = glob.glob("/workspace/datasets/hartwig/20230914/somatic/"+sample+"*/purple/" + sample + "*.purple.somatic.vcf.gz")[0]
        if os.path.isfile("../results/multi_sites_by_sample_" + subset + ".chemo_alkyl_immuno/" + sample + "_multi_sites.txt") == False: 
            print("Creating file")
            outfile = open("../results/multi_sites_by_sample_" + subset + ".chemo_alkyl_immuno/" + sample + "_multi_sites.txt", "a")
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
            if subset == "enriched":
                germline_file_name = glob.glob("/workspace/datasets/hartwig/20230914/germline/*" + sample +".annotated.vcf.gz")
            else:
                germline_file_name = glob.glob("/workspace/datasets/hartwig/20230914/germline/2nd_batch/*" + sample +".annotated.vcf.gz")

            if len(germline_file_name) > 0:
                germline_dict = load_germline(germline_file_name[0])
            else:
                print(sample)
                print("Germline file is absent!")
                germline_dict = make_empty_dict_chroms()

            for chromosome in chromosomes:
                result = {'A':[],'T':[],'G':[],'C':[]}
                for nucleotide in nucleotides:
                    for pos in dict[nucleotide][chromosome]:
                        if len(dict[nucleotide][chromosome][pos]) > 1:
                            n_multi += 1
                            germline_list = []
                            for i in germline_dict[chromosome]:
                                if i >= pos-150 and i <= pos + 150:
                                    germline_list.append(i)
                            
                            context = hg19(chromosome, int(pos)-1, 3)
                            outfile.write(sample + "\t" + chromosome + "\t" + nucleotide + "\t" + str(pos)+ "\t" + context)
                            for mut in dict[nucleotide][chromosome][pos]:
                                outfile.write("\t" + mut)
                            #print(germline_list)
                            if len(germline_list) > 0:
                                n_multi_with_germline += 1
                                summary_outfile.write(sample + "\t" + chromosome + ":" + str(pos-200) + "-" + str(pos+200) + "\n")
                                n_germline = 1
                                for germline in germline_list:
                                    if n_germline == 1:
                                        outfile.write("\t")
                                    else:
                                        outfile.write(";")
                                    outfile.write(str(germline) + ":" + germline_dict[chromosome][germline])
                                    n_germline +=1
                            else:
                                outfile.write("\tNA")

                            outfile.write("\n")
    print("All multi sites = " + str(n_multi))
    print("Multi sites with germline = " + str(n_multi_with_germline))
        

  

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)
