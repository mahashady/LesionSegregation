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
import pysam



""" 
"""
def main(args):
    subset = "enriched"
    phased = 0
    anti = 0
    no_phasing = 0
    nucleotides = ['A','T','G','C']
    chromosomes = [str(i) for i in range(1,23)]
    files = glob.glob("../results/multi_sites_by_sample_"+ subset + "_chemo.alkyl.immuno/*_multi_sites.txt")
    outfile = open("../results/phased_multi_" + subset + "_chemo.alkyl.immuno.txt","w")
    outfile_summary = open("../results/phased_summary_multi_" + subset + "_chemo.alkyl.immuno.txt","w")
    all_multi = 0
    no_germline = 0
    for file_name in files:
        file = open(file_name)
        sample = file_name.split("/")[7].split("_multi")[0]
        print(sample)
        for line in file:
            data = line.strip("\r\n").split("\t")
            #print(data)
            all_multi += 1
            outfile.write(sample)
            outfile.write("\t")
            outfile.write(data[1]+"_"+data[3])
            if data[6] == "NA":
                no_germline += 1
                outfile.write("\tNA")
                no_phasing += 1
                outfile.write("\tNA")
            else:
                chrom = data[1]
                ref_multi = data[2]
                pos_multi = int(data[3])
                alt1_multi = data[4].split(":")[0]
                alt2_multi = data[5].split(":")[0]
                germlines = data[6].split(";")
                n_germlines = len(germlines)
                germline_indels = 0
                germline_snvs = []
                print("Multi pos=", str(pos_multi))
                for germline in range(n_germlines):
                    pos_germline = int(germlines[germline].split(":")[0])
                    print("Germline pos=", str(pos_germline))
                    ref_germline = germlines[germline].split(":")[1]
                    alt_germline = germlines[germline].split(":")[2]
                    if len(ref_germline) != 1 or len(alt_germline) != 1:
                        germline_indels += 1
                    else:
                        try:
                            result = phase(sample, chrom, pos_multi, ref_multi, alt1_multi, alt2_multi, pos_germline, ref_germline, alt_germline, subset)
                            germline_snvs.append(result)
                            print(result)
                        except:
                            print("Some problems with bam files")
                if germline_snvs == []:
                    outfile.write("\tonly_germline_indels")
                    no_phasing += 1
                    outfile.write("\tNA")

                else:
                    outfile.write("\t"+str(germline_snvs))
                    germline1 = germline_snvs[0]
                    if (germline1[2] != 0 and germline1[3] == 0 and germline1[4] != 0 and germline1[5] == 0):
                        phased += 1
                        outfile.write("\tPhased")
                    elif (germline1[2] == 0 and germline1[3] != 0 and germline1[4] == 0 and germline1[5] != 0):
                        phased += 1
                        outfile.write("\tPhased")
                    elif (germline1[2] == 0 and germline1[3] != 0 and germline1[4] != 0 and germline1[5] == 0):
                        anti += 1
                        outfile.write("\tAnti")
                    elif (germline1[2] != 0 and germline1[3] == 0 and germline1[4] == 0 and germline1[5] != 0):
                        anti += 1
                        outfile.write("\tAnti")
                    else:
                        no_phasing += 1
                        outfile.write("\tNo_phasing")

            outfile.write("\n")
    print(all_multi)
    print(no_germline)
    outfile_summary.write("Phased\t" + str(phased) + "\n" + "Anti\t" + str(anti) + "\n" + "No phasing\t" + str(no_phasing))



        

def phase(sample, chrom, pos_multi, ref_multi, alt1_multi, alt2_multi, pos_germline, ref_germline, alt_germline, subset):
    if subset == "enriched":
        samfile = pysam.AlignmentFile("/workspace/datasets/hartwig/20230914/scripts/minibam/results/" + sample + "T_dedup.realigned.minisorted.bam", "rb")
        output1 = pysam.AlignmentFile("/workspace/datasets/hartwig/20230914/scripts/minibam/results/" + sample + "T_" + chrom + ":"+str(pos_multi) + "-" + str(pos_germline) + "_dedup.realigned.mini.intersect1.bam", "wb", template=samfile)
        output2 = pysam.AlignmentFile("/workspace/datasets/hartwig/20230914/scripts/minibam/results/" + sample + "T_" + chrom + ":" + str(pos_multi) + "-" + str(pos_germline) + "_dedup.realigned.mini.intersect2.bam", "wb", template=samfile)
    else:
        samfile = pysam.AlignmentFile("/workspace/datasets/hartwig/20230914/scripts/minibam/results/batch2/" + sample + "T_dedup.realigned.minisorted.bam", "rb")
        output1 = pysam.AlignmentFile("/workspace/datasets/hartwig/20230914/scripts/minibam/results/batch2/" + sample + "T_" + chrom + ":"+str(pos_multi) + "-" + str(pos_germline) + "_dedup.realigned.mini.intersect1.bam", "wb", template=samfile)
        output2 = pysam.AlignmentFile("/workspace/datasets/hartwig/20230914/scripts/minibam/results/batch2/" + sample + "T_" + chrom + ":" + str(pos_multi) + "-" + str(pos_germline) + "_dedup.realigned.mini.intersect2.bam", "wb", template=samfile)
    for read in samfile.fetch(chrom, pos_multi - 1, pos_multi):
        output1.write(read)
    output1.close()
    if subset == "enriched":
        file2index =  "/workspace/datasets/hartwig/20230914/scripts/minibam/results/" + sample + "T_" + chrom + ":"+str(pos_multi) + "-" + str(pos_germline) + "_dedup.realigned.mini.intersect1.bam"
    else:
        file2index =  "/workspace/datasets/hartwig/20230914/scripts/minibam/results/batch2/" + sample + "T_" + chrom + ":"+str(pos_multi) + "-" + str(pos_germline) + "_dedup.realigned.mini.intersect1.bam"
    
    pysam.index(file2index)

    if subset == "enriched":    
        output1 = pysam.AlignmentFile("/workspace/datasets/hartwig/20230914/scripts/minibam/results/" + sample + "T_" + chrom + ":"+str(pos_multi) + "-" + str(pos_germline) + "_dedup.realigned.mini.intersect1.bam")
    else:
        output1 = pysam.AlignmentFile("/workspace/datasets/hartwig/20230914/scripts/minibam/results/batch2/" + sample + "T_" + chrom + ":"+str(pos_multi) + "-" + str(pos_germline) + "_dedup.realigned.mini.intersect1.bam")

    for read in output1.fetch(chrom, pos_germline - 1, pos_germline):
        #print(read)
        output2.write(read)

    samfile.close()
    output2.close()
    if subset == "enriched":
        file2index2 = "/workspace/datasets/hartwig/20230914/scripts/minibam/results/" + sample + "T_" + chrom + ":" + str(pos_multi) + "-" + str(pos_germline) + "_dedup.realigned.mini.intersect2.bam"
    else:
        file2index2 = "/workspace/datasets/hartwig/20230914/scripts/minibam/results/batch2/" + sample + "T_" + chrom + ":" + str(pos_multi) + "-" + str(pos_germline) + "_dedup.realigned.mini.intersect2.bam"

    pysam.index(file2index2)

    if subset == "enriched":
        samfile = pysam.AlignmentFile("/workspace/datasets/hartwig/20230914/scripts/minibam/results/" + sample + "T_" + chrom + ":" + str(pos_multi) + "-" + str(pos_germline) + "_dedup.realigned.mini.intersect2.bam", "rb")
    else:
        samfile = pysam.AlignmentFile("/workspace/datasets/hartwig/20230914/scripts/minibam/results/batch2/" + sample + "T_" + chrom + ":" + str(pos_multi) + "-" + str(pos_germline) + "_dedup.realigned.mini.intersect2.bam", "rb")

    dict = {}
    if pos_multi < pos_germline:
        start = pos_multi - 1
        end = pos_germline - 1
    else:
        end = pos_multi - 1
        start = pos_germline - 1

    for pileupcolumn in samfile.pileup(chrom, start, end):
        if pileupcolumn.pos == (pos_multi - 1) or pileupcolumn.pos == (pos_germline - 1):
            print(pileupcolumn.n)
            for pileupread in pileupcolumn.pileups:
                if pileupread.alignment.query_name not in dict:
                    dict[pileupread.alignment.query_name] = ["N","N"]
                if pileupcolumn.pos == (pos_multi - 1):
                    dict[pileupread.alignment.query_name][0] = pileupread.alignment.query_sequence[pileupread.query_position]
                if pileupcolumn.pos == (pos_germline - 1):
                    dict[pileupread.alignment.query_name][1] = pileupread.alignment.query_sequence[pileupread.query_position]

    dict_filtered = {}
    for el in dict:
        if dict[el][0] != "N" and dict[el][1] != "N":
            dict_filtered[el] = dict[el]

    print("Len dict = " + str(len(dict)))
    print("Len dict = " + str(len(dict_filtered)))
    ref_list = []
    alt1_list = []
    alt2_list = []
    for i in dict_filtered:
        if dict_filtered[i][0] == alt1_multi:
            alt1_list.append(dict_filtered[i][1])
        elif dict_filtered[i][0] == alt2_multi:
            alt2_list.append(dict_filtered[i][1])
        elif dict_filtered[i][0] == ref_multi:    
            ref_list.append(dict_filtered[i][1])
    n_ref_ref_list = ref_list.count(ref_germline)
    n_alt_ref_list = ref_list.count(alt_germline)
    n_ref_alt1_list = alt1_list.count(ref_germline)
    n_alt_alt1_list = alt1_list.count(alt_germline)
    n_ref_alt2_list = alt2_list.count(ref_germline)
    n_alt_alt2_list = alt2_list.count(alt_germline)
    return([n_ref_ref_list, n_alt_ref_list, n_ref_alt1_list, n_alt_alt1_list, n_ref_alt2_list, n_alt_alt2_list])

        


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)


			
















