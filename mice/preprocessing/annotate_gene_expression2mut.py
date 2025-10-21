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
def subset_genes(subset_file_name):
    genes_list = []
    file_subset = open(subset_file_name)
    for line in file_subset:
        line = line.strip("\r\n").split("\t")
        genes_list.append(line[0])
    return genes_list




def main(args):
    genes_MAPK = subset_genes("../../data/gene_sets/MAPK_pathway_genes.txt")
    genes_census = subset_genes("../../data/gene_sets/Census_tier1.txt")
    genes_file = open("../../data/expression/all_expression.csv")
    genes_dict = {}
    for line in genes_file:
        line = line.strip("\r\n").split(";")
        if line[0] != "chr":
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            gene_name = line[6].upper()
            tranId = line[7]
            gene_strand = line[5]
            gene_expr = line[12]
            gene_exprq = line[13]
            if chrom not in genes_dict:
                genes_dict[chrom] = []
            genes_dict[chrom].append([start, end, gene_name, gene_strand, tranId ,gene_expr, gene_exprq])


    df_MRCA = pd.read_csv("../../MRCA/results/Summary_divisions_with_symmetrical_no_tetra.txt", sep=",")
    print(df_MRCA.head())
    samples = df_MRCA["sample"]
    for sample in samples:
        mut_file = open("../../data/mutations/" + sample + ".nodMat")
        print(sample)
        outname = "../../data/mutations_vs_genes/" + sample + ".with_gene_annot.nodMat"
        outfile = open(outname, "a")
        for line in mut_file:
            line=line.strip("\r\n")
            data = line.split(",")
            if data[0] == 'chr':
                outfile.write(line + "," + "geneName,geneStrand,tranId,expression,expression_q,inCensus,inMAPK\n")
            else:
                mut_chr = data[0]
                mut_coord = int(data[1])
                gene_name=""
                gene_strand=""
                transcript_id=""
                gene_expression=""
                gene_expression_q=""
                gene_in_census=0
                gene_in_MAPK=0
                for gene in genes_dict[mut_chr]:
                    if (mut_coord >= gene[0] and mut_coord <= gene[1]):
                        gene_name=gene[2]
                        gene_strand=gene[3]
                        transcript_id=gene[4]
                        gene_expression=gene[5]
                        gene_expression_q=gene[6]
                        if gene[2] in genes_census:
                            gene_in_census = 1
                        if gene[2] in genes_MAPK:
                            gene_in_MAPK = 1

                outfile.write(line + "," + gene_name + "," + gene_strand + "," + transcript_id  + "," + gene_expression  + "," + gene_expression_q + "," + str(gene_in_census) + "," + str(gene_in_MAPK) + "\n")    

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)


			
