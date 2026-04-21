#$ -S /usr/bin/python
import argparse
import glob
import gzip

from collections import Counter
from pathlib import Path


""" 
"""
nucleotides = ['A','T','G','C']
chromosomes = [str(i) for i in range(1,23)]

root_folder = "/home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors"

def make_empty_dict():
    counts_dict = {}
    for nucleotide in nucleotides:
        counts_dict[nucleotide] = {}
        for chromosome in chromosomes:
            counts_dict[nucleotide][chromosome] = []
    return counts_dict


def main():
    #PATH CONTAINING SOMATIC VCFs 
    files = glob.glob("/data/bbg/datasets/hartwig/20230914/somatic/*/purple/*.purple.somatic.vcf.gz")
    outfile_full = Path(root_folder) / "results/ALL_Hartwig_bi_multi.by_chrom.txt"
    outfile_short = Path(root_folder) / "results/ALL_Hartwig_bi_multi.txt"

    outfile_full.parent.mkdir(parents=True, exist_ok=True)
    outfile_short.parent.mkdir(parents=True, exist_ok=True)

    with open(outfile_full, "w") as full_f, open(outfile_short, "w") as short_f:
        full_f.write(
            "sample\tchr\tA>N\tT>N\tG>N\tC>N\t"
            "A>N_triallelic\tT>N_triallelic\tG>N_triallelic\tC>N_triallelic\t"
            "A>N_cuatro\tT>N_cuatro\tG>N_cuatro\tC>N_cuatro\n"
        )  
        short_f.write("sample\tn_biallelic\tn_multi\n")  

        for i, file in enumerate(files, start=1):
            sample = Path(file).parts[-3]
            print(i)
            print(sample)
            counts_dict = make_empty_dict()
            with gzip.open(file, "rt") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    data = line.strip("\r\n").split("\t")
                    chrom, pos, ref, alt = data[0], data[1], data[3], data[4]
                    if chrom in chromosomes and ref in nucleotides and alt in nucleotides:
                        counts_dict[ref][chrom].append(pos)
        
            full_result = [0,0]
            for chromosome in chromosomes:
                result = {'A':[],'T':[],'G':[],'C':[]}
                for nucleotide in nucleotides:
                    n_bi_multi = Counter(Counter(counts_dict[nucleotide][chromosome]).values())
                    n_biallelic = n_bi_multi[1]
                    n_triallelic = n_bi_multi[2] 
                    n_cuatro = n_bi_multi[3]
                    result[nucleotide] = [n_biallelic, n_triallelic, n_cuatro]
                full_f.write(
                    sample + "\t" + chromosome + "\t"
                    + str(result["A"][0]) + "\t" + str(result["T"][0]) + "\t"
                    + str(result["G"][0]) + "\t" + str(result["C"][0]) + "\t"
                    + str(result["A"][1]) + "\t" + str(result["T"][1]) + "\t"
                    + str(result["G"][1]) + "\t" + str(result["C"][1]) + "\t"
                    + str(result["A"][2]) + "\t" + str(result["T"][2]) + "\t"
                    + str(result["G"][2]) + "\t" + str(result["C"][2]) + "\n"
                )
                all_bi = result['A'][0] + result['T'][0] + result['G'][0] + result['C'][0]
                all_multi = result['A'][1] + result['T'][1] + result['G'][1] + result['C'][1] + result['A'][2] + result['T'][2] + result['G'][2] + result['C'][2]
                full_result[0] += all_bi
                full_result[1] += all_multi
            short_f.write(sample + "\t" + str(full_result[0]) + "\t" + str(full_result[1]) + "\n")
        


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main()


			
















