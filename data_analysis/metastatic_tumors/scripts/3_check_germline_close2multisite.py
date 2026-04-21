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
from utils.utils import *
from pathlib import Path


def main(args):
    subset = "enriched" #NONenriched or enriched

    #define path to HARTWIG data
    HARTWIG_PATH = "/data/bbg/datasets/hartwig/"

    # define results root
    SCRIPT_DIR = Path(__file__).resolve().parent
    PROJECT_ROOT = SCRIPT_DIR.parent
    RESULTS_DIR = PROJECT_ROOT / "results" / "multi-allelic_sites"
    
    samples_table = PROJECT_ROOT / "results" / "enrichment" / f"all_{subset}.chemo_alkyl_immuno.txt"
    output_dir = RESULTS_DIR / f"multi_sites_by_sample_{subset}.chemo_alkyl_immuno"
    summary_file = RESULTS_DIR / f"summary_multi_sites_with_germline_{subset}_chemo.alkyl.immuno.txt"

    os.makedirs(output_dir, exist_ok=True)

    multi_samples_list = get_multi_sample_list(samples_table, subset=subset)

    total_multi = 0
    total_multi_with_germline = 0

    with open(summary_file, "w") as summary_outfile:
        for sample in multi_samples_list:
            print(sample)

            somatic_pattern = (
                f"{HARTWIG_PATH}20230914/somatic/{sample}T/purple/"
                f"{sample}*.purple.somatic.vcf.gz"
            )
            somatic_vcf = find_single_file(somatic_pattern, label="somatic VCF")
            if somatic_vcf is None:
                print(f"Somatic VCF is absent for {sample}")
                continue

            outfile_path = os.path.join(output_dir, f"{sample}_multi_sites.txt")
            if os.path.isfile(outfile_path):
                print(f"Output already exists, skipping: {outfile_path}")
                continue

            print("Creating file")
            mut_dict = read_vcf_to_mut_dict(somatic_vcf, sample_name = sample+"T")
            print("Mutation dictionary created")

            if subset == "enriched":
                germline_pattern = f"{HARTWIG_PATH}20230914/germline/*{sample}.annotated.vcf.gz"
            else:
                germline_pattern = f"{HARTWIG_PATH}20230914/germline/2nd_batch/*{sample}.annotated.vcf.gz"

            germline_vcf = find_single_file(germline_pattern, label="germline VCF")
            if germline_vcf is None:
                print(f"Germline file is absent for {sample}")
                germline_dict = make_empty_chrom_dict()
            else:
                germline_dict = load_germline(germline_vcf)

            with open(outfile_path, "w") as outfile:
                n_multi, n_multi_with_germline = write_multiallelic_sites(
                    sample=sample,
                    mut_dict=mut_dict,
                    germline_dict=germline_dict,
                    outfile_handle=outfile,
                    summary_handle=summary_outfile,
                )

            total_multi += n_multi
            total_multi_with_germline += n_multi_with_germline

    print(f"All multi sites = {total_multi}")
    print(f"Multi sites with germline = {total_multi_with_germline}")

        

  

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)
