import pysam
import pandas as pd
import glob 
import os 
import gzip
from bgreference import hg19

NUCLEOTIDES = {"A", "T", "G", "C"}
CHROMOSOMES = {str(i) for i in range(1, 23)}

def make_empty_dict():
    return {
        ref: {chrom: {} for chrom in CHROMOSOMES}
        for ref in NUCLEOTIDES
    }


def make_empty_chrom_dict(chromosomes=CHROMOSOMES):
    return {chrom: {} for chrom in chromosomes}


def read_vcf_to_mut_dict(vcf_file, sample_name=None):
    """
    Reads a VCF using pysam and returns mut_dict:
    dict[ref][chrom][pos] = [mutations]
    """

    mut_dict = make_empty_dict()

    vcf = pysam.VariantFile(vcf_file)

    # auto-detect sample if not provided
    if sample_name is None:
        samples = list(vcf.header.samples)
        if len(samples) == 0:
            raise ValueError("No samples found in VCF")
        sample_name = samples[0]

    for rec in vcf:
        chrom = rec.contig
        pos = rec.pos
        ref = rec.ref

        # skip multiallelic rows if needed
        if len(rec.alts) != 1:
            continue

        alt = rec.alts[0]

        if chrom not in CHROMOSOMES:
            continue
        if ref not in NUCLEOTIDES or alt not in NUCLEOTIDES:
            continue

        sample_data = rec.samples[sample_name]

        tumour_ad = sample_data.get("AD", "")
        tumour_af = sample_data.get("AF", "")

        mutation = f"{alt}:{tumour_ad}:{tumour_af}"

        mut_dict[ref][chrom].setdefault(pos, []).append(mutation)

    return mut_dict

def load_germline(file):
    germline_dict = make_empty_chrom_dict()
    if os.path.isfile(file):
        f=gzip.open(file,'rb')
        file_content=f.readlines()
        for line in file_content:
            line = line.decode("utf-8")
            if line[0] != "#":
                data = line.strip("\r\n").split("\t")
                chrom, pos, ref, alt, gt_n = data[0], int(data[1]), data[3], data[4], data[9]                
                if chrom in CHROMOSOMES and gt_n.split(":")[0] == "0/1":
                    germline_dict[chrom][pos] = ref + ":" + alt + ":" + gt_n
    else:
        print("Germline file is absent")
    return(germline_dict)    

def get_multi_sample_list(samples_table, subset):
    df = pd.read_csv(samples_table, sep="\t")
    df["n_multi"] = pd.to_numeric(df["n_multi"], errors="coerce")

    sample_ids = (
        df.loc[df["n_multi"] > 0, "patientIdentifier"]
        .dropna()
        .drop_duplicates()
        .tolist()
    )

    print(f"Number of samples with multiallelic sites in {subset} subset: {len(sample_ids)}")
    return sample_ids 

def nearby_germline_positions(germline_dict_chr, pos, window=150):
    return [
        gpos for gpos in germline_dict_chr
        if pos - window <= gpos <= pos + window
    ]


def write_multiallelic_sites(
    sample,
    mut_dict,
    germline_dict,
    outfile_handle,
    summary_handle,
    summary_window=200,
    germline_window=150,
):
    n_multi = 0
    n_multi_with_germline = 0

    for chromosome in CHROMOSOMES:
        for nucleotide in NUCLEOTIDES:
            for pos, mutations in sorted(mut_dict[nucleotide][chromosome].items()):
                if len(mutations) <= 1:
                    continue

                n_multi += 1
                context = hg19(chromosome, pos - 1, 3)

                outfile_handle.write(
                    f"{sample}\t{chromosome}\t{nucleotide}\t{pos}\t{context}"
                )

                for mut in mutations:
                    outfile_handle.write(f"\t{mut}")

                germline_list = nearby_germline_positions(
                    germline_dict[chromosome],
                    pos=pos,
                    window=germline_window,
                )

                if len(germline_list) > 0:
                    n_multi_with_germline += 1
                    summary_handle.write(
                        f"{sample}\t{chromosome}:{pos-summary_window}-{pos+summary_window}\n"
                    )

                    germline_strings = [
                        f"{gpos}:{germline_dict[chromosome][gpos]}"
                        for gpos in sorted(germline_list)
                    ]
                    outfile_handle.write("\t" + ";".join(germline_strings))
                else:
                    outfile_handle.write("\tNA")

                outfile_handle.write("\n")

    return n_multi, n_multi_with_germline

def find_single_file(pattern, label):
    matches = glob.glob(pattern)
    if len(matches) == 0:
        return None
    if len(matches) > 1:
        print(f"Warning: multiple {label} files found for pattern {pattern}. Using first one.")
    return matches[0]    