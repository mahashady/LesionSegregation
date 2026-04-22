#!/usr/bin/env python3

import argparse
import glob
import os
import pysam
 from utils.phasing_utils import *

# =========================
# User-editable paths
# =========================
MULTI_SITES_DIRS = {
    "enriched": "../results/multi-allelic_sites/sites/multi_sites_by_sample_enriched_chemo.alkyl.immuno",
    "NONenriched": "../results/multi-allelic_sites/sites/multi_sites_by_sample_NONenriched_chemo.alkyl.immuno",
}

BAM_BASE_DIRS = {
    "enriched": "/workspace/datasets/hartwig/20230914/scripts/minibam/results",
    "NONenriched": "/workspace/datasets/hartwig/20230914/scripts/minibam/results/batch2",
}

PHASING_OUT_DIR = "../results/phasing"




def main(args):
    subset = args.subset
    phased = 0
    anti = 0
    no_phasing = 0
    all_multi = 0
    no_germline = 0

    input_pattern = os.path.join(
        MULTI_SITES_DIRS[subset],
        "*_multi_sites.txt",
    )
    files = glob.glob(input_pattern)

    os.makedirs(PHASING_OUT_DIR, exist_ok=True)

    out_path = os.path.join(
        PHASING_OUT_DIR,
        f"phased_multi_{subset}_chemo.alkyl.immuno.txt",
    )
    out_summary_path = os.path.join(
        PHASING_OUT_DIR,
        f"phased_summary_multi_{subset}_chemo.alkyl.immuno.txt",
    )

    with open(out_path, "w") as outfile, open(out_summary_path, "w") as outfile_summary:
        for file_name in files:
            sample = os.path.basename(file_name).split("_multi")[0]
            print(sample)

            with open(file_name) as infile:
                for line in infile:
                    data = line.strip("\r\n").split("\t")
                    all_multi += 1

                    outfile.write(sample)
                    outfile.write("\t")
                    outfile.write(data[1] + "_" + data[3])

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

                        germline_snvs = []
                        print("Multi pos=", str(pos_multi))

                        for germline_idx in range(n_germlines):
                            pos_germline = int(germlines[germline_idx].split(":")[0])
                            print("Germline pos=", str(pos_germline))

                            ref_germline = germlines[germline_idx].split(":")[1]
                            alt_germline = germlines[germline_idx].split(":")[2]

                            if len(ref_germline) == 1 and len(alt_germline) == 1:
                                try:
                                    result = phase(
                                        sample,
                                        chrom,
                                        pos_multi,
                                        ref_multi,
                                        alt1_multi,
                                        alt2_multi,
                                        pos_germline,
                                        ref_germline,
                                        alt_germline,
                                        subset,
                                    )
                                    germline_snvs.append(result)
                                    print(result)
                                except Exception as e:
                                    print(
                                        f"Problem with BAM files for sample={sample}, "
                                        f"chrom={chrom}, pos_multi={pos_multi}, "
                                        f"germline={germlines[germline_idx]}: {e}"
                                    )

                        if germline_snvs == []:
                            outfile.write("\tonly_germline_indels")
                            no_phasing += 1
                            outfile.write("\tNA")
                        else:
                            outfile.write("\t" + str(germline_snvs))
                            germline1 = germline_snvs[0]
                            status = classify_phasing(germline1)

                            if status == "Phased":
                                phased += 1
                            elif status == "Anti":
                                anti += 1
                            else:
                                no_phasing += 1

                            outfile.write("\t" + status)

                    outfile.write("\n")

        print(all_multi)
        print(no_germline)
        outfile_summary.write(
            "Phased\t" + str(phased) + "\n"
            + "Anti\t" + str(anti) + "\n"
            + "No phasing\t" + str(no_phasing)
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--subset",
        required=True,
        choices=["enriched", "NONenriched"],
        help="Which subset to analyze",
    )
    args = parser.parse_args()
    main(args)