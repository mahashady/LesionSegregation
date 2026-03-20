#$ -S /usr/bin/python

import glob
import pandas as pd
from pathlib import Path

def main():
    in_files = glob.glob("../../data/raw_data_all/*.nodMat")
    out_dir = Path("../../data/encoded_data_all")
    out_dir.mkdir(parents=True, exist_ok=True)

    for file in in_files:
        sample = Path(file).stem
        print(sample)
        out_file = out_dir / f"{sample}.hmm"

        df = pd.read_csv(file, sep=",")
        df = df[df['isIndel'] == 0].copy()  #select only non-Indels


        ##############################
        ## 1. Multiallelic coding
        ##############################

        # create Boolean DataFrame showing if each of nucleotides has count >= 2
        ge2 = df[["au","cu","gu","tu"]] >= 2
        # Count how many nucleotides per row have values >=2
        count_ge2 = ge2.sum(axis=1)
        # Condition: at least 3 of 4 nucleotides have more than 2 alt reads - multiallelic sites
        cond_multi = count_ge2 >= 3

        # check sum of counts != ref+alt
        sum_check = (df["au"] + df["cu"] + df["gu"] + df["tu"]) != (df["refCount"] + df["altCount"])

        # Final condition for code_multi = "M"
        df["code_multi"] = "B"  # default biallelic
        df.loc[sum_check & cond_multi, "code_multi"] = "M" #multiallelic

        ##############################
        ## 2. Asymmetry coding
        ##############################
        df["code_as"] = "other"
        df.loc[df["refN"] == "T", "code_as"] = "T>N"
        df.loc[df["refN"] == "A", "code_as"] = "A>N"
                                
        df.to_csv(out_file, index=False)


if __name__=="__main__":
    main()