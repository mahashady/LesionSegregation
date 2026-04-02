from pathlib import Path
import pandas as pd
from utils.LAD_simul_utils import build_likelihood_df

def main():
    out_file = Path("/home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/LAD_inference/likelihood_matrix_for_LAD_inference.csv")

    likelihood_df = build_likelihood_df(
        M_values=range(2, 16),   
        D_values=[1, 2, 3, 4, 5],
        n_iter=20000,
        seed=42,
        pseudocount=0.0,
        complete_grid=True,
    )

    likelihood_df.to_csv(out_file, sep="\t", index=False)
    print(f"Saved {len(likelihood_df)} rows to {out_file}")

if __name__ == "__main__":
    main()