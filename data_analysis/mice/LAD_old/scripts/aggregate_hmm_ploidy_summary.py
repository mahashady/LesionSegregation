import os
import pandas as pd

# Input list of samples
samples_file = "../../data/list_of_samples.txt"
output_file = "../results/HMM/HMM_ploidy/Ploidy_summary.tsv"

with open(samples_file) as f:
    samples = [line.strip() for line in f]

results = []

for sample in samples:
    vaf_file = f"../results/HMM/HMM_ploidy/{sample}.hmm.Clonesize"
    bw_file = f"../results/HMM/HMM_ploidy/{sample}.hmm.BW"
    as_file = f"../results/HMM/HMM_ploidy/{sample}.hmm.CloneAS"

    # If files missing, no mixture fitting was possible -> 1 clone
    if not all(os.path.exists(f) for f in [vaf_file, bw_file, as_file]):
        results.append([sample, '1_clone', '1_clone', '1_clone', '1_clone'])
        continue

    # Read VAF properties
    vaf_df = pd.read_csv(vaf_file, sep=' ')
    vaf_counts = vaf_df['vafMed_vafboundaries_Distr_size']  
    print(vaf_counts)

    # Determine number of clones
    Ncl = 2
    if vaf_counts[4] < 800 or vaf_counts[5] < 800:
        Ncl = 1
    elif Ncl == 2 and vaf_counts[3] > vaf_counts[4]:
        Ncl = 'ER'

    # --- Read BW HMM ---
    bw_df = pd.read_csv(bw_file, sep=' ', header=None)
    bwp = 'Good'
    S1 = {}
    S2 = {}

    # Check emissions from first row 
    row1 = bw_df.iloc[0]
    print(row1)
    if (row1[1] < 0.67 or row1[4] < 0.67 or row1[3] > 0.33 or row1[6] > 0.33
        or row1[2] > 0.6 or row1[5] > 0.6 or row1[2] < 0.4 or row1[5] < 0.4):
        bwp = 'Bad'

    # Compute states proportions from second row
    row2 = bw_df.iloc[1]
    print(row2)

    all1 = sum(row2[1:4])
    all2 = sum(row2[4:7])

    if all1 < 800 and all2 < 800:
        Ncl = 1
    for i in range(3):
        S1[i+1] = row2[i+1] / all1
        S2[i+1] = row2[i+4] / all2
        if row2[i+1] < 100 or row2[i+4] < 100:
            bwp = 'Bad'
    print(S1)
    print(S2)
    # --- Read CloneAS ---
    as_df = pd.read_csv(as_file, sep=' ', header=None)
    combined_labels = as_df.iloc[0].apply(lambda x: x[1]+x[4])  # first row, combine characters
    print(as_df.iloc[0])
    counts_as = as_df.iloc[1].astype(int)
    total_as = counts_as.sum()
    S12 = {combined_labels[i]: counts_as[i]/total_as for i in range(len(counts_as))}
    print("combined")
    print(S12)

    # --- Compute observed/expected ratios ---
    o_e = {}
    exp = {}
    for f in range(1,4):
        for s in range(1,4):
            fs = f*10 + s  # create label like '20'
            print(fs)
            fs_str = f"{f}{s}"
            exp[fs_str] = S1[f]*S2[s]
            o_e[fs_str] = S12.get(fs_str, 0)/exp[fs_str] if exp[fs_str] > 0 else 0

    print(exp)
    # --- Assign asymmetry class ---
    ASS = 'tetra'
    if o_e['11']>1.5 and o_e['33']>1.5 and o_e['31']<0.6 and o_e['13']<0.6:
        ASS = 'Dip*'
    if o_e['11']>1.8 and o_e['33']>1.8 and o_e['31']<0.6 and o_e['13']<0.6:
        ASS = 'Dip'
    if o_e['31']>2 and o_e['13']>2:
        ASS = 'Symmetric'

    results.append([sample, bwp, o_e['11'],o_e['33'],o_e['31'], o_e['13'], ASS])

# Save summary
pd.DataFrame(results, columns=['Sample', 'BW_status', 'o_e_11', 'o_e_33' , 'o_e_31', 'o_e_13', 'ASS']) \
    .to_csv(output_file, sep='\t', index=False)