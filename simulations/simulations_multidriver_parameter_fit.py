import numpy as np
import pandas as pd
import os
from simulation_plots import get_mrca_perc, get_mrca_and_all_4_list
from analytics import mrca_mse

conf_interval_df = pd.read_csv('C3H_confidence_intervals_poisson.csv')
conf_interval_groups = conf_interval_df.groupby('Gene_name')
for gene, grp in conf_interval_groups:
    srted_grp = grp.sort_values(by='division')
    grp_total = grp['N_samples'].astype(float).sum()
    CF = 95
    for gene, grp in conf_interval_groups:
        srted_grp = grp.sort_values(by='division')
        grp_total = grp['N_samples'].astype(float).sum()
        if gene == 'Egfr':
            egfr_low = srted_grp[f'Lower_CF_{CF}'].astype(float).values / grp_total
            egfr_high = srted_grp[f'Upper_CF_{CF}'].astype(float).values / grp_total
        elif gene == 'Braf':
            braf_low = srted_grp[f'Lower_CF_{CF}'].astype(float).values / grp_total
            braf_high = srted_grp[f'Upper_CF_{CF}'].astype(float).values / grp_total
        elif gene == 'Hras1':
            hras1_low = srted_grp[f'Lower_CF_{CF}'].astype(float).values / grp_total
            hras1_high = srted_grp[f'Upper_CF_{CF}'].astype(float).values / grp_total
        elif gene == 'Hras0':
            hras0_low = srted_grp[f'Lower_CF_{CF}'].astype(float).values / grp_total
            hras0_high = srted_grp[f'Upper_CF_{CF}'].astype(float).values / grp_total

egfr_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(egfr_low, egfr_high)]
braf_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(braf_low, braf_high)]
hras1_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(hras1_low, hras1_high)]
hras0_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(hras0_low, hras0_high)]



outcols_full = ['# drivers', 's', 'efr', 'repair', 'duplex rate',
           'MRCA 0', 'MRCA 1', 'MRCA 2', 'MRCA 3', 'MRCA 4', 'MRCA 5-300',
           'egfr_c3h_sd_mse', 'braf_c3h_sd_mse', 'hras1_c3h_sd_mse', 'hras0_c3h_sd_mse']

outcols_params = ['# drivers', 's', 'e', 'r', 'u']


c = 0

params_fit_hras0 = [] # fill with (r, u, e, min s, max s)
params_fit_hras1 = []
params_fit_braf = []
params_fit_egfr = []
params_fit_egfr_shape = []

params_not_run = []
params_failed = []

rows_match_hras0_sd = []
rows_match_hras1_sd = []
rows_match_egfr_sd = []
rows_match_egfr_shape_sd = []
rows_match_braf_sd = []

outdir = 'simulation_grid_multidriver/'
if not os.path.isdir(outdir):
    os.mkdir(outdir)
cell_path = 'cells/grcm38.p6_rounded_driver_loci_freq_5/2025-10-10_11-22/ '

driver_count_range = np.arange(2, 6, 1)
s_range = np.arange(0.4, 1, 0.1)
efr_range = np.arange(0.05, 0.55, 0.1)
r_range = np.arange(0.2, 0.9, 0.1)
u_range = np.arange(0.05, 0.8, 0.1)


for nd in driver_count_range:
    for s in s_range:
        other_s = 0.2
        if nd == 2:
            if round(s, 1) in [0.4, 0.5, 0.6]:
                other_s = 0.25
            elif round(s, 1) in [0.7, 0.8]:
                other_s = 0.3
            elif round(s, 1) in [0.9, 1]:
                other_s = 0.5

        lesions_path = f'{cell_path}/lesions_{nd}_drivers_0_other_smax_{round(s, 1)}_s_{other_s}/simulation_results/'

        all_rows = []
        print(lesions_path)
        for e in efr_range:
            for r in r_range:
                for u in u_range:
                    params_row = [nd, round(s, 1), round(e, 2), round(r, 1), round(u, 2)]
                    run = True
                    simulation_path = f'{lesions_path}/epistasis_syn_driver_on_nontx_strands_efr_{round(e, 2)}_eft_0.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
                    if not os.path.isdir(simulation_path):
                        simulation_path = f'{lesions_path}/epistasis_syn_driver_on_tx_strands_efr_{round(e, 2)}_eft_0.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
                        if not os.path.isdir(simulation_path):
                            run = False
                            params_not_run.append(params_row)
                    if run:
                        subdirs = sorted(os.listdir(f'{simulation_path}'), reverse=True)
                        completed_found = False
                        for d in subdirs:
                            subdir = f'{simulation_path}/{d}/'
                            if os.path.isdir(subdir) and 'combined_clone_info_no_filter.csv' in os.listdir(subdir):
                                print(subdir)
                                completed_found = True
                                row = []

                                row.extend(params_row)


                                clones_df = pd.read_csv(f'{subdir}/combined_clone_info_no_filter.csv')
                                mrca_list, all_4_list = get_mrca_and_all_4_list(clones_df)
                                mrca_row = get_mrca_perc(mrca_list)
                                # all_4_perc_row = get_all_4_perc(mrca_list, all_4_list)
                                row.extend(mrca_row)

    
                                egfr_mse_sd = mrca_mse(mrca_row, egfr_range_sd)
                                braf_mse_sd = mrca_mse(mrca_row, braf_range_sd)
                                hras_1_mse_sd = mrca_mse(mrca_row, hras1_range_sd)
                                hras_0_mse_sd = mrca_mse(mrca_row, hras0_range_sd)
                                row.extend([egfr_mse_sd, braf_mse_sd, hras_1_mse_sd, hras_0_mse_sd])

                                all_rows.append(row)

                                if mrca_row[0] >= egfr_range_sd[0][0] and mrca_row[0] <= egfr_range_sd[0][1]:
                                    if mrca_row[1] >= egfr_range_sd[1][0] and mrca_row[1] <= egfr_range_sd[1][1]:
                                        if mrca_row[2] >= egfr_range_sd[2][0] and mrca_row[2] <= egfr_range_sd[2][1]:
                                            if mrca_row[3] >= egfr_range_sd[3][0] and mrca_row[3] <= egfr_range_sd[3][1]:
                                                if mrca_row[4] >= egfr_range_sd[4][0] and mrca_row[4] <= egfr_range_sd[4][
                                                    1]:
                                                    if mrca_row[-1] >= egfr_range_sd[5][0] and mrca_row[-1] <= \
                                                            egfr_range_sd[5][1]:
                                                        rows_match_egfr_sd.append(row)
                                                        params_fit_egfr.append(params_row)
                                                        if mrca_row[2] > mrca_row[1]:
                                                            rows_match_egfr_shape_sd.append(row)
                                                            params_fit_egfr_shape.append(params_row)

                                if mrca_row[0] >= braf_range_sd[0][0] and mrca_row[0] <= braf_range_sd[0][1]:
                                    if mrca_row[1] >= braf_range_sd[1][0] and mrca_row[1] <= braf_range_sd[1][1]:
                                        if mrca_row[2] >= braf_range_sd[2][0] and mrca_row[2] <= braf_range_sd[2][1]:
                                            if mrca_row[3] >= braf_range_sd[3][0] and mrca_row[3] <= braf_range_sd[3][1]:
                                                if mrca_row[4] >= braf_range_sd[4][0] and mrca_row[4] <= braf_range_sd[4][
                                                    1]:
                                                    if mrca_row[-1] >= braf_range_sd[5][0] and mrca_row[-1] <= \
                                                            braf_range_sd[5][1]:
                                                        rows_match_braf_sd.append(row)
                                                        params_fit_braf.append(params_row)


                                if mrca_row[0] >= hras1_range_sd[0][0] and mrca_row[0] <= hras1_range_sd[0][1]:
                                    if mrca_row[1] >= hras1_range_sd[1][0] and mrca_row[1] <= hras1_range_sd[1][1]:
                                        if mrca_row[2] >= hras1_range_sd[2][0] and mrca_row[2] <= hras1_range_sd[2][1]:
                                            if mrca_row[3] >= hras1_range_sd[3][0] and mrca_row[3] <= hras1_range_sd[3][1]:
                                                if mrca_row[4] >= hras1_range_sd[4][0] and mrca_row[4] <= hras1_range_sd[4][
                                                    1]:
                                                    if mrca_row[-1] >= hras1_range_sd[5][0] and mrca_row[-1] <= \
                                                            hras1_range_sd[5][1]:
                                                        rows_match_hras1_sd.append(row)
                                                        params_fit_hras1.append(params_row)

                                if mrca_row[0] >= hras0_range_sd[0][0] and mrca_row[0] <= hras0_range_sd[0][1]:
                                    if mrca_row[1] >= hras0_range_sd[1][0] and mrca_row[1] <= hras0_range_sd[1][1]:
                                        if mrca_row[2] >= hras0_range_sd[2][0] and mrca_row[2] <= hras0_range_sd[2][1]:
                                            if mrca_row[3] >= hras0_range_sd[3][0] and mrca_row[3] <= hras0_range_sd[3][1]:
                                                if mrca_row[4] >= hras0_range_sd[4][0] and mrca_row[4] <= hras0_range_sd[4][
                                                    1]:
                                                    if mrca_row[-1] >= hras0_range_sd[-1][0] and mrca_row[-1] <= \
                                                            hras0_range_sd[-1][1]:
                                                        rows_match_hras0_sd.append(row)
                                                        params_fit_hras0.append(params_row)

                                break
                        if not completed_found:
                            params_failed.append(params_row)

        df = pd.DataFrame(all_rows, columns=outcols_full)
        df.to_csv(f'{outdir}/s_{round(s, 2)}_{nd}_drivers.csv', index=False)



df = pd.DataFrame(rows_match_braf_sd, columns=outcols_full)
df.to_csv(f'{outdir}/matching_braf_sd_multidriver_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_braf, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_braf_multidriver_cf_{CF}.csv', index=False)


df = pd.DataFrame(rows_match_hras1_sd, columns=outcols_full)
df.to_csv(f'{outdir}/matching_hras1_multidriver_sd_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_hras1, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_hras1_multidriver_cf_{CF}.csv', index=False)


df = pd.DataFrame(rows_match_hras0_sd, columns=outcols_full)
df.to_csv(f'{outdir}/matching_hras0_multidriver_sd_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_hras0, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_hras0_multidriver_cf_{CF}.csv', index=False)


df = pd.DataFrame(rows_match_egfr_sd, columns=outcols_full)
df.to_csv(f'{outdir}/matching_egfr_multidriver_sd_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_egfr, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_egfr_multidriver_cf_{CF}.csv', index=False)


df = pd.DataFrame(rows_match_egfr_shape_sd, columns=outcols_full)
df.to_csv(f'{outdir}/matching_egfr_shape_multidriver_sd_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_egfr_shape, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_egfr_shape_multidriver_cf_{CF}.csv', index=False)

df = pd.DataFrame(params_failed, columns=outcols_params)
df.to_csv(f'{outdir}/params_failed_runs_multidriver.csv', index=False)
df = pd.DataFrame(params_not_run, columns=outcols_params)
df.to_csv(f'{outdir}/params_missing_multidriver.csv', index=False)















