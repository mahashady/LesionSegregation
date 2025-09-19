import numpy as np
import pandas as pd
import os
from simulation_plots import get_mrca_perc, get_mrca_list, get_mrca_and_all_4_list, get_all_4_perc
from analytics import mrca_mse



conf_interval_df = pd.read_csv('/n/data2/hms/dbmi/sunyaev/lab/maha/drivers/cancer_simulations/C3H_confidence_intervals_poisson.csv')
conf_interval_groups = conf_interval_df.groupby('Gene_name')
for gene, grp in conf_interval_groups:
    srted_grp = grp.sort_values(by='division')
    grp_total = grp['N_samples'].astype(float).sum()
    CF = 95

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

conf_interval_all_4_df = pd.read_csv('/n/data2/hms/dbmi/sunyaev/lab/maha/drivers/cancer_simulations/all_4_survived_mrca_1_cf.csv')
conf_interval_all_4_groups = conf_interval_all_4_df.groupby('drivers')
for gene, grp in conf_interval_all_4_groups:
    grp_total = grp['n'].astype(float).sum()
    subgrp = grp.loc[grp['type'] == 'all_cells_survived']
    low = subgrp['lower_cf'].values[0]/grp_total
    hi = min(1, subgrp['upper_cf'].values[0]/grp_total)
    if gene == 'Egfr':
        egfr_all_4_cf = (low, hi)
    elif gene == 'Braf':
        braf_all_4_cf = (low, hi)
    elif gene == 'Hras1':
        hras1_all_4_cf = (low, hi)
    elif gene == 'Hras0':
        hras0_all_4_cf = (low, hi)


outcols_full = ['s', 'efr', 'repair', 'duplex rate', 'r-u',
           'MRCA 0', 'MRCA 1', 'MRCA 2', 'MRCA 3', 'MRCA 4', 'MRCA 5-300',
           'egfr_c3h_sd_mse', 'braf_c3h_sd_mse', 'hras1_c3h_sd_mse', 'hras0_c3h_sd_mse',
                'MRCA 0 all 4', 'MRCA 1 all 4', 'MRCA 2 all 4', 'MRCA 3 all 4']
outcols_params = ['s', 'e', 'r', 'u', 'r-u']



params_fit_hras0 = [] # fill with (r, u, e, min s, max s)
params_fit_hras1 = []
params_fit_braf = []
params_fit_egfr = []
params_fit_egfr_shape = []

params_fit_hras0_all_4 = []
params_fit_hras1_all_4 = []
params_fit_braf_all_4 = []
params_fit_egfr_all_4 = []

params_fit_hras0_all_4_only = []
params_fit_hras1_all_4_only = []
params_fit_braf_all_4_only = []
params_fit_egfr_all_4_only = []

params_not_run = []
params_failed = []

rows_match_hras0_sd = []
rows_match_hras1_sd = []
rows_match_egfr_sd = []
rows_match_egfr_shape_sd = []
rows_match_braf_sd = []


outdir = '/n/data2/hms/dbmi/sunyaev/lab/maha/drivers/simulation_grid_2driver/'
if not os.path.isdir(outdir):
    os.mkdir(outdir)
cell_path = '/n/data2/hms/dbmi/sunyaev/lab/maha/drivers/cells/grcm38.p6_rounded_autosomes_driver_loci_freq_5/2024-06-02_13-23/'

s_list = [(0.2, 0.1) ,(0.3, 0.1), (0.4, 0.25), (0.5, 0.25), (0.6, 0.25), (0.7, 0.3), (0.8, 0.3)]
efr_range = np.arange(0.05, 0.55, 0.1)
r_range = np.arange(0.2, 0.9, 0.1)
u_range = np.arange(0.05, 0.8, 0.1)


for s, s1 in s_list:
    lesions_path = f'{cell_path}/lesions_2_drivers_0_other_smax_{round(s, 1)}_s_{round(s1, 2)}/simulation_results/'
    all_rows = []
    for e in efr_range:
        for r in r_range:
            for u in u_range:
                run = True

                simulation_path = f'{lesions_path}/epistasis_syn_driver_on_nontx_strands_efr_{round(e, 2)}_eft_0.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
                params_row = [round(s, 1), round(e, 2), round(r, 1), round(u, 2), round(r-u, 2)]

                if not os.path.isdir(simulation_path):
                    simulation_path = f'{lesions_path}/epistasis_syn_driver_on_tx_strands_efr_{round(e, 2)}_eft_0.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
                    if not os.path.isdir(simulation_path):
                        print(f'not run: {simulation_path}')
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
                            all_4_perc_row = get_all_4_perc(mrca_list, all_4_list)
                            row.extend(mrca_row)

                            egfr_mse_sd = mrca_mse(mrca_row, egfr_range_sd)
                            braf_mse_sd = mrca_mse(mrca_row, braf_range_sd)
                            hras_1_mse_sd = mrca_mse(mrca_row, hras1_range_sd)
                            hras_0_mse_sd = mrca_mse(mrca_row, hras0_range_sd)
                            row.extend([egfr_mse_sd, braf_mse_sd, hras_1_mse_sd, hras_0_mse_sd])

                            row.extend(all_4_perc_row)

                            all_rows.append(row)

                            if all_4_perc_row[1] >= egfr_all_4_cf[0] and all_4_perc_row[1] <= egfr_all_4_cf[1]:
                                params_fit_egfr_all_4_only.append(params_row)
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
                                                    if all_4_perc_row[1] >= egfr_all_4_cf[0] and all_4_perc_row[1] <= egfr_all_4_cf[1]:
                                                        params_fit_egfr_all_4.append(params_row)

                            if all_4_perc_row[1] >= braf_all_4_cf[0] and all_4_perc_row[1] <= braf_all_4_cf[1]:
                                params_fit_braf_all_4_only.append(params_row)
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
                                                    if all_4_perc_row[1] >= braf_all_4_cf[0] and all_4_perc_row[1] <= braf_all_4_cf[1]:
                                                        params_fit_braf_all_4.append(params_row)

                            if all_4_perc_row[1] >= hras1_all_4_cf[0] and all_4_perc_row[1] <= hras1_all_4_cf[1]:
                                params_fit_hras1_all_4_only.append(params_row)
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
                                                    if all_4_perc_row[1] >= hras1_all_4_cf[0] and all_4_perc_row[1] <= hras1_all_4_cf[1]:
                                                        params_fit_hras1_all_4.append(params_row)

                            if all_4_perc_row[1] >= hras0_all_4_cf[0] and all_4_perc_row[1] <= hras0_all_4_cf[1]:
                                params_fit_hras0_all_4_only.append(params_row)
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
                                                    if all_4_perc_row[1] >= hras0_all_4_cf[0] and all_4_perc_row[1] <= hras0_all_4_cf[1]:
                                                        params_fit_hras0_all_4.append(params_row)
                            break
                    if not completed_found:
                        params_failed.append(params_row)

    df = pd.DataFrame(all_rows, columns=outcols_full)
    df.to_csv(f'{outdir}/s_{round(s, 2)}_{round(s1, 2)}.csv', index=False)



df = pd.DataFrame(rows_match_braf_sd, columns=outcols_full)
df.to_csv(f'{outdir}/matching_braf_sd_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_braf, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_braf_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_braf_all_4, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_braf_all_4_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_braf_all_4_only, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_braf_all_4_only_cf_{CF}.csv', index=False)

df = pd.DataFrame(rows_match_hras1_sd, columns=outcols_full)
df.to_csv(f'{outdir}/matching_hras1_sd_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_hras1, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_hras1_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_hras1_all_4, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_hras1_all_4_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_hras1_all_4_only, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_hras1_all_4_only_cf_{CF}.csv', index=False)

df = pd.DataFrame(rows_match_hras0_sd, columns=outcols_full)
df.to_csv(f'{outdir}/matching_hras0_sd_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_hras0, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_hras0_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_hras0_all_4, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_hras0_all_4_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_hras0_all_4_only, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_hras0_all_4_only_cf_{CF}.csv', index=False)

df = pd.DataFrame(rows_match_egfr_sd, columns=outcols_full)
df.to_csv(f'{outdir}/matching_egfr_sd_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_egfr, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_egfr_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_egfr_all_4, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_egfr_all_4_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_egfr_all_4_only, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_egfr_all_4_only_cf_{CF}.csv', index=False)

df = pd.DataFrame(rows_match_egfr_shape_sd, columns=outcols_full)
df.to_csv(f'{outdir}/matching_egfr_shape_sd_cf_{CF}.csv', index=False)
df = pd.DataFrame(params_fit_egfr_shape, columns=outcols_params)
df.to_csv(f'{outdir}/params_fit_egfr_shape_cf_{CF}.csv', index=False)

df = pd.DataFrame(params_failed, columns=outcols_params)
df.to_csv(f'{outdir}/params_failed_runs.csv', index=False)
df = pd.DataFrame(params_not_run, columns=outcols_params)
df.to_csv(f'{outdir}/params_missing.csv', index=False)















