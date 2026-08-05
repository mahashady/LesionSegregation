# Goal: assess what fits under different assumptions, compare fits
# Compare:
#   1) baseline fitness 0.02 vs no baseline fitness (all drivers)
#   2) lw has fitness vs not (HRAS)

# Recording mechanism:

# Counts: two dataframes:
#   1) driver, s, e, r, u, fit original, fit 0.02 baseline
#   2) driver, s, e, r, u, fit original, fit lw fitness (should only include HRAS and fewer parameters)
##  Options for fit columns are 0/1/NA
##  When finished:
##      Group by driver
##      Count of each combination of fits (1,1), (1,0), (0,1), (0,0). Will have 4 counts per driver.
##      Check that totals reflect expected number of parameters sets testes, and is the same across drivers

# MSE: two lists:
#   1) Squared error per parameter set per driver for original vs 0.02 baseline fitness
#   2) Squared error per parameter set per driver for original vs lw has fitness
##  When finished
##      Check that first list matches parameter set counts * 4 from counts analysis
##      Check that second list matches parameter set counts * 2 from counts analysis
##      Compute overall mean for each list (mean across parameters and across drivers)

# Missing params dataframe: s, e, r, u, type missing

import pandas as pd
import numpy as np
import os
from simulation_plots import get_mrca_perc, get_mrca_and_all_4_list

# Steps:
##  Define ranges for each driver
conf_interval_df = pd.read_csv('C3H_confidence_intervals_poisson.csv')
conf_interval_groups = conf_interval_df.groupby('Gene_name')
for gene, grp in conf_interval_groups:
    srted_grp = grp.sort_values(by='division')
    grp_total = grp['N_samples'].astype(float).sum()
    CF = 95
    if gene == 'Hras1':
        hras1_low = srted_grp[f'Lower_CF_{CF}'].astype(float).values / grp_total
        hras1_high = srted_grp[f'Upper_CF_{CF}'].astype(float).values / grp_total
    elif gene == 'Hras0':
        hras0_low = srted_grp[f'Lower_CF_{CF}'].astype(float).values / grp_total
        hras0_high = srted_grp[f'Upper_CF_{CF}'].astype(float).values / grp_total
    elif gene == 'Braf':
        braf_low = srted_grp[f'Lower_CF_{CF}'].astype(float).values / grp_total
        braf_high = srted_grp[f'Upper_CF_{CF}'].astype(float).values / grp_total
    elif gene == 'Egfr':
        egfr_low = srted_grp[f'Lower_CF_{CF}'].astype(float).values / grp_total
        egfr_high = srted_grp[f'Upper_CF_{CF}'].astype(float).values / grp_total

hras1_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(hras1_low, hras1_high)]
hras0_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(hras0_low, hras0_high)]
braf_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(braf_low, braf_high)]
egfr_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(egfr_low, egfr_high)]

outdate = '_06132026'
outdir = f'assumptions_test{outdate}/'
if not os.path.isdir(outdir):
    os.mkdir(outdir)

cell_path = 'cells/grcm38.p6_rounded_driver_loci_freq_5/2025-10-10_11-22/'

s_range = np.arange(0.2, 0.9, 0.1)
parameters_to_test = [[0.15, 0.2, 0.05], [0.15, 0.2, 0.15],  #e, r, u
                      [0.15, 0.3, 0.05], [0.15, 0.3, 0.15], [0.15, 0.3, 0.25],
                      [0.15, 0.4, 0.05], [0.15, 0.4, 0.15], [0.15, 0.4, 0.25], [0.15, 0.4, 0.35],
                      [0.15, 0.5, 0.05], [0.15, 0.5, 0.15], [0.15, 0.5, 0.25], [0.15, 0.5, 0.35],
                      [0.25, 0.2, 0.05], [0.25, 0.2, 0.15],  # e, r, u
                      [0.25, 0.3, 0.05], [0.25, 0.3, 0.15], [0.25, 0.3, 0.25],
                      [0.25, 0.4, 0.05], [0.25, 0.4, 0.15], [0.25, 0.4, 0.25], [0.25, 0.4, 0.35],
                      [0.25, 0.5, 0.05], [0.25, 0.5, 0.15], [0.25, 0.5, 0.25], [0.25, 0.5, 0.35],
                      [0.35, 0.2, 0.05], [0.35, 0.3, 0.05], [0.35, 0.4, 0.05], [0.35, 0.5, 0.05],
                      [0.45, 0.2, 0.05], [0.45, 0.3, 0.05], [0.45, 0.4, 0.05], [0.45, 0.5, 0.05]]

baseline_fitness_comp = [] #driver, s, e, r, u, fit original, fit 0.02 baseline
baseline_05_comp = []
baseline_10_comp = []
baseline_15_comp = []
lw_fitness_comp = [] #driver, s, e, r, u, fit original, fit lw fitness
lm_fitness_comp = []

baseline_fitness_mse_list = []
baseline_05_mse_list = []
baseline_10_mse_list = []
baseline_15_mse_list = []
lw_fitness_mse_list = []
lm_fitness_mse_list = []

missing = [] #type missing (original, baseline fitness, lw fitness), s, e, r, u

def check(conf_intervals, measured_mrca_list):
    if measured_mrca_list[0] >= conf_intervals[0][0] and measured_mrca_list[0] <= conf_intervals[0][1]:
        if measured_mrca_list[1] >= conf_intervals[1][0] and measured_mrca_list[1] <= conf_intervals[1][1]:
            if measured_mrca_list[2] >= conf_intervals[2][0] and measured_mrca_list[2] <= conf_intervals[2][1]:
                if measured_mrca_list[3] >= conf_intervals[3][0] and measured_mrca_list[3] <= conf_intervals[3][1]:
                    if measured_mrca_list[4] >= conf_intervals[4][0] and measured_mrca_list[4] <= conf_intervals[4][1]:
                        if measured_mrca_list[-1] >= conf_intervals[5][0] and measured_mrca_list[-1] <= conf_intervals[5][1]:
                            return 1
    return 0
##  For each set of parameters:
for s in s_range:
    lesions_path = f'{cell_path}/lesions_1_drivers_0_other_s_{round(s, 1)}/simulation_results/'
    for param_ind in range(len(parameters_to_test)):
        params = parameters_to_test[param_ind]
        e = params[0]
        r = params[1]
        u = params[2]
        params_row = [round(s, 1), round(e, 2), round(r, 1), round(u, 2), round(r-u, 2)]
###     Does it fit original assumptions for each driver?
        orig_run = True
        orig_path = f'{lesions_path}/epistasis_add_driver_on_nontx_strands_efr_{round(e, 2)}_eft_0.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
        if not os.path.isdir(orig_path):
            print(f'lm fitness looking for tx: {orig_path}')
            orig_path = f'{lesions_path}/epistasis_add_driver_on_tx_strands_efr_{round(e, 2)}_eft_0.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
            if not os.path.isdir(orig_path):
                print(f'tx also missing')
                orig_run = False
                missing.append(['original'] + params_row)
                orig_hras0_fit, orig_hras1_fit, orig_braf_fit, orig_egfr_fit = None, None, None, None
                # if original is missing, no comparison possible, go to parameter set
                continue
        if orig_run:
            subdirs = sorted(os.listdir(f'{orig_path}'), reverse=True)
            orig_completed_found = False
            for d in subdirs:
                subdir = f'{orig_path}/{d}/'
                if os.path.isdir(subdir) and 'combined_clone_info_no_filter.csv' in os.listdir(subdir):
                    orig_completed_found = True
                    clones_df = pd.read_csv(f'{subdir}/combined_clone_info_no_filter.csv')
                    orig_mrca_list, _ = get_mrca_and_all_4_list(clones_df)
                    orig_mrca_row = get_mrca_perc(orig_mrca_list)
                    orig_hras0_fit = check(hras0_range_sd, orig_mrca_row)
                    orig_hras1_fit = check(hras1_range_sd, orig_mrca_row)
                    orig_braf_fit = check(braf_range_sd, orig_mrca_row)
                    orig_egfr_fit = check(egfr_range_sd, orig_mrca_row)
                    break # no need to look at other subdirectories
            if not orig_completed_found:
                orig_hras0_fit, orig_hras1_fit, orig_braf_fit, orig_egfr_fit = None, None, None, None
                print(f'no completed found {orig_path}')
                missing.append(['original'] + params_row)
                # if original is missing, no comparison possible, go to parameter set
                continue
###     Does it fit 0.02 baseline fitness for each driver?
        baseline_run = True
        baseline_path = f'{lesions_path}/epistasis_add_driver_on_tx_strands_efr_{round(e, 2)}_eft_0.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut_baseline_fitness_0.02/'
        if not os.path.isdir(baseline_path):
            baseline_hras0_fit, baseline_hras1_fit, baseline_braf_fit, baseline_egfr_fit = None, None, None, None
            baseline_run = False
            missing.append(['baseline fitness'] + params_row)
        if baseline_run:
            subdirs = sorted(os.listdir(f'{baseline_path}'), reverse=True)
            baseline_completed_found = False
            for d in subdirs:
                subdir = f'{baseline_path}/{d}/'
                if os.path.isdir(subdir) and 'combined_clone_info_no_filter.csv' in os.listdir(subdir):
                    baseline_completed_found = True
                    clones_df = pd.read_csv(f'{subdir}/combined_clone_info_no_filter.csv')
                    baseline_mrca_list, _ = get_mrca_and_all_4_list(clones_df)
                    baseline_mrca_row = get_mrca_perc(baseline_mrca_list)
                    baseline_hras0_fit = check(hras0_range_sd, baseline_mrca_row)
                    baseline_hras1_fit = check(hras1_range_sd, baseline_mrca_row)
                    baseline_braf_fit = check(braf_range_sd, baseline_mrca_row)
                    baseline_egfr_fit = check(egfr_range_sd, baseline_mrca_row)
                    break
            if not baseline_completed_found:
                baseline_hras0_fit, baseline_hras1_fit, baseline_braf_fit, baseline_egfr_fit = None, None, None, None
                print(f'no completed found {baseline_path}')
                missing.append(['baseline fitness'] + params_row)
###     Does it fit 0.05 baseline fitness till P30 for each driver?
        baseline_05_run = True
        baseline_05_path = f'{lesions_path}/epistasis_add_driver_on_tx_strands_efr_{round(e, 2)}_eft_0.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut_baseline_fitness_0.05_P30/'
        if not os.path.isdir(baseline_05_path):
            baseline_05_hras0_fit, baseline_05_hras1_fit, baseline_05_braf_fit, baseline_05_egfr_fit = None, None, None, None
            baseline_05_run = False
            missing.append(['baseline 0.05 fitness'] + params_row)
        if baseline_05_run:
            subdirs = sorted(os.listdir(f'{baseline_05_path}'), reverse=True)
            baseline_05_completed_found = False
            for d in subdirs:
                subdir = f'{baseline_05_path}/{d}/'
                if os.path.isdir(subdir) and 'combined_clone_info_no_filter.csv' in os.listdir(subdir):
                    baseline_05_completed_found = True
                    clones_df = pd.read_csv(f'{subdir}/combined_clone_info_no_filter.csv')
                    baseline_05_mrca_list, _ = get_mrca_and_all_4_list(clones_df)
                    baseline_05_mrca_row = get_mrca_perc(baseline_05_mrca_list)
                    baseline_05_hras0_fit = check(hras0_range_sd, baseline_05_mrca_row)
                    baseline_05_hras1_fit = check(hras1_range_sd, baseline_05_mrca_row)
                    baseline_05_braf_fit = check(braf_range_sd, baseline_05_mrca_row)
                    baseline_05_egfr_fit = check(egfr_range_sd, baseline_05_mrca_row)
                    break
            if not baseline_05_completed_found:
                baseline_05_hras0_fit, baseline_05_hras1_fit, baseline_05_braf_fit, baseline_05_egfr_fit = None, None, None, None
                print(f'no completed found {baseline_05_path}')
                missing.append(['baseline 0.05 fitness'] + params_row)
###     Does it fit 0.1 baseline fitness till P30 for each driver?
        baseline_10_run = True
        baseline_10_path = f'{lesions_path}/epistasis_add_driver_on_tx_strands_efr_{round(e, 2)}_eft_0.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut_baseline_fitness_0.1_P30/'
        if not os.path.isdir(baseline_10_path):
            baseline_10_hras0_fit, baseline_10_hras1_fit, baseline_10_braf_fit, baseline_10_egfr_fit = None, None, None, None
            baseline_10_run = False
            missing.append(['baseline 0.1 fitness'] + params_row)
        if baseline_10_run:
            subdirs = sorted(os.listdir(f'{baseline_10_path}'), reverse=True)
            baseline_10_completed_found = False
            for d in subdirs:
                subdir = f'{baseline_10_path}/{d}/'
                if os.path.isdir(subdir) and 'combined_clone_info_no_filter.csv' in os.listdir(subdir):
                    baseline_10_completed_found = True
                    clones_df = pd.read_csv(f'{subdir}/combined_clone_info_no_filter.csv')
                    baseline_10_mrca_list, _ = get_mrca_and_all_4_list(clones_df)
                    baseline_10_mrca_row = get_mrca_perc(baseline_10_mrca_list)
                    baseline_10_hras0_fit = check(hras0_range_sd, baseline_10_mrca_row)
                    baseline_10_hras1_fit = check(hras1_range_sd, baseline_10_mrca_row)
                    baseline_10_braf_fit = check(braf_range_sd, baseline_10_mrca_row)
                    baseline_10_egfr_fit = check(egfr_range_sd, baseline_10_mrca_row)
                    break
            if not baseline_10_completed_found:
                baseline_10_hras0_fit, baseline_10_hras1_fit, baseline_10_braf_fit, baseline_10_egfr_fit = None, None, None, None
                print(f'no completed found {baseline_10_path}')
                missing.append(['baseline 0.1 fitness'] + params_row)
###     Does it fit 0.15 baseline fitness till P30 for each driver?
        baseline_15_run = True
        baseline_15_path = f'{lesions_path}/epistasis_add_driver_on_tx_strands_efr_{round(e, 2)}_eft_0.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut_baseline_fitness_0.15_P30/'
        if not os.path.isdir(baseline_15_path):
            baseline_15_hras0_fit, baseline_15_hras1_fit, baseline_15_braf_fit, baseline_15_egfr_fit = None, None, None, None
            baseline_15_run = False
            missing.append(['baseline 0.15 fitness'] + params_row)
        if baseline_15_run:
            subdirs = sorted(os.listdir(f'{baseline_15_path}'), reverse=True)
            baseline_15_completed_found = False
            for d in subdirs:
                subdir = f'{baseline_15_path}/{d}/'
                if os.path.isdir(subdir) and 'combined_clone_info_no_filter.csv' in os.listdir(subdir):
                    baseline_15_completed_found = True
                    clones_df = pd.read_csv(f'{subdir}/combined_clone_info_no_filter.csv')
                    baseline_15_mrca_list, _ = get_mrca_and_all_4_list(clones_df)
                    baseline_15_mrca_row = get_mrca_perc(baseline_15_mrca_list)
                    baseline_15_hras0_fit = check(hras0_range_sd, baseline_15_mrca_row)
                    baseline_15_hras1_fit = check(hras1_range_sd, baseline_15_mrca_row)
                    baseline_15_braf_fit = check(braf_range_sd, baseline_15_mrca_row)
                    baseline_15_egfr_fit = check(egfr_range_sd, baseline_15_mrca_row)
                    break
            if not baseline_15_completed_found:
                baseline_15_hras0_fit, baseline_15_hras1_fit, baseline_15_braf_fit, baseline_15_egfr_fit = None, None, None, None
                print(f'no completed found {baseline_15_path}')
                missing.append(['baseline 0.15 fitness'] + params_row)
###     If applicable: Does it fit lw fitness for each driver?
        if e <= 0.25:
            lw_run = True
            lw_path = f'{lesions_path}/epistasis_add_driver_on_nontx_strands_efr_{round(e, 2)}_eft_0.0_2strands_False_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
            if not os.path.isdir(lw_path):
                print(f'all fitness looking for tx {lw_path}')
                lw_path = f'{lesions_path}/epistasis_add_driver_on_tx_strands_efr_{round(e, 2)}_eft_0.0_2strands_False_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
                if not os.path.isdir(lw_path):
                    print('tx also missing')
                    lw_path = f'{lesions_path}/epistasis_add_driver_on_tx_strands_efr_{round(e, 2)}_eft_0.0_2strands_False_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut_baseline_fitness_0.0/'
                    if not os.path.isdir(lw_path):
                        lw_hras0_fit, lw_hras1_fit = None, None
                        lw_run = False
                        missing.append(['lw fitness'] + params_row)
            if lw_run:
                subdirs = sorted(os.listdir(f'{lw_path}'), reverse=True)
                lw_completed_found = False
                for d in subdirs:
                    subdir = f'{lw_path}/{d}/'
                    if os.path.isdir(subdir) and 'combined_clone_info_no_filter.csv' in os.listdir(subdir):
                        lw_completed_found = True
                        clones_df = pd.read_csv(f'{subdir}/combined_clone_info_no_filter.csv')
                        lw_mrca_list, _ = get_mrca_and_all_4_list(clones_df)
                        lw_mrca_row = get_mrca_perc(lw_mrca_list)
                        lw_hras0_fit = check(hras0_range_sd, lw_mrca_row)
                        lw_hras1_fit = check(hras1_range_sd, lw_mrca_row)
                        break
                if not lw_completed_found:
                    lw_hras0_fit, lw_hras1_fit = None, None
                    print(f'no completed found {lw_path}')
                    missing.append(['lw fitness'] + params_row)

            lm_run = True
            lm_path = f'{lesions_path}/epistasis_add_driver_on_nontx_strands_efr_{round(e, 2)}_eft_1.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
            if not os.path.isdir(lm_path):
                print(f'only mm fitness looking for tx {lm_path}')
                lm_path = f'{lesions_path}/epistasis_add_driver_on_tx_strands_efr_{round(e, 2)}_eft_1.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
                if not os.path.isdir(lm_path):
                    print('tx also missing')
                    lm_path = f'{lesions_path}/epistasis_add_driver_on_tx_strands_efr_{round(e, 2)}_eft_1.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut_baseline_fitness_0.0/'
                    if not os.path.isdir(lm_path):
                        lm_hras0_fit, lm_hras1_fit = None, None
                        lm_run = False
                        missing.append(['lm fitness'] + params_row)
            if lm_run:
                subdirs = sorted(os.listdir(f'{lm_path}'), reverse=True)
                lm_completed_found = False
                for d in subdirs:
                    subdir = f'{lm_path}/{d}/'
                    if os.path.isdir(subdir) and 'combined_clone_info_no_filter.csv' in os.listdir(subdir):
                        lm_completed_found = True
                        clones_df = pd.read_csv(f'{subdir}/combined_clone_info_no_filter.csv')
                        lm_mrca_list, _ = get_mrca_and_all_4_list(clones_df)
                        lm_mrca_row = get_mrca_perc(lm_mrca_list)
                        lm_hras0_fit = check(hras0_range_sd, lm_mrca_row)
                        lm_hras1_fit = check(hras1_range_sd, lm_mrca_row)
                        break
                if not lm_completed_found:
                    lm_hras0_fit, lm_hras1_fit = None, None
                    print(f'no completed found {lm_path}')
                    missing.append(['lm fitness'] + params_row)
###     Is it missing from any of those?
###     If not:
###         Fill counts in dataframes accordingly
###         Compute each squared error and add to list
        if orig_hras0_fit == None:
            raise ValueError('Loop ran through without finding original assumptions simulation')
        if baseline_hras0_fit != None: # no need to check all of them
            baseline_fitness_comp.append(['Hras0', s, e, r, u, r-u, orig_hras0_fit, baseline_hras0_fit])
            baseline_fitness_comp.append(['Hras1', s, e, r, u, r - u, orig_hras1_fit, baseline_hras1_fit])
            baseline_fitness_comp.append(['Braf', s, e, r, u, r - u, orig_braf_fit, baseline_braf_fit])
            baseline_fitness_comp.append(['Egfr', s, e, r, u, r - u, orig_egfr_fit, baseline_egfr_fit])
            baseline_fitness_mse_list.append(sum([(x-y)**2 for x, y in zip(orig_mrca_row, baseline_mrca_row)]))
        if baseline_05_hras0_fit != None: # no need to check all of them
            baseline_05_comp.append(['Hras0', s, e, r, u, r-u, orig_hras0_fit, baseline_05_hras0_fit])
            baseline_05_comp.append(['Hras1', s, e, r, u, r - u, orig_hras1_fit, baseline_05_hras1_fit])
            baseline_05_comp.append(['Braf', s, e, r, u, r - u, orig_braf_fit, baseline_05_braf_fit])
            baseline_05_comp.append(['Egfr', s, e, r, u, r - u, orig_egfr_fit, baseline_05_egfr_fit])
            baseline_05_mse_list.append(sum([(x-y)**2 for x, y in zip(orig_mrca_row, baseline_05_mrca_row)]))
        if baseline_10_hras0_fit != None: # no need to check all of them
            baseline_10_comp.append(['Hras0', s, e, r, u, r-u, orig_hras0_fit, baseline_10_hras0_fit])
            baseline_10_comp.append(['Hras1', s, e, r, u, r - u, orig_hras1_fit, baseline_10_hras1_fit])
            baseline_10_comp.append(['Braf', s, e, r, u, r - u, orig_braf_fit, baseline_10_braf_fit])
            baseline_10_comp.append(['Egfr', s, e, r, u, r - u, orig_egfr_fit, baseline_10_egfr_fit])
            baseline_10_mse_list.append(sum([(x-y)**2 for x, y in zip(orig_mrca_row, baseline_10_mrca_row)]))
        if baseline_15_hras0_fit != None: # no need to check all of them
            baseline_15_comp.append(['Hras0', s, e, r, u, r-u, orig_hras0_fit, baseline_15_hras0_fit])
            baseline_15_comp.append(['Hras1', s, e, r, u, r - u, orig_hras1_fit, baseline_15_hras1_fit])
            baseline_15_comp.append(['Braf', s, e, r, u, r - u, orig_braf_fit, baseline_15_braf_fit])
            baseline_15_comp.append(['Egfr', s, e, r, u, r - u, orig_egfr_fit, baseline_15_egfr_fit])
            baseline_15_mse_list.append(sum([(x-y)**2 for x, y in zip(orig_mrca_row, baseline_15_mrca_row)]))
        if e <= 0.25:
            if lw_hras0_fit != None:
                lw_fitness_comp.append(['Hras0', s, e, r, u, r-u, orig_hras0_fit, lw_hras0_fit])
                lw_fitness_comp.append(['Hras1', s, e, r, u, r - u, orig_hras1_fit, lw_hras1_fit])
                lw_fitness_mse_list.append(sum([(x - y) ** 2 for x, y in zip(orig_mrca_row, lw_mrca_row)]))
            if lm_hras0_fit != None:
                lm_fitness_comp.append(['Hras0', s, e, r, u, r - u, orig_hras0_fit, lm_hras0_fit])
                lm_fitness_comp.append(['Hras1', s, e, r, u, r - u, orig_hras1_fit, lm_hras1_fit])
                lm_fitness_mse_list.append(sum([(x - y) ** 2 for x, y in zip(orig_mrca_row, lm_mrca_row)]))

### When finished, do final checks and record results.
##      Group counts df by driver
##      Count of each combination of fits (1,1), (1,0), (0,1), (0,0). Will have 4 counts per driver.
##      Check that totals reflect expected number of parameters sets testes, and is the same across drivers
baseline_fitness_comp_df = pd.DataFrame(baseline_fitness_comp, columns=['driver', 's', 'e', 'r', 'u', 'r-u', 'original fit', 'baseline fit'])
print('baseline fitness comparison df: ')
print(baseline_fitness_comp_df)
baseline_fitness_comp_counts = baseline_fitness_comp_df.groupby(by=['driver', 'original fit', 'baseline fit']).size().reset_index(name='count')
baseline_fitness_comp_df.to_csv(f'{outdir}/baseline_fitness_comparison.csv', index=False)
baseline_fitness_comp_counts.to_csv(f'{outdir}/baseline_fitness_comparison_counts.csv', index=False)
print(baseline_fitness_comp_counts.shape, 7*len(parameters_to_test)*4*4) # s counts * param counts * drivers * fit combos

baseline_05_comp_df = pd.DataFrame(baseline_05_comp, columns=['driver', 's', 'e', 'r', 'u', 'r-u', 'original fit', 'baseline 0.05 fit'])
print('baseline fitness comparison df: ')
print(baseline_05_comp_df)
baseline_05_comp_counts = baseline_05_comp_df.groupby(by=['driver', 'original fit', 'baseline 0.05 fit']).size().reset_index(name='count')
baseline_05_comp_df.to_csv(f'{outdir}/baseline_fitness_05_comparison.csv', index=False)
baseline_05_comp_counts.to_csv(f'{outdir}/baseline_fitness_05_comparison_counts.csv', index=False)

baseline_10_comp_df = pd.DataFrame(baseline_10_comp, columns=['driver', 's', 'e', 'r', 'u', 'r-u', 'original fit', 'baseline 0.1 fit'])
print('baseline fitness comparison df: ')
print(baseline_10_comp_df)
baseline_10_comp_counts = baseline_10_comp_df.groupby(by=['driver', 'original fit', 'baseline 0.1 fit']).size().reset_index(name='count')
baseline_10_comp_df.to_csv(f'{outdir}/baseline_fitness_10_comparison.csv', index=False)
baseline_10_comp_counts.to_csv(f'{outdir}/baseline_fitness_10_comparison_counts.csv', index=False)

baseline_15_comp_df = pd.DataFrame(baseline_15_comp, columns=['driver', 's', 'e', 'r', 'u', 'r-u', 'original fit', 'baseline 0.15 fit'])
print('baseline fitness comparison df: ')
print(baseline_15_comp_df)
baseline_15_comp_counts = baseline_15_comp_df.groupby(by=['driver', 'original fit', 'baseline 0.15 fit']).size().reset_index(name='count')
baseline_15_comp_df.to_csv(f'{outdir}/baseline_fitness_15_comparison.csv', index=False)
baseline_15_comp_counts.to_csv(f'{outdir}/baseline_fitness_15_comparison_counts.csv', index=False)

lw_fitness_comp_df = pd.DataFrame(lw_fitness_comp, columns=['driver', 's', 'e', 'r', 'u', 'r-u', 'original fit', 'lw fit'])
print('lw fitness comparison df: ')
print(lw_fitness_comp_df)
lw_fitness_comp_counts = lw_fitness_comp_df.groupby(by=['driver', 'original fit', 'lw fit']).size().reset_index(name='count')
lw_fitness_comp_df.to_csv(f'{outdir}/lw_fitness_comparison.csv', index=False)
lw_fitness_comp_counts.to_csv(f'{outdir}/lw_fitness_comparison_counts.csv', index=False)
print(lw_fitness_comp_counts.shape, 7*26*2*4) #s counts * param counts * drivers * fit combos

lm_fitness_comp_df = pd.DataFrame(lm_fitness_comp, columns=['driver', 's', 'e', 'r', 'u', 'r-u', 'original fit', 'lm no fitness fit'])
print('lm fitness comparison df: ')
print(lm_fitness_comp_df)
lm_fitness_comp_counts = lm_fitness_comp_df.groupby(by=['driver', 'original fit', 'lw no fitness fit']).size().reset_index(name='count')
lm_fitness_comp_df.to_csv(f'{outdir}/lm_fitness_comparison.csv', index=False)
lm_fitness_comp_counts.to_csv(f'{outdir}/lm_fitness_comparison_counts.csv', index=False)
print(lm_fitness_comp_counts.shape, 7*26*2*4) #s counts * param counts * drivers * fit combos

##      Compute overall mean for each squared error list (mean across parameters and across drivers)
baseline_fitness_mse = sum(baseline_fitness_mse_list)/len(baseline_fitness_mse_list)
baseline_05_mse = sum(baseline_05_mse_list)/len(baseline_05_mse_list)
baseline_10_mse = sum(baseline_10_mse_list)/len(baseline_10_mse_list)
baseline_15_mse = sum(baseline_15_mse_list)/len(baseline_15_mse_list)
lw_fitness_mse = sum(lw_fitness_mse_list)/len(lw_fitness_mse_list)
##      Check that first mse list matches parameter set counts from counts analysis
print(len(baseline_fitness_mse_list), baseline_fitness_comp_df.shape)
##      Check that second mselist matches parameter set counts from counts analysis
print(len(lw_fitness_mse_list), lw_fitness_comp_df.shape)
mse_output = [[baseline_fitness_mse, len(baseline_fitness_mse_list)],
            [baseline_05_mse, len(baseline_05_mse_list)],
            [baseline_10_mse, len(baseline_10_mse_list)],
            [baseline_15_mse, len(baseline_15_mse_list)],
              [lw_fitness_mse, len(lw_fitness_mse_list)]]
mse_output = pd.DataFrame(mse_output, columns=['MSE', 'count'], index=['baseline fitness comparison', 'baseline 5% fitness comparison',
                                                                       'baseline 10% fitness comparison', 'baseline 15% fitness comparison',
                                                                       'lw fitness comparison'])
mse_output.to_csv(f'{outdir}/mse_results.csv')

## write missing
outcols_params = ['assumption', 's', 'e', 'r', 'u', 'r-u']
df = pd.DataFrame(missing, columns=outcols_params)
df.to_csv(f'{outdir}/params_missing.csv', index=False)


