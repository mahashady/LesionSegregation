import pandas as pd
import numpy as np
import os
from simulation_plots import get_mrca_perc, get_mrca_and_all_4_list


conf_interval_df = pd.read_csv('simulations/C3H_confidence_intervals_poisson.csv')
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

hras1_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(hras1_low, hras1_high)]
hras0_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(hras0_low, hras0_high)]



outdir = 'assumptions_test/'
if not os.path.isdir(outdir):
    os.mkdir(outdir)
cell_path = 'cells/grcm38.p6_rounded_driver_loci_freq_5/2025-10-10_11-22/'

s_range = np.arange(0.2, 0.9, 0.1)
parameters_to_test = [[0.15, 0.2, 0.05], [0.15, 0.2, 0.15],
                      [0.15, 0.3, 0.05], [0.15, 0.3, 0.15],
                      [0.25, 0.2, 0.05], [0.25, 0.2, 0.15],
                      [0.25, 0.3, 0.05], [0.25, 0.3, 0.15],
                      [0.15, 0.4, 0.25], [0.15, 0.4, 0.25],
                      [0.15, 0.5, 0.35], [0.15, 0.5, 0.35],
                      [0.25, 0.4, 0.25], [0.25, 0.4, 0.35],
                      [0.25, 0.5, 0.25], [0.25, 0.5, 0.35]]
hras1_fit = [[0 for i in parameters_to_test] for j in range(3)] # lm fitness, lm no fitness, all fitness
hras0_fit = [[0 for i in parameters_to_test] for j in range(3)]
missing = []

hras1_lm_fitness = []
hras1_lm_no_fitness = []
hras1_all_fitness = []
hras1_lm_fitness_no_fitness_comp = [0, 0] # match, no match
hras1_lm_fitness_all_fitness_comp = [0, 0]

hras0_lm_fitness = []
hras0_lm_no_fitness = []
hras0_all_fitness = []
hras0_lm_fitness_no_fitness_comp = [0, 0] # match, no match
hras0_lm_fitness_all_fitness_comp = [0, 0]

def check(conf_intervals, measured_mrca_list):
    if measured_mrca_list[0] >= conf_intervals[0][0] and measured_mrca_list[0] <= conf_intervals[0][1]:
        if measured_mrca_list[1] >= conf_intervals[1][0] and measured_mrca_list[1] <= conf_intervals[1][1]:
            if measured_mrca_list[2] >= conf_intervals[2][0] and measured_mrca_list[2] <= conf_intervals[2][1]:
                if measured_mrca_list[3] >= conf_intervals[3][0] and measured_mrca_list[3] <= conf_intervals[3][1]:
                    if measured_mrca_list[4] >= conf_intervals[4][0] and measured_mrca_list[4] <= conf_intervals[4][1]:
                        if measured_mrca_list[-1] >= conf_intervals[5][0] and measured_mrca_list[-1] <= conf_intervals[5][1]:
                            return 1
    return 0

for s in s_range:
    lesions_path = f'{cell_path}/lesions_1_drivers_0_other_s_{round(s, 1)}/simulation_results/'
    for param_ind in range(len(parameters_to_test)):
        params = parameters_to_test[param_ind]
        e = params[0]
        r = params[1]
        u = params[2]
        params_row = [round(s, 1), round(e, 2), round(r, 1), round(u, 2), round(r-u, 2)]
        run = True
        lm_fitness_path = f'{lesions_path}/epistasis_add_driver_on_nontx_strands_efr_{round(e, 2)}_eft_0.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
        if not os.path.isdir(lm_fitness_path):
            print(f'lm fitness looking for tx: {lm_fitness_path}')
            lm_fitness_path = f'{lesions_path}/epistasis_add_driver_on_tx_strands_efr_{round(e, 2)}_eft_0.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
            if not os.path.isdir(lm_fitness_path):
                print(f'tx also missing')
                run = False
                missing.append(['lm_fitness'] + params_row)
        if run:
            subdirs = sorted(os.listdir(f'{lm_fitness_path}'), reverse=True)
            completed_found = False
            for d in subdirs:
                subdir = f'{lm_fitness_path}/{d}/'
                if os.path.isdir(subdir) and 'combined_clone_info_no_filter.csv' in os.listdir(subdir):
                    completed_found = True
                    clones_df = pd.read_csv(f'{subdir}/combined_clone_info_no_filter.csv')
                    lm_fitness_mrca_list, all_4_list = get_mrca_and_all_4_list(clones_df)
                    lm_fitness_mrca_row = get_mrca_perc(lm_fitness_mrca_list)
                    lm_fitness_hras1_fit = check(hras1_range_sd, lm_fitness_mrca_row)
                    # hras1_fit[0][param_ind] = lm_fitness_hras1_fit
                    hras1_lm_fitness.append(params_row + [lm_fitness_hras1_fit])
                    lm_fitness_hras0_fit = check(hras0_range_sd, lm_fitness_mrca_row)
                    # hras0_fit[0][param_ind] = lm_fitness_hras0_fit
                    hras0_lm_fitness.append(params_row + [lm_fitness_hras0_fit])
                    break
            if not completed_found:
                print(f'no completed found {lm_fitness_path}')
                missing.append(['lm_fitness'] + params_row)
        run = True
        lm_no_fitness_path = f'{lesions_path}/epistasis_add_driver_on_nontx_strands_efr_{round(e, 2)}_eft_1.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
        if not os.path.isdir(lm_no_fitness_path):
            print(f'lm no fitness looking for tx {lm_no_fitness_path}')
            lm_no_fitness_path = f'{lesions_path}/epistasis_add_driver_on_tx_strands_efr_{round(e, 2)}_eft_1.0_2strands_True_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
            if not os.path.isdir(lm_no_fitness_path):
                print('tx also missing')
                run = False
                missing.append(['lm_no_fitness'] + params_row)
        if run:
            subdirs = sorted(os.listdir(f'{lm_no_fitness_path}'), reverse=True)
            completed_found = False
            for d in subdirs:
                subdir = f'{lm_no_fitness_path}/{d}/'
                if os.path.isdir(subdir) and 'combined_clone_info_no_filter.csv' in os.listdir(subdir):
                    completed_found = True
                    clones_df = pd.read_csv(f'{subdir}/combined_clone_info_no_filter.csv')
                    lm_no_fitness_mrca_list, all_4_list = get_mrca_and_all_4_list(clones_df)
                    lm_no_fitness_mrca_row = get_mrca_perc(lm_no_fitness_mrca_list)
                    lm_no_fitness_hras1_fit = check(hras1_range_sd, lm_no_fitness_mrca_row)
                    # hras1_fit[1][param_ind] = lm_no_fitness_hras1_fit
                    hras1_lm_no_fitness.append(params_row + [lm_no_fitness_hras1_fit])
                    lm_no_fitness_hras0_fit = check(hras0_range_sd, lm_no_fitness_mrca_row)
                    hras0_lm_no_fitness.append(params_row + [lm_no_fitness_hras0_fit])
                    break
            if not completed_found:
                print(f'no completed found {lm_no_fitness_path}')
                missing.append(['lm_no_fitness'] + params_row)


        run = True
        all_fitness_path = f'{lesions_path}/epistasis_add_driver_on_nontx_strands_efr_{round(e, 2)}_eft_0.0_2strands_False_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
        if not os.path.isdir(all_fitness_path):
            print(f'all fitness looking for tx {all_fitness_path}')
            all_fitness_path = f'{lesions_path}/epistasis_add_driver_on_tx_strands_efr_{round(e, 2)}_eft_0.0_2strands_False_repair_{round(r, 1)}_{round(u, 2)}_lesion_mut/'
            if not os.path.isdir(all_fitness_path):
                print('tx also missing')
                run = False
                missing.append(['all_fitness'] + params_row)
        if run:
            subdirs = sorted(os.listdir(f'{all_fitness_path}'), reverse=True)
            completed_found = False
            for d in subdirs:
                subdir = f'{all_fitness_path}/{d}/'
                if os.path.isdir(subdir) and 'combined_clone_info_no_filter.csv' in os.listdir(subdir):
                    completed_found = True
                    clones_df = pd.read_csv(f'{subdir}/combined_clone_info_no_filter.csv')
                    all_fitness_mrca_list, all_4_list = get_mrca_and_all_4_list(clones_df)
                    all_fitness_mrca_row = get_mrca_perc(all_fitness_mrca_list)
                    all_fitness_hras1_fit = check(hras1_range_sd, all_fitness_mrca_row)
                    # hras1_fit[2][param_ind] = all_fitness_hras1_fit
                    hras1_all_fitness.append(params_row + [all_fitness_hras1_fit])
                    all_fitness_hras0_fit = check(hras0_range_sd, all_fitness_mrca_row)
                    # hras0_fit[2][param_ind] = all_fitness_hras0_fit
                    hras0_all_fitness.append(params_row + [all_fitness_hras0_fit])
                    break
            if not completed_found:
                print(f'no completed found {all_fitness_path}')
                missing.append(['all_fitness'] + params_row)
        if lm_fitness_hras1_fit == lm_no_fitness_hras1_fit: hras1_lm_fitness_no_fitness_comp[0] += 1
        else: hras1_lm_fitness_no_fitness_comp[1] += 1
        if lm_fitness_hras1_fit == all_fitness_hras1_fit: hras1_lm_fitness_all_fitness_comp[0] += 1
        else: hras1_lm_fitness_all_fitness_comp[1] += 1

        if lm_fitness_hras0_fit == lm_no_fitness_hras0_fit: hras0_lm_fitness_no_fitness_comp[0] += 1
        else: hras0_lm_fitness_no_fitness_comp[1] += 1
        if lm_fitness_hras0_fit == all_fitness_hras0_fit: hras0_lm_fitness_all_fitness_comp[0] += 1
        else: hras0_lm_fitness_all_fitness_comp[1] += 1


outcols_fit = ['s', 'e', 'r', 'u', 'r-u', 'fit']
hras1_lm_fitness = pd.DataFrame(hras1_lm_fitness, columns=outcols_fit)
hras1_lm_fitness.to_csv(f'{outdir}/hras1_lm_fitness.csv')
hras1_lm_no_fitness = pd.DataFrame(hras1_lm_no_fitness, columns=outcols_fit)
hras1_lm_no_fitness.to_csv(f'{outdir}/hras1_lm_no_fitness.csv')
hras1_all_fitness = pd.DataFrame(hras1_all_fitness, columns=outcols_fit)
hras1_all_fitness.to_csv(f'{outdir}/hras1_all_fitness.csv')

hras0_lm_fitness = pd.DataFrame(hras0_lm_fitness, columns=outcols_fit)
hras0_lm_fitness.to_csv(f'{outdir}/hras0_lm_fitness.csv')
hras0_lm_no_fitness = pd.DataFrame(hras0_lm_no_fitness, columns=outcols_fit)
hras0_lm_no_fitness.to_csv(f'{outdir}/hras0_lm_no_fitness.csv')
hras0_all_fitness = pd.DataFrame(hras0_all_fitness, columns=outcols_fit)
hras0_all_fitness.to_csv(f'{outdir}/hras0_all_fitness.csv')



outcols_params = ['assumption', 's', 'e', 'r', 'u', 'r-u']
df = pd.DataFrame(missing, columns=outcols_params)
df.to_csv(f'{outdir}/params_missing.csv', index=False)

fit_array = np.array([['', 'match', "don't match"],
                      ['hras1 lm fitness - no fitness comparison', hras1_lm_fitness_no_fitness_comp[0], hras1_lm_fitness_no_fitness_comp[1]],
                      ['hras1 lm fitness - all fitness comparison', hras1_lm_fitness_all_fitness_comp[0], hras1_lm_fitness_all_fitness_comp[1]],
                      ['hras0 lm fitness - no fitness comparison', hras0_lm_fitness_no_fitness_comp[0], hras0_lm_fitness_no_fitness_comp[1]],
                      ['hras0 lm fitness - all fitness comparison', hras0_lm_fitness_all_fitness_comp[0], hras0_lm_fitness_all_fitness_comp[1]]
                      ])
np.savetxt(f'{outdir}/fit_array.txt', fit_array, fmt='%s')




