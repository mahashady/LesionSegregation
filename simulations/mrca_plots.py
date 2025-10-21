import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib import rcParams
custom_rcParams = {'font.family': 'Arial',
                   'font.size': 14,
                   'pdf.fonttype': 42,
                   'ps.fonttype': 42,
                   'svg.fonttype': 'none'
                   }

# update rcParams
rcParams.update(custom_rcParams)

title_maps = {'braf': 'BRAF V637E', 'egfr': 'EGFR F254I', 'hras0': 'HRAS Q61L', 'hras1': 'HRAS Q61R'}


def get_mrca_perc(mrca_list):
    mrca = np.array(mrca_list).astype(int)
    mrca_val, counts = np.unique(mrca, return_counts=True)
    perc_to_annotate = [0 for i in range(6)]
    for m, c in zip(mrca_val, counts):
        perc_to_annotate[m] = c / sum(counts)
    return perc_to_annotate

def get_mrca_from_ids(clones):
    if len(clones) == 1:
        return 5
    clone_ids = clones
    for gen in range(len(clone_ids[0])):
        first = clone_ids[0][gen]
        for other_clone in clone_ids[1:]:
            # print(first, other_clone[gen])
            if other_clone[gen] != first:
                return min(gen, 5)
    return min(gen, 5)

def mrca_recursive_filter(clones, proportions):
    found = False
    current_id_tuples = [(id_, p) for id_, p in zip(clones, proportions)]
    threshold = 0.1
    counter = 0
    while not found:
        counter += 1
        mrca = get_mrca_from_ids([i[0] for i in current_id_tuples])
        if mrca == 5:
            found = True
            break
        subclone1_prop = 0
        subclone2_prop = 0
        for id_, p in current_id_tuples:
            if id_[mrca] == 'A':
                subclone1_prop += p
            elif id_[mrca] == 'B':
                subclone2_prop += p
        if subclone1_prop < threshold:
            current_id_tuples = [(id_, p) for id_, p in current_id_tuples if id_[mrca] != 'A']
        elif subclone2_prop < threshold:
            current_id_tuples = [(id_, p) for id_, p in current_id_tuples if id_[mrca] != 'B']
        else:
            found = True
            break
    return mrca

def get_mrca_list(clones_df):
    tumor_size_min = 1000000
    mrca_filtered = []
    sim_groups = clones_df.groupby('simulation')
    for _, sim in sim_groups:
        tumor_size = sum(sim['clone size'].values)
        if tumor_size >= tumor_size_min:
            mrca = mrca_recursive_filter(sim['clone id'].values, sim['clone proportion'].values)
            mrca_filtered.append(mrca)
    return mrca_filtered

def initialize_background(low, high, color, title):
    x = np.arange(6)
    fig, ax = plt.subplots()
    ax.plot(x, low, '.', color=color, alpha=0.1)
    ax.plot(x, high, '.', color=color, alpha=0.1)
    ax.fill_between(x, low, high, color=color, alpha=0.2)
    ax.set_xlabel('MRCA generation', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.set_xticks(range(6))
    ax.set_xticklabels(['0', '1', '2', '3', '4', '5+'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(title)
    return fig, ax

def add_simulation_line_plot(fig, ax, sim_dir, color, name):
    x = np.arange(6)
    if os.path.isdir(sim_dir):
        # print(sim_dir)
        if 'combined_clone_info_no_filter.csv' in os.listdir(sim_dir):
            # print(sim_dir)
            # mrca_filtered = []
            clones_df = pd.read_csv(f'{sim_dir}/combined_clone_info_no_filter.csv')
            mrca_filtered = get_mrca_list(clones_df)
            mrca_row = get_mrca_perc(mrca_filtered)
            if mrca_filtered:
                print('plotting:', mrca_row)
            ax.plot(x, mrca_row, 'o--', color=color, label=name)
    return fig, ax

def main():
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

    colors = ['deepskyblue', 'lightcoral', 'dodgerblue', 'darkorange', 'steelblue', 'olivedrab', 'red', 'firebrick',
              'darkolivegreen', 'brown', 'lightseagreen', 'rosybrown', 'indianred']

    outdir = 'figures/final_mrca_figures/'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    cell_path = 'cells/grcm38.p6_rounded_autosomes_driver_loci_freq_5/2024-06-02_13-23/'

    lesion_20 = f'{cell_path}/lesions_1_drivers_0_other_s_0.2/simulation_results/'
    lesion_30 = f'{cell_path}/lesions_1_drivers_0_other_s_0.3/simulation_results/'
    lesion_40 = f'{cell_path}/lesions_1_drivers_0_other_s_0.4/simulation_results/'
    lesion_50 = f'{cell_path}/lesions_1_drivers_0_other_s_0.5/simulation_results/'
    lesion_60 = f'{cell_path}/lesions_1_drivers_0_other_s_0.6/simulation_results/'
    lesion_70 = f'{cell_path}/lesions_1_drivers_0_other_s_0.7/simulation_results/'
    lesion_80 = f'{cell_path}/lesions_1_drivers_0_other_s_0.8/simulation_results/'

    two_lesions_hi_lo_5 = f'{cell_path}/lesions_2_drivers_0_other_smax_0.5_s_0.3/simulation_results/'
    two_lesions_eq_5 = f'{cell_path}/lesions_2_drivers_0_other_smax_0.5_s_0.25/simulation_results/'
    three_lesions_5 = f'{cell_path}/lesions_3_drivers_0_other_smax_0.5_s_0.2/simulation_results/'
    four_lesions_5 = f'{cell_path}/lesions_4_drivers_0_other_smax_0.5_s_0.2/simulation_results/'
    five_lesions_5 = f'{cell_path}/lesions_5_drivers_0_other_smax_0.5_s_0.2/simulation_results/'

    two_lesions_eq_6 = f'{cell_path}/lesions_2_drivers_0_other_smax_0.6_s_0.25/simulation_results/'
    three_lesions_6 = f'{cell_path}/lesions_3_drivers_0_other_smax_0.6_s_0.2/simulation_results/'
    four_lesions_6 = f'{cell_path}/lesions_4_drivers_0_other_smax_0.6_s_0.2/simulation_results/'
    five_lesions_6 = f'{cell_path}/lesions_5_drivers_0_other_smax_0.6_s_0.2/simulation_results/'

    def plot_r15_u05_vary_s(background_low, background_high, background_color, background_name):
        fig, ax = initialize_background(background_low, background_high, background_color, background_name)

        tx_paths = [
            f'{lesion_20}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_11-08_rs_0_500k_s_0.2/',
            f'{lesion_30}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_15-34_rs_0_500k_s_0.3/',
            f'{lesion_40}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_17-12_rs_0_500k_s_0.4/',
            f'{lesion_50}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_11-07_rs_0_500k_s_0.5/',
            f'{lesion_60}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_15-36_rs_0_500k_s_0.6/',
            f'{lesion_70}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_17-04_rs_0_500k_s_0.7/',
            f'{lesion_80}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_11-06_rs_0_500k_s_0.8/']
        names = ['s=0.2', 's=0.3', 's=0.4', 's=0.5', 's=0.6', 's=0.7', 's=0.8']
        for i in range(len(tx_paths)):
            fig, ax = add_simulation_line_plot(fig, ax, tx_paths[i], colors[i], names[i])
        ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.05, 1))
        # ax.set_title('r = 0.15, u = 0.05')
        # ax.set_title(f'{title_maps[background_name]} e=0.15, r=0.15, u=0.05', fontsize=14)
        ax.set_title(f'e=0.15, r-u=0.15, u=0.05', fontsize=14)
        filename = f'{outdir}/{background_name}_vary_s_r_0.15_u_0.05_cf_{CF}'
        for suffix in ['png']: #, 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)


    def plot_r20_u15_vary_s(background_low, background_high, background_color, background_name):
        fig, ax = initialize_background(background_low, background_high, background_color, background_name)

        tx_paths = [
            f'{lesion_20}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.35_0.15_lesion_mut/2024-08-13_12-06_rs_0_500k//',
            f'{lesion_30}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.35_0.15_lesion_mut/2024-08-13_12-12_rs_0_500k/',
            f'{lesion_40}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.35_0.15_lesion_mut/2024-08-13_12-15_rs_0_500k//',
            f'{lesion_50}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.35_0.15_lesion_mut/2024-08-13_12-17_rs_0_500k/',
            f'{lesion_60}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.35_0.15_lesion_mut/2024-08-13_12-23_rs_0_500k/',
            f'{lesion_70}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.35_0.15_lesion_mut/2024-08-13_12-22_rs_0_500k/',
            f'{lesion_80}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.35_0.15_lesion_mut/2024-08-13_12-23_rs_0_500k/']
        names = ['s=0.2', 's=0.3', 's=0.4', 's=0.5', 's=0.6', 's=0.7', 's=0.8']
        for i in range(len(tx_paths)):
            fig, ax = add_simulation_line_plot(fig, ax, tx_paths[i], colors[i],
                                                           names[i])
        ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.05, 1))
        # ax.set_title(f'{title_maps[background_name]} e=0.15, r=0.2, u=0.15', fontsize=14)
        ax.set_title(f'e=0.15, r-u=0.2, u=0.15', fontsize=14)
        filename = f'{outdir}/{background_name}_vary_s_r_0.2_u_0.15_cf_{CF}'
        for suffix in ['png']: #, 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)

    def plot_r20_u05_vary_s(background_low, background_high, background_color, background_name):
        fig, ax = initialize_background(background_low, background_high, background_color, background_name)
        tx_paths = [
            f'{lesion_20}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.25_0.05_lesion_mut/2024-08-13_11-57_rs_0_500k//',
            f'{lesion_30}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.25_0.05_lesion_mut/2024-08-13_12-12_rs_0_500k/',
            f'{lesion_40}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.25_0.05_lesion_mut/2024-08-13_12-13_rs_0_500k//',
            f'{lesion_50}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.25_0.05_lesion_mut/2024-08-13_12-17_rs_0_500k//',
            f'{lesion_60}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.25_0.05_lesion_mut/2024-08-13_12-19_rs_0_500k//',
            f'{lesion_70}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.25_0.05_lesion_mut/2024-08-13_12-22_rs_0_500k/',
            f'{lesion_80}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.25_0.05_lesion_mut/2024-08-13_12-22_rs_0_500k/']
        names = ['s=0.2', 's=0.3', 's=0.4', 's=0.5', 's=0.6', 's=0.7', 's=0.8']
        for i in range(len(tx_paths)):
            fig, ax = add_simulation_line_plot(fig, ax, tx_paths[i], colors[i],
                                                           names[i])
        ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.05, 1))
        ax.set_title(f'e=0.15, r-u=0.2, u=0.05', fontsize=14)
        filename = f'{outdir}/{background_name}_vary_s_r_0.2_u_0.05_cf_{CF}'
        for suffix in ['png']: #, 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)


    def plot_vary_drivers_s60_r_15_u05_e_15(background_low, background_high, background_color, background_name):
        fig, ax = initialize_background(background_low, background_high, background_color, background_name)
        path_list = [
            f'{lesion_60}/epistasis_add_driver_on_nontx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-07-02_10-51_rs_0_500k_0.6//',
            f'{two_lesions_eq_6}/epistasis_syn_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-21_rs_0_100k//',
            f'{three_lesions_6}/epistasis_syn_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-22_rs_0_100k//',
            f'{four_lesions_6}/epistasis_syn_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-24_rs_0_100k//',
            f'{five_lesions_6}/epistasis_syn_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-42_rs_0_100k//']
        names = ['1 driver s=0.6', '2 drivers s=0.6', '3 drivers s=0.6', '4 drivers s=0.6',
                         '5 drivers s=0.6']
        for i in range(len(names)):
            fig, ax = add_simulation_line_plot(fig, ax, path_list[i], colors[i],
                                                           names[i])
        ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.05, 1))
        ax.set_title(f'e=0.15, r-u=0.15, u=0.05', fontsize=14)
        filename = f'{outdir}/{background_name}_vary_driver_count_s_0.6_syn_r_0.15_u_0.05_e_0.15_cf_{CF}'
        for suffix in ['png']: #, 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)


    def plot_vary_drivers_s60_r_15_u05_e_25(background_low, background_high, background_color, background_name):
        fig, ax = initialize_background(background_low, background_high, background_color, background_name)
        path_list = [
            f'{lesion_60}/epistasis_add_driver_on_nontx_strands_efr_0.25_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-07-02_10-53_rs_0_500k_0.6/',
            f'{two_lesions_eq_6}/epistasis_syn_driver_on_tx_strands_efr_0.25_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-18_rs_0_100k/',
            f'{three_lesions_6}/epistasis_syn_driver_on_tx_strands_efr_0.25_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-22_rs_0_100k//',
            f'{four_lesions_6}/epistasis_syn_driver_on_tx_strands_efr_0.25_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-24_rs_0_100k/',
            f'{five_lesions_6}/epistasis_syn_driver_on_tx_strands_efr_0.25_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-46_rs_0_100k/']
        names = ['1 driver s=0.6', '2 drivers s=0.6', '3 drivers s=0.6', '4 drivers s=0.6',
                         '5 drivers s=0.6']
        for i in range(len(names)):
            fig, ax = add_simulation_line_plot(fig, ax, path_list[i], colors[i],
                                                           names[i])
        ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.05, 1))
        ax.set_title(f'e=0.25, r-u=0.15, u=0.05', fontsize=14)
        filename = f'{outdir}/{background_name}_vary_driver_count_s_0.6_syn_r_0.15_u_0.05_e_0.25_cf_{CF}'
        for suffix in ['png']: #, 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)

    def plot_vary_drivers_s60_r_15_u05_e_35(background_low, background_high, background_color, background_name):
        fig, ax = initialize_background(background_low, background_high, background_color, background_name)
        path_list = [
            f'{lesion_60}/epistasis_add_driver_on_nontx_strands_efr_0.35_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-07-02_10-54_rs_0_500k_0.6/',
            f'{two_lesions_eq_6}/epistasis_syn_driver_on_tx_strands_efr_0.35_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-19_rs_0_100k/',
            f'{three_lesions_6}/epistasis_syn_driver_on_tx_strands_efr_0.35_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-21_rs_0_100k/',
            f'{four_lesions_6}/epistasis_syn_driver_on_tx_strands_efr_0.35_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-26_rs_0_100k/',
            f'{five_lesions_6}/epistasis_syn_driver_on_tx_strands_efr_0.35_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-48_rs_0_100k/']
        names = ['1 driver s=0.6', '2 drivers s=0.6', '3 drivers s=0.6', '4 drivers s=0.6',
                         '5 drivers s=0.6']
        for i in range(len(names)):
            fig, ax = add_simulation_line_plot(fig, ax, path_list[i], colors[i],
                                                           names[i])
        ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.05, 1))
        ax.set_title(f'e=0.35, r-u=0.15, u=0.05', fontsize=14)
        filename = f'{outdir}/{background_name}_vary_driver_count_s_0.6_syn_r_0.15_u_0.05_e_0.35_cf_{CF}'
        for suffix in ['png']: #, 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)

    def plot_vary_drivers_s60_r_15_u05_e_45(background_low, background_high, background_color, background_name):
        fig, ax = initialize_background(background_low, background_high, background_color, background_name)
        path_list = [
            f'{lesion_60}/epistasis_add_driver_on_nontx_strands_efr_0.45_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-07-02_10-55_rs_0_500k_0.6/',
            f'{two_lesions_eq_6}/epistasis_syn_driver_on_tx_strands_efr_0.45_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-19_rs_0_100k/',
            f'{three_lesions_6}/epistasis_syn_driver_on_tx_strands_efr_0.45_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-24_rs_0_100k/',
            f'{four_lesions_6}/epistasis_syn_driver_on_tx_strands_efr_0.45_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-38_rs_0_100k//',
            f'{five_lesions_6}/epistasis_syn_driver_on_tx_strands_efr_0.45_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-09-04_19-48_rs_0_100k//']
        names = ['1 driver s=0.6', '2 drivers s=0.6', '3 drivers s=0.6', '4 drivers s=0.6',
                         '5 drivers s=0.6']
        for i in range(len(names)):
            fig, ax = add_simulation_line_plot(fig, ax, path_list[i], colors[i],
                                                           names[i])
        ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.05, 1))
        ax.set_title(f'e=0.45, r-u=0.15, u=0.05', fontsize=14)
        filename = f'{outdir}/{background_name}_vary_driver_count_s_0.6_syn_r_0.15_u_0.05_e_0.45_cf_{CF}'
        for suffix in ['png']: #, 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)


    def plot_vary_drivers_s50_r_15_u05(background_low, background_high, background_color, background_name):
        fig, ax = initialize_background(background_low, background_high, background_color, background_name)
        path_list = [
            f'{lesion_50}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-08-13_12-17_rs_0_500k/',
            f'{two_lesions_hi_lo_5}/epistasis_syn_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-08-13_15-26_rs_0_100k//',
            f'{three_lesions_5}/epistasis_syn_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-08-27_11-36_rs_0_100k//',
            f'{four_lesions_5}/epistasis_syn_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-08-27_11-38_rs_0_100k//',
            f'{five_lesions_5}/epistasis_syn_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-08-27_11-42_rs_0_100k/']
        names = ['1 driver s=0.5', '2 drivers s=0.5', '3 drivers s=0.5', '4 drivers s=0.5',
                         '5 drivers s=0.5']
        for i in range(len(names)):
            fig, ax = add_simulation_line_plot(fig, ax, path_list[i], colors[i],
                                                           names[i])
        ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.05, 1))
        ax.set_title(f'e=0.15, r-u=0.15, u=0.05', fontsize=14)
        filename = f'{outdir}/{background_name}_vary_driver_count_s_0.5_syn_r_0.15_u_0.05_e_0.15_cf_{CF}'
        for suffix in ['png']: #, 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)

    def plot_vary_drivers_s50_r_15_u35_e15(background_low, background_high, background_color, background_name):
        fig, ax = initialize_background(background_low, background_high, background_color, background_name)
        path_list = [
            f'{lesion_50}/epistasis_add_driver_on_nontx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.5_0.35_lesion_mut/2024-07-02_09-16_rs_0_500k_0.5//',
            f'{two_lesions_eq_5}/epistasis_syn_driver_on_nontx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.5_0.35_lesion_mut/2025-02-25_11-44_rs_0_500k_smax_0.5_hi_lo/',
            f'{three_lesions_5}/epistasis_syn_driver_on_nontx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.5_0.35_lesion_mut/2025-04-28_11-49_rs_0_100k_s0.5/',
            f'{four_lesions_5}/epistasis_syn_driver_on_nontx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.5_0.35_lesion_mut/2025-04-28_11-49_rs_0_100k_s0.5/',
            f'{five_lesions_5}/epistasis_syn_driver_on_nontx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.5_0.35_lesion_mut/2025-04-28_11-49_rs_0_100k_s0.5/']
        names = ['1 driver s=0.5', '2 drivers s=0.5', '3 drivers s=0.5', '4 drivers s=0.5',
                         '5 drivers s=0.5']
        for i in range(len(names)):
            fig, ax = add_simulation_line_plot(fig, ax, path_list[i], colors[i],
                                                           names[i])
        ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.05, 1))
        ax.set_title(f'e=0.15, r-u=0.15, u=0.35', fontsize=14)
        filename = f'{outdir}/{background_name}_vary_driver_count_s_0.5_syn_r_0.15_u_0.35_e_0.15_cf_{CF}'
        for suffix in ['png']: #, 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)


    def plot_s_50_r15_vary_u(background_low, background_high, background_color, background_name):
        fig, ax = initialize_background(background_low, background_high, background_color, background_name)
        path_list = [
            f'{lesion_50}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_11-07_rs_0_500k_s_0.5/',
            f'{lesion_50}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.3_0.15_lesion_mut/2024-06-03_11-07_rs_0_500k_s_0.5/',
            f'{lesion_50}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.4_0.25_lesion_mut/2024-06-03_11-07_rs_0_500k_s_0.5/',
            f'{lesion_50}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.5_0.35_lesion_mut/2024-06-03_11-07_rs_0_500k_s_0.5/',
            f'{lesion_50}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.6_0.45_lesion_mut/2024-06-03_11-07_rs_0_500k_s_0.5/',
            f'{lesion_50}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.7_0.55_lesion_mut/2024-06-14_14-55_rs_0_500k/']
        names = ['u=0.05', 'u=0.15', 'u=0.25', 'u=0.35', 'u=0.45', 'u=0.55']
        for i in range(len(names)):
            fig, ax = add_simulation_line_plot(fig, ax, path_list[i],
                                               colors[i], names[i])
        ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.05, 1))
        # ax.set_title(f'{title_maps[background_name]} s=0.5, e=0.15, r-u=0.15', fontsize=14)
        ax.set_title(f's=0.5, e=0.15, r-u=0.15', fontsize=14)
        filename = f'{outdir}/{background_name}_vary_u_s0.5_r0.15_cf_{CF}'
        for suffix in ['png']:#, 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)

    # Combined braf hras plot (vary s, r 0.2, u 0.05, e0.15)
    x = np.arange(6)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    s_tx_paths = [
        f'{lesion_20}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_11-08_rs_0_500k_s_0.2/',
        f'{lesion_30}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_15-34_rs_0_500k_s_0.3/',
        f'{lesion_40}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_17-12_rs_0_500k_s_0.4/',
        f'{lesion_50}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_11-07_rs_0_500k_s_0.5/',
        f'{lesion_60}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_15-36_rs_0_500k_s_0.6/',
        f'{lesion_70}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_17-04_rs_0_500k_s_0.7/',
        f'{lesion_80}/epistasis_add_driver_on_tx_strands_efr_0.15_eft_0.0_2strands_True_repair_0.2_0.05_lesion_mut/2024-06-03_11-06_rs_0_500k_s_0.8/']
    names = ['s=0.2', 's=0.3', 's=0.4', 's=0.5', 's=0.6', 's=0.7', 's=0.8']
    # hras 1
    ax1.plot(x, hras1_low, '.', color='dodgerblue', alpha=0.1)
    ax1.plot(x, hras1_high, '.', color='dodgerblue', alpha=0.1)
    ax1.fill_between(x, hras1_low, hras1_high, color='dodgerblue', alpha=0.2)
    ax1.set_xlabel('MRCA generation')
    ax1.set_ylabel('Frequency')
    ax1.set_xticks(range(6))
    ax1.set_xticklabels(['0', '1', '2', '3', '4', '5+'])
    # ax1.set_ylim(0, 1)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    # braf
    ax2.plot(x, braf_low, '.', color='palevioletred', alpha=0.1)
    ax2.plot(x, braf_high, '.', color='palevioletred', alpha=0.1)
    ax2.fill_between(x, braf_low, braf_high, color='palevioletred', alpha=0.2)
    ax2.set_xlabel('MRCA generation')
    # ax2.set_yticks([])
    ax2.set_xticks(range(6))
    ax2.set_xticklabels(['0', '1', '2', '3', '4', '5+'])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    for i in range(len(s_tx_paths)):
        sim_dir = s_tx_paths[i]
        if os.path.isdir(sim_dir):
            # print(sim_dir)
            if 'combined_clone_info_no_filter.csv' in os.listdir(sim_dir):
                clones_df = pd.read_csv(f'{sim_dir}/combined_clone_info_no_filter.csv')
                mrca_filtered = get_mrca_list(clones_df)

                mrca_row = get_mrca_perc(mrca_filtered)
                if mrca_filtered:
                    print('plotting:', mrca_row)
                ax1.plot(x, mrca_row, 'o--', color=colors[i], label=names[i])
                ax2.plot(x, mrca_row, 'o--', color=colors[i], label=names[i])
    ax2.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.05, 1))
    # ax.set_title('r = 0.15, u = 0.05')
    ax1.set_title(f'{title_maps["hras1"]}')
    ax2.set_title(f'{title_maps["braf"]}')
    fig.suptitle('Parameters: e=0.15, r=0.15, u=0.05')
    filename = f'{outdir}/combined_hras1_braf_vary_s_r_0.15_u_0.05_cf_{CF}'
    for suffix in ['png']: #, 'svg', 'pdf']:
        outname = f'{filename}.{suffix}'
        plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)


    plot_r15_u05_vary_s(hras1_low, hras1_high, 'dodgerblue', 'hras1') # main figure, r=0.2, r-u=0.15

    plot_r15_u05_vary_s(braf_low, braf_high, 'palevioletred', 'braf')
    plot_s_50_r15_vary_u(hras0_low, hras0_high, 'dodgerblue', 'hras0')
    #
    plot_r20_u05_vary_s(hras1_low, hras1_high, 'dodgerblue', 'hras1') #supplement, r=0.25 (r-u=0.2)
    plot_r20_u05_vary_s(braf_low, braf_high, 'palevioletred', 'braf')
    plot_r20_u15_vary_s(hras0_low, hras0_high, 'dodgerblue', 'hras0') # r=0.35, r-u=0.2
    #
    plot_vary_drivers_s50_r_15_u05(hras1_low, hras1_high, 'dodgerblue', 'hras1')
    plot_vary_drivers_s50_r_15_u05(braf_low, braf_high, 'palevioletred', 'braf')
    plot_vary_drivers_s50_r_15_u35_e15(hras0_low, hras0_high, 'dodgerblue', 'hras0')
    #
    plot_vary_drivers_s60_r_15_u05_e_15(egfr_low, egfr_high, 'seagreen', 'egfr')
    plot_vary_drivers_s60_r_15_u05_e_25(egfr_low, egfr_high, 'seagreen', 'egfr')
    plot_vary_drivers_s60_r_15_u05_e_35(egfr_low, egfr_high, 'seagreen', 'egfr')
    plot_vary_drivers_s60_r_15_u05_e_45(egfr_low, egfr_high, 'seagreen', 'egfr')


if __name__ == '__main__':
    main()