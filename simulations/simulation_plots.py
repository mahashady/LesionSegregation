import numpy as np
import pandas as pd
import matplotlib as mplt
import matplotlib.pyplot as plt
import os

custom_rcParams = {'font.family': 'Arial',
                   'font.size': 14,
                   'pdf.fonttype': 42,
                   'ps.fonttype': 42,
                   'svg.fonttype': 'none'
                   }


def get_mrca_perc(mrca_list):
    mrca = np.array(mrca_list).astype(int)
    mrca_val, counts = np.unique(mrca, return_counts=True)
    perc_to_annotate = [0 for i in range(6)]
    for m, c in zip(mrca_val, counts):
        perc_to_annotate[m] = c / sum(counts)
    return perc_to_annotate

def get_all_4_perc(mrca_list, all_4_list):
    out_list = []
    mrca_all_4_df = pd.DataFrame({'mrca': mrca_list, 'all_4': all_4_list})
    for mrca in range(4):
        subdf = mrca_all_4_df.loc[mrca_all_4_df['mrca'] == mrca]
        perc = sum(subdf['all_4'])/subdf.shape[0]
        out_list.append(perc)
    return out_list



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


def get_mrca_all_4(clones, mrca):
    two_gens_after = set()
    for clone in clones:
        two_after = clone[mrca:mrca+2]
        two_gens_after.add(two_after)
    return int(len(two_gens_after) == 4)


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
        # norm_factor = sum([i[1] for i in current_id_tuples])
        # current_id_tuples = [(id_, p / norm_factor) for id_, p in current_id_tuples]
    return mrca

def mrca_recursive_filter_all_4(clones, proportions):
    found = False
    current_id_tuples = [(id_, p) for id_, p in zip(clones, proportions)]
    threshold = 0.1
    counter = 0
    while not found:
        counter += 1
        current_clones = [i[0] for i in current_id_tuples]
        mrca = get_mrca_from_ids(current_clones)
        all_4 = get_mrca_all_4(current_clones, mrca)
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
        # norm_factor = sum([i[1] for i in current_id_tuples])
        # current_id_tuples = [(id_, p / norm_factor) for id_, p in current_id_tuples]
    return mrca, all_4


def initialize_background(low, high, color, title):
    x = np.arange(6)
    fig, ax = plt.subplots()
    ax.plot(x, low, '.', color=color, alpha=0.1)
    ax.plot(x, high, '.', color=color, alpha=0.1)
    ax.fill_between(x, low, high, color=color, alpha=0.2)
    ax.set_xlabel('MRCA generation')
    ax.set_ylabel('Frequency')
    ax.set_xticks(range(6))
    ax.set_xticklabels(['0', '1', '2', '3', '4', '5+'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(title)
    return fig, ax


def get_mrca_list_old(clones_df):
    tumor_size_min = 1000000
    mrca_filtered = []
    sim_groups = clones_df.groupby('simulation')
    for _, sim in sim_groups:
        tumor_size = sum(sim['clone size'].values)
        if tumor_size >= tumor_size_min:
            sim_filtered = sim.loc[sim['clone proportion'] >= 0.1]
            if sim_filtered.shape[0] > 0:
                mrca = get_mrca_from_ids(sim_filtered['clone id'].values)
                mrca_filtered.append(mrca)
    return mrca_filtered


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

def get_mrca_and_all_4_list(clones_df):
    tumor_size_min = 1000000
    mrca_filtered = []
    all_4_list = []
    sim_groups = clones_df.groupby('simulation')
    for _, sim in sim_groups:
        tumor_size = sum(sim['clone size'].values)
        if tumor_size >= tumor_size_min:
            mrca, all_4 = mrca_recursive_filter_all_4(sim['clone id'].values, sim['clone proportion'].values)
            mrca_filtered.append(mrca)
            all_4_list.append(all_4)
    return mrca_filtered, all_4_list


def add_simulation_line_plot(fig, ax, sim_dir, color, name):
    x = np.arange(6)
    if os.path.isdir(sim_dir):
        # print(sim_dir)
        if 'combined_clone_info_no_filter.csv' in os.listdir(sim_dir):
            print(sim_dir)
            # mrca_filtered = []
            clones_df = pd.read_csv(f'{sim_dir}/combined_clone_info_no_filter.csv')
            mrca_filtered = get_mrca_list(clones_df)

            mrca_row = get_mrca_perc(mrca_filtered)
            if mrca_filtered:
                print('plotting:', mrca_row)
            ax.plot(x, mrca_row, 'o--', color=color, label=name)
    return fig, ax

