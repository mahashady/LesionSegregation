# read file with min/max s
# get min s for each r, e pair (across u)
# record min s
# create r list, e list, s grid
# make plot
import pandas as pd
import numpy as np
import matplotlib as mplt
import matplotlib.pyplot as plt
import seaborn as sns
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

indir = os.path.expanduser('~/Desktop/drivers/analytics_results/')
outdir = os.path.expanduser('~/Desktop/drivers/figures/parameter_space_plots/')


CF = 95
files = [f'min_max_s_per_param_braf_cf_{CF}.csv', f'min_max_s_per_param_egfr_cf_{CF}.csv',
         f'min_max_s_per_param_hras0_cf_{CF}.csv', f'min_max_s_per_param_hras1_cf_{CF}.csv']
names = ['braf', 'egfr', 'hras0', 'hras1']
title_maps = {'braf': 'BRAF V637E', 'egfr': 'EGFR F254I', 'hras0': 'HRAS Q61L', 'hras1': 'HRAS Q61R'}

format = np.vectorize(lambda x: round(x, 2))
rounding = lambda x: round(x, 2)
for i in range(len(files)):
    df = pd.read_csv(f'{indir}/{files[i]}').dropna(subset='min s')
    df['r-u'] = df['r'] - df['u']
    df = df.loc[df['r-u'] >= 0]

    # r vs e plots
    r = np.arange(0, 1.05, 0.05)
    e = np.arange(0, 1.05, 0.05)
    s_min = []
    s_max = []
    anything_fits = False
    for e_val in e:
        e_row_min = []
        e_row_max = []
        e_val = round(e_val, 2)
        for r_val in r:
            r_val = round(r_val, 2)
            subdf = df.loc[(df['r'].apply(rounding) == r_val) & (df['e'].apply(rounding) == e_val)]
            print(e_val, r_val, round(df['e'], 2), round(df['r'], 2))
            print(subdf)
            print()
            if not subdf.empty:
                anything_fits = True
            if len(np.unique(subdf['r'])) > 1 or len(np.unique(subdf['e'])) > 1:
                print('TROUBLE!!')
            e_row_min.append(round(subdf['min s'].min(), 2))
            e_row_max.append(round(subdf['max s'].max(), 2))
        s_min.append(e_row_min)
        s_max.append(e_row_max)
    # reverse so e=0 is the bottom row
    s_min.reverse()
    s_max.reverse()
    # e.reverse()
    s_min_df = pd.DataFrame(s_min, index=format(np.flip(e)), columns=format(r))
    s_max_df = pd.DataFrame(s_max, index=format(np.flip(e)), columns=format(r))
    if anything_fits:
        fig, ax = plt.subplots()
        ax = sns.heatmap(s_min_df, cmap='YlGnBu', vmin=0, vmax=1)
        ax.set_xlabel('Repair Rate', fontsize=14)
        ax.set_ylabel('Error-free Replication Rate', fontsize=14)
        # ax.set_title(f'{title_maps[names[i]]}', fontsize=14)
        cbar = ax.collections[0].colorbar
        cbar.set_label('Minimum Selection Coefficient', fontsize=14)
        filename = f'{outdir}/heatmap_analytics_r_vs_e_{names[i]}_cf_{CF}'
        for suffix in ['png', 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)
        plt.close(fig)

        fig, ax = plt.subplots()
        ax = sns.heatmap(s_max_df, cmap='YlGnBu', vmin=0, vmax=1)
        ax.set_xlabel('Repair Rate', fontsize=14)
        ax.set_ylabel('Error-free Replication Rate', fontsize=14)
        # ax.set_title(f'{title_maps[names[i]]}', fontsize=14)
        cbar = ax.collections[0].colorbar
        cbar.set_label('Maximum Selection Coefficient', fontsize=14)
        filename = f'{outdir}/heatmap_analytics_r_vs_e_{names[i]}_max_cf_{CF}'
        for suffix in ['png', 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)
        plt.close(fig)

    # r vs u plots
    r = np.arange(0, 1.05, 0.05)
    u = np.arange(0, 1.05, 0.05)
    s_min = []
    s_max = []
    anything_fits = False
    for u_val in u:
        u_row_min = []
        u_row_max = []
        u_val = round(u_val, 2)
        for r_val in r:
            r_val = round(r_val, 2)
            subdf = df.loc[(df['r'].apply(rounding) == r_val) & (df['u'].apply(rounding) == u_val)]
            if not subdf.empty:
                anything_fits = True
            if len(np.unique(subdf['r'])) > 1 or len(np.unique(subdf['u'])) > 1:
                print('TROUBLE!!')
            u_row_min.append(round(subdf['min s'].min(), 2))
            u_row_max.append(round(subdf['max s'].max(), 2))
        s_min.append(u_row_min)
        s_max.append(u_row_max)
    # reverse so e=0 is the bottom row
    s_min.reverse()
    s_max.reverse()
    # e.reverse()
    s_min_df = pd.DataFrame(s_min, index=format(np.flip(u)), columns=format(r))
    s_max_df = pd.DataFrame(s_max, index=format(np.flip(u)), columns=format(r))
    if anything_fits:
        fig, ax = plt.subplots()
        ax = sns.heatmap(s_min_df, cmap='YlGnBu', vmin=0, vmax=1)
        ax.set_xlabel('Repair Rate', fontsize=14)
        ax.set_ylabel('Mutagenic Repair Rate', fontsize=14)
        # ax.set_title(f'{title_maps[names[i]]}', fontsize=14)
        cbar = ax.collections[0].colorbar
        cbar.set_label('Minimum Selection Coefficient', fontsize=14)
        filename = f'{outdir}/heatmap_analytics_r_vs_u_{names[i]}_cf_{CF}'
        for suffix in ['png', 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)
        plt.close(fig)

        fig, ax = plt.subplots()
        ax = sns.heatmap(s_max_df, cmap='YlGnBu', vmin=0, vmax=1)
        ax.set_xlabel('Repair Rate', fontsize=14)
        ax.set_ylabel('Mutagenic Repair Rate', fontsize=14)
        # ax.set_title(f'{title_maps[names[i]]}', fontsize=14)
        cbar = ax.collections[0].colorbar
        cbar.set_label('Maximum Selection Coefficient', fontsize=14)
        filename = f'{outdir}/heatmap_analytics_r_vs_u_{names[i]}_max_cf_{CF}'
        for suffix in ['png', 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)
        plt.close(fig)


