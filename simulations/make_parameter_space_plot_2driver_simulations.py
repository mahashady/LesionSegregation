
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

indir = os.path.expanduser('~/Desktop/drivers/simulation_grid_2driver/')
outdir = os.path.expanduser('~/Desktop/drivers/figures/parameter_space_plots/')


CF = 95

files = [f'params_fit_braf_cf_{CF}.csv', f'params_fit_egfr_cf_{CF}.csv', f'params_fit_hras0_cf_{CF}.csv', f'params_fit_hras1_cf_{CF}.csv']
names = ['braf', 'egfr', 'hras0', 'hras1']
title_maps = {'braf': 'BRAF V637E', 'egfr': 'EGFR F254I', 'hras0': 'HRAS Q61L', 'hras1': 'HRAS Q61R'}


format = np.vectorize(lambda x: round(x, 2))

e_range = np.arange(0.05, 0.55, 0.1)
r_range = np.arange(0.2, 0.9, 0.1)
u_range = np.arange(0.05, 0.8, 0.1)



for i in range(len(files)):
    df = pd.read_csv(f'{indir}/{files[i]}') #.dropna(subset='min s')
    print('here')
    s_min = []
    s_max = []
    anything_fits = False
    for e_val in e_range:
        e_val = round(e_val, 2)
        e_row_min = []
        e_row_max = []
        for r_val in r_range:
            r_val = round(r_val, 2)
            subdf = df.loc[(df['r'] == r_val) & (df['e'] == e_val)]
            # print(r_val, e_val, subdf)

            if not subdf.empty:
                anything_fits = True
            if len(np.unique(subdf['r'])) > 1 or len(np.unique(subdf['e'])) > 1:
                print('TROUBLE!!')
            e_row_min.append(round(subdf['s'].min(), 2))
            e_row_max.append(round(subdf['s'].max(), 2))
        s_min.append(e_row_min)
        s_max.append(e_row_max)

    # reverse so e=0 is the bottom row
    s_min.reverse()
    s_max.reverse()
    s_min_df = pd.DataFrame(s_min, index=format(np.flip(e_range)), columns=format(r_range))
    s_max_df = pd.DataFrame(s_max, index=format(np.flip(e_range)), columns=format(r_range))
    if anything_fits:
        print(anything_fits)
        fig, ax = plt.subplots()
        ax = sns.heatmap(s_min_df, cmap='YlGnBu', vmin=0, vmax=1)
        ax.set_xlabel('Repair Rate', fontsize=14)
        ax.set_ylabel('Error-free Replication Rate', fontsize=14)
        # ax.set_title(f'{title_maps[names[i]]}', fontsize=14)
        cbar = ax.collections[0].colorbar
        cbar.set_label('Minimum Selection Coefficient', fontsize=14)
        filename = f'{outdir}/heatmap_2driver_simulations_r_vs_e_{names[i]}_cf_{CF}'
        for suffix in ['png', 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)
        plt.close(fig)


        # fig, ax = plt.subplots()
        # ax = sns.heatmap(s_max_df, cmap='YlGnBu', vmin=0, vmax=1)
        # ax.set_xlabel('repair rate')
        # ax.set_ylabel('error-free replication rate')
        # fig.savefig(f'{outdir}/heatmap_2driver_simulations_r_vs_e_{names[i]}_max_cf_{CF}.pdf', bbox_inches='tight')
        # plt.close(fig)

for i in range(len(files)):
    df = pd.read_csv(f'{indir}/{files[i]}')  # .dropna(subset='min s')
    # r vs u plots
    s_min = []
    s_max = []
    anything_fits = False
    for u_val in u_range:
        u_row_min = []
        u_row_max = []
        u_val = round(u_val, 2)
        for r_val in r_range:
            r_val = round(r_val, 2)
            subdf = df.loc[(df['r'] == r_val) & (df['u'] == u_val)]
            if not subdf.empty:
                anything_fits = True
            if len(np.unique(subdf['r'])) > 1 or len(np.unique(subdf['u'])) > 1:
                print('TROUBLE!!')
            u_row_min.append(round(subdf['s'].min(), 2))
            u_row_max.append(round(subdf['s'].max(), 2))
        s_min.append(u_row_min)
        s_max.append(u_row_max)
    # reverse so e=0 is the bottom row
    s_min.reverse()
    s_max.reverse()
    # e.reverse()
    s_min_df = pd.DataFrame(s_min, index=format(np.flip(u_range)), columns=format(r_range))
    s_max_df = pd.DataFrame(s_max, index=format(np.flip(u_range)), columns=format(r_range))
    if anything_fits:
        print(anything_fits)
        fig, ax = plt.subplots()
        ax = sns.heatmap(s_min_df, cmap='YlGnBu', vmin=0, vmax=1)
        ax.set_xlabel('Repair Rate', fontsize=14)
        ax.set_ylabel('Mutagenic Repair Rate', fontsize=14)
        # ax.set_title(f'{title_maps[names[i]]}', fontsize=14)
        cbar = ax.collections[0].colorbar
        cbar.set_label('Minimum Selection Coefficient', fontsize=14)
        filename = f'{outdir}/heatmap_2driver_simulations_r_vs_u_{names[i]}_cf_{CF}'
        for suffix in ['png', 'svg', 'pdf']:
            outname = f'{filename}.{suffix}'
            plt.savefig(outname, format=suffix, bbox_inches='tight', dpi=300)
        plt.close(fig)

        # fig, ax = plt.subplots()
        # ax = sns.heatmap(s_max_df, cmap='YlGnBu', vmin=0, vmax=1)
        # ax.set_xlabel('repair rate')
        # ax.set_ylabel('mutagenic repair rate')
        # fig.savefig(f'{outdir}/heatmap_2driver_simulations_r_vs_u_{names[i]}_max_cf_{CF}.pdf', bbox_inches='tight')
        # plt.close(fig)


