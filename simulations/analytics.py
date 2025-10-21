import numpy as np
import pandas as pd
import os

if not os.path.isdir('analytics_results/'):
    os.mkdir('analytics_results/')

cf_mult = 2.58

conf_interval_df = pd.read_csv('C3H_confidence_intervals_poisson.csv')
conf_interval_groups = conf_interval_df.groupby('Gene_name')

CF = 95
for gene, grp in conf_interval_groups:
    srted_grp = grp.sort_values(by='division')
    grp_total = grp['N_samples'].astype(float).sum()
    if gene == 'Egfr':
        egfr_low = srted_grp[f'Lower_CF_{CF}'].astype(float).values/grp_total
        egfr_high = srted_grp[f'Upper_CF_{CF}'].astype(float).values/grp_total
    elif gene == 'Braf':
        braf_low = srted_grp[f'Lower_CF_{CF}'].astype(float).values/grp_total
        braf_high = srted_grp[f'Upper_CF_{CF}'].astype(float).values/grp_total
    elif gene == 'Hras1':
        hras1_low = srted_grp[f'Lower_CF_{CF}'].astype(float).values/grp_total
        hras1_high = srted_grp[f'Upper_CF_{CF}'].astype(float).values/grp_total
    elif gene == 'Hras0':
        hras0_low = srted_grp[f'Lower_CF_{CF}'].astype(float).values/grp_total
        hras0_high = srted_grp[f'Upper_CF_{CF}'].astype(float).values/grp_total

egfr_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(egfr_low, egfr_high)]
braf_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(braf_low, braf_high)]
hras1_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(hras1_low, hras1_high)]
hras0_range_sd = [(round(lo, 2), round(hi, 2)) for lo, hi in zip(hras0_low, hras0_high)]


def mrca_mse(data, obs_range):
    to_return = 0
    for i in range(len(data)):
        error = max(obs_range[i][0]-data[i], data[i]-obs_range[i][1])
        if error < 0: # in range
            error = 0
        to_return += error**2
        # print(data[i], obs_range[i], error)
    return to_return**0.5

def main():
    df_columns = ['fitness', 'error free replication', 'repair', 'duplex rate',
                  'MRCA 0', 'MRCA 1', 'MRCA 2', 'MRCA 3', 'MRCA 4', 'MRCA 5-10', 'MRCA 5-300',
                  'egfr_c3h_sd_mse', 'braf_c3h_sd_mse', 'hras1_c3h_sd_mse', 'hras0_c3h_sd_mse',
                  'b', 'p', 'q', 'q/p', 'w_prime', 'w_prime - e*p_prime - (1-e)*q_prime']
    rows = []

    rows_match_hras0_sd = []
    rows_match_hras1_sd = []
    rows_match_egfr_sd = []
    rows_match_braf_sd = []

    c = 0

    min_s_per_param_hras0 = [] # fill with (r, u, e, min s, max s)
    min_s_per_param_hras1 = []
    min_s_per_param_braf = []
    min_s_per_param_egfr = []

    for r in np.arange(0, 1.05, 0.05):
        rows = []
        # for u in np.arange(0, 1.05, 0.05):
        for u in np.arange(0, r+0.05, 0.05):
            for e in np.arange(0, 1.05, 0.05):
                min_s_per_param_egfr.append([r, u, e, None, None])
                min_s_per_param_braf.append([r, u, e, None, None])
                min_s_per_param_hras1.append([r, u, e, None, None])
                min_s_per_param_hras0.append([r, u, e, None, None])
                for s in np.arange(0.05, 1.05, 0.05):
                    parameter_row = [s, e, r, u]

                    b = (1 - c + s) / 2

                    w_num = (1+s)*((1-e)*(1-r) + e*u) + 2*r*(1-e)
                    w_denom = 2 - (1-e)*(1-r)*(1-s) - e*(1-r) - e*u*(1-s)
                    w_prime = w_num/w_denom

                    P_prime = 0.5*w_prime*(1-r+u*(1-s)) + 0.5*u*(1+s)
                    P = P_prime*(s/b)
                    Q_prime = 0.5*w_prime*(1-r)*(1-s) + 0.5*(1-r)*(1+s) + r
                    Q = Q_prime*(s/b)

                    # assert w_prime == e*P_prime + (1-e)*Q_prime
                    # assert P <= Q

                    prob_row = [b, P, Q, Q/P, w_prime, w_prime-(e*P_prime+(1-e)*Q_prime)]

                    prob_lw_lw = 0.5*e*(1 - r + u*(1-s))*P
                    prob_lm_lw = 0.5*(1-e)*(1 - r + u*(1-s))*Q
                    prob_mm_lw = u*((1+s)/2 - s*w_prime)*(s/b)

                    prob_lw_lm = 0.5*e*(1-r)*(1-s)*P
                    prob_lm_lm = 0.5*(1-e)*(1-r)*(1-s)*Q
                    prob_mm_lm = (1-r)*((1+s)/2 - s*w_prime)*(s/b) + r*(1-s)*(s/b)

                    prob_mm_mm = (1-s)*(s/b)

                    prob_mrca_lw = u*s*(e + (1-e)*Q_prime/P_prime)
                    prob_mrca_lm = r*s/Q_prime + (1-r)*s*(1 - e + e*P_prime/Q_prime)

                    prob_mrca_mm = s

                    # mrca recursion
                    prob_mrca_mat = [[prob_mrca_lw, prob_mrca_lm, prob_mrca_mm]] # index 0: [prob(mrca 0|lw), p(mrca 0|lm), prob(mrca 0|mm)]
                    for i in range(1, 301):
                        p_mrca_i_lw = (1/P)*(prob_lw_lw*prob_mrca_mat[i-1][0]
                                             + prob_lm_lw*prob_mrca_mat[i-1][1]
                                             + prob_mm_lw*prob_mrca_mat[i-1][2])
                        p_mrca_i_lm = (1/Q)*(prob_lw_lm*prob_mrca_mat[i-1][0]
                                             + prob_lm_lm*prob_mrca_mat[i-1][1]
                                             + prob_mm_lm*prob_mrca_mat[i-1][2])
                        p_mrca_i_mm = (b/s)*(prob_mm_mm*prob_mrca_mat[i-1][2])
                        prob_mrca_mat.append([p_mrca_i_lw, p_mrca_i_lm, p_mrca_i_mm])

                    prob_mrca_mat = np.array(prob_mrca_mat)
                    mrca_0_4 = prob_mrca_mat[:5, 0]
                    prob_mrca_5_10 = np.sum(prob_mrca_mat[5:11, 0])
                    prob_mrca_5_300 = np.sum(prob_mrca_mat[5:300, 0])
                    # print(mrca_0_4)
                    mrca_row = np.concatenate([mrca_0_4 ,[prob_mrca_5_10, prob_mrca_5_300]])
                    # print(mrca_row)


                    mrca_row_for_mse = np.concatenate([mrca_0_4, [prob_mrca_5_300]])

                    egfr_mse_sd = mrca_mse(mrca_row_for_mse, egfr_range_sd)
                    braf_mse_sd = mrca_mse(mrca_row_for_mse, braf_range_sd)
                    hras_1_mse_sd = mrca_mse(mrca_row_for_mse, hras1_range_sd)
                    hras_0_mse_sd = mrca_mse(mrca_row_for_mse, hras0_range_sd)
                    mse_row = [egfr_mse_sd, braf_mse_sd, hras_1_mse_sd, hras_0_mse_sd]

                    full_row = np.concatenate([parameter_row, mrca_row, mse_row, prob_row])
                    rows.append(full_row)

                    # Based on std dev
                    if mrca_row[0] >= egfr_range_sd[0][0] and mrca_row[0] <= egfr_range_sd[0][1]:
                        if mrca_row[1] >= egfr_range_sd[1][0] and mrca_row[1] <= egfr_range_sd[1][1]:
                            if mrca_row[2] >= egfr_range_sd[2][0] and mrca_row[2] <= egfr_range_sd[2][1]:
                                if mrca_row[3] >= egfr_range_sd[3][0] and mrca_row[3] <= egfr_range_sd[3][1]:
                                    if mrca_row[4] >= egfr_range_sd[4][0] and mrca_row[4] <= egfr_range_sd[4][1]:
                                        if mrca_row[-1] >= egfr_range_sd[5][0] and mrca_row[-1] <= egfr_range_sd[5][1]:
                                            rows_match_egfr_sd.append(full_row)
                                            if min_s_per_param_egfr[-1][-2] is None:
                                                min_s_per_param_egfr[-1][-2] = s
                                            min_s_per_param_egfr[-1][-1] = s
                    if mrca_row[0] >= braf_range_sd[0][0] and mrca_row[0] <= braf_range_sd[0][1]:
                        if mrca_row[1] >= braf_range_sd[1][0] and mrca_row[1] <= braf_range_sd[1][1]:
                            if mrca_row[2] >= braf_range_sd[2][0] and mrca_row[2] <= braf_range_sd[2][1]:
                                if mrca_row[3] >= braf_range_sd[3][0] and mrca_row[3] <= braf_range_sd[3][1]:
                                    if mrca_row[4] >= braf_range_sd[4][0] and mrca_row[4] <= braf_range_sd[4][1]:
                                        if mrca_row[-1] >= braf_range_sd[5][0] and mrca_row[-1] <= braf_range_sd[5][1]:
                                            rows_match_braf_sd.append(full_row)
                                            if min_s_per_param_braf[-1][-2] is None:
                                                min_s_per_param_braf[-1][-2] = s
                                            min_s_per_param_braf[-1][-1] = s
                    if mrca_row[0] >= hras1_range_sd[0][0] and mrca_row[0] <= hras1_range_sd[0][1]:
                        if mrca_row[1] >= hras1_range_sd[1][0] and mrca_row[1] <= hras1_range_sd[1][1]:
                            if mrca_row[2] >= hras1_range_sd[2][0] and mrca_row[2] <= hras1_range_sd[2][1]:
                                if mrca_row[3] >= hras1_range_sd[3][0] and mrca_row[3] <= hras1_range_sd[3][1]:
                                    if mrca_row[4] >= hras1_range_sd[4][0] and mrca_row[4] <= hras1_range_sd[4][1]:
                                        if mrca_row[-1] >= hras1_range_sd[5][0] and mrca_row[-1] <= hras1_range_sd[5][1]:
                                            rows_match_hras1_sd.append(full_row)
                                            if min_s_per_param_hras1[-1][-2] is None:
                                                min_s_per_param_hras1[-1][-2] = s
                                            min_s_per_param_hras1[-1][-1] = s
                    if mrca_row[0] >= hras0_range_sd[0][0] and mrca_row[0] <= hras0_range_sd[0][1]:
                        if mrca_row[1] >= hras0_range_sd[1][0] and mrca_row[1] <= hras0_range_sd[1][1]:
                            if mrca_row[2] >= hras0_range_sd[2][0] and mrca_row[2] <= hras0_range_sd[2][1]:
                                if mrca_row[3] >= hras0_range_sd[3][0] and mrca_row[3] <= hras0_range_sd[3][1]:
                                    if mrca_row[4] >= hras0_range_sd[4][0] and mrca_row[4] <= hras0_range_sd[4][1]:
                                        if mrca_row[-1] >= hras0_range_sd[-1][0] and mrca_row[-1] <= hras0_range_sd[-1][1]:
                                            rows_match_hras0_sd.append(full_row)
                                            if min_s_per_param_hras0[-1][-2] is None:
                                                min_s_per_param_hras0[-1][-2] = s
                                            min_s_per_param_hras0[-1][-1] = s

        df = pd.DataFrame(rows, columns=df_columns)
        df.to_csv(f'analytics_results/r_{round(r, 2)}.csv', index=False)
        rows = []


    df_columns_params = ['r', 'u', 'e', 'min s', 'max s']

    df = pd.DataFrame(rows_match_braf_sd, columns=df_columns)
    df.to_csv(f'analytics_results/matching_braf_sd_cf_{CF}.csv', index=False)
    df = pd.DataFrame(min_s_per_param_braf, columns=df_columns_params)
    df.to_csv(f'analytics_results/min_max_s_per_param_braf_cf_{CF}.csv', index=False)

    df = pd.DataFrame(rows_match_hras1_sd, columns=df_columns)
    df.to_csv(f'analytics_results/matching_hras1_sd_cf_{CF}.csv', index=False)
    df = pd.DataFrame(min_s_per_param_hras1, columns=df_columns_params)
    df.to_csv(f'analytics_results/min_max_s_per_param_hras1_cf_{CF}.csv', index=False)

    df = pd.DataFrame(rows_match_hras0_sd, columns=df_columns)
    df.to_csv(f'analytics_results/matching_hras0_sd_cf_{CF}.csv', index=False)
    df = pd.DataFrame(min_s_per_param_hras0, columns=df_columns_params)
    df.to_csv(f'analytics_results/min_max_s_per_param_hras0_cf_{CF}.csv', index=False)

    df = pd.DataFrame(rows_match_egfr_sd, columns=df_columns)
    df.to_csv(f'analytics_results/matching_egfr_sd_cf_{CF}.csv', index=False)
    df = pd.DataFrame(min_s_per_param_egfr, columns=df_columns_params)
    df.to_csv(f'analytics_results/min_max_s_per_param_egfr_cf_{CF}.csv', index=False)







if __name__ == '__main__':
    main()
