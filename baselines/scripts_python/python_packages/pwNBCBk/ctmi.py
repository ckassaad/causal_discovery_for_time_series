#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Mutual information and conditional mutual information between time series: script implementing
the TMI and CTMI methods.

Date: Dec 2019
Author: Karim Assaad, karimassaad3@gmail.com, karim.assaad@univ.grenoble.alpes.fr, karim.assaad@coservit.com
"""

import numpy as np
import pandas as pd

from baselines.scripts_python.python_packages.pwNBCBk.tigramite.tigramite.independence_tests import CMIknn
from baselines.scripts_python.python_packages.pwNBCBk.tigramite.tigramite.independence_tests import ParCorr

# from tools.functions import ts_order

# from scipy.special import digamma, erfinv
# from sklearn.neighbors import KDTree

# try:
#     from tools import fast_functions
# except:
#     print("Could not import cython packages")


def indep_test(x, y, z=None, sig_samples=10000, p_value=True, measure="cmiknn", k=10):
    if measure == "cmiknn":
        cd = CMIknn(mask_type=None, significance='shuffle_test', fixed_thres=None, sig_samples=sig_samples,
                    sig_blocklength=3, knn=k, confidence='bootstrap', conf_lev=0.9, conf_samples=10000,
                    conf_blocklength=1, verbosity=0)
    elif measure == "parcorr":
        cd = ParCorr(mask_type=None, significance='shuffle_test', fixed_thres=None, sig_samples=sig_samples,
                     sig_blocklength=3, confidence='bootstrap', conf_lev=0.9, conf_samples=10000, conf_blocklength=1,
                     verbosity=0)
    else:
        cd = None
        print("Independence measure '" + str(measure) + "' do not exist.")
        exit(0)
    dim_x = x.shape[1]
    dim_y = y.shape[1]
    if z is not None:
        dim_z = z.shape[1]
        X = np.concatenate((x.values, y.values, z.values), axis=1)
        xyz = np.array([0] * dim_x + [1] * dim_y + [2] * dim_z)
    else:
        X = np.concatenate((x, y), axis=1)
        xyz = np.array([0] * dim_x + [1] * dim_y)

    value = cd.get_dependence_measure(X.T, xyz)
    if p_value:
        pvalue = cd.get_shuffle_significance(X.T, xyz, value)
    else:
        pvalue = value

    # N = x.shape[0]
    # dim_x = x.shape[1]
    # dim_y = y.shape[1]
    # v = np.concatenate((x, y), axis=1)
    # # if gamma==1:
    # #     print(v)
    # #     print("----------------")
    # kdt = KDTree(v, leaf_size=30, metric='infinity')
    # dist, ind = kdt.query(v, k=k+1, return_distance=True)
    # epsarray = dist[:, k]
    #
    # nx, ny = fast_functions._num_neighbors_xy_cython(v, N, dim_x, dim_y, epsarray, k)
    #
    # # nx = np.zeros(N, dtype='int32')
    # # ny = np.zeros(N, dtype='int32')
    # # for i in range(N):
    # #     for j in range(N):
    # #         dx = 0.
    # #         for d in range(0, dim_y):
    # #             dx = max(abs(v[i, d] - v[j, d]), dx)
    # #         dy = 0.
    # #         for d in range(dim_y, dim_x+dim_y):
    # #             dy = max(abs(v[i, d] - v[j, d]), dy)
    # #
    # #         # For no conditions, kz is counted up to T
    # #         if dx < epsarray[i]:
    # #             nx[i] += 1
    # #         if dy < epsarray[i]:
    # #             ny[i] += 1
    #
    # res = digamma(k) + digamma(N) - np.mean(digamma(nx) + digamma(ny))

    return pvalue, value
    # return value


def get_alpha(mts, k=10):
    mi_list = []
    for i in range(mts.shape[1]):
        for t in range(100):
            ts_i = mts[mts.columns[i]].dropna().to_frame()
            ts_j = 0.05*ts_i + 0.95*np.random.randn(ts_i.shape[0], ts_i.shape[1])
            pval, val = tmi(ts_i, ts_j, k=k, p_value=False)
            mi_list.append(val)
    alpha = abs(max(mi_list))
    return alpha


def window_size(ts, alpha, lag_max=5, k=10, p_value=True):
    # def get_opt_lambda(l_list):
    #     candidate = np.min(l_list)
    #     if candidate + 1 in l_list:
    #         l_list.remove(candidate)
    #         return get_opt_lambda(l_list)
    #     else:
    #         return candidate + 2
    #
    # mi_list = []
    # for i in range(2, lag_max+1):
    #     wts = window_representation(ts, windows_size=i)
    #     i_data = wts[wts.columns[0]]
    #     j_data = wts[wts.columns[i-1]]
    #
    #     mi_pval, mi_val = tmi(i_data, j_data, k=k, p_value=False)
    #     if mi_val == np.inf:
    #         mi_val = 1
    #     mi_list.append(mi_val)
    #
    # mi_array = np.array(mi_list)
    # test = mi_array <= alpha
    # if sum(test) == len(test):
    #     window = 1
    # else:
    #     upper = np.max(mi_list)
    #     lower = np.min(mi_list)
    #     lower = 0.9 * (upper - lower) + lower
    #
    #     mi_array = np.array(mi_list)
    #     l_set = list(np.where(mi_array >= lower)[0])
    #     window = get_opt_lambda(l_set)

    mi_list = []
    for i in range(2, lag_max+1):
        wts = window_representation(ts, windows_size=i)
        i_data = wts[wts.columns[0]]
        j_data = wts[wts.columns[i-1]]

        mi_pval, mi_val = tmi(i_data, j_data, k=k, p_value=False)
        if mi_val == np.inf:
            mi_val = 1
        mi_list.append(mi_val)

    mi_array = np.array(mi_list)

    j = mi_array.argmax()
    if p_value:
        wts = window_representation(ts, windows_size=j+2)
        i_data = wts[wts.columns[0]]
        j_data = wts[wts.columns[j+1]]

        mi_pval, mi_val = tmi(i_data, j_data, k=k, p_value=True)

        if mi_pval<=alpha:
            window = j+2
        else:
            window = 1
    else:
        if mi_array[j]>alpha:
            window = j + 2
        else:
            window = 1
    return window


def window_representation(ts, windows_size=4, overlap=True):
    ts = ts.dropna()
    if windows_size == 0:
        return ts.to_frame()
    else:
        ts_windows = pd.DataFrame()
        for i in range(windows_size):
            i_data = ts[i:(ts.shape[0]-windows_size+i+1)].values
            ts_windows.loc[:, str(ts.name)+"_"+str(i+1)] = i_data
        if not overlap:
            ts_windows = ts_windows.iloc[::windows_size, :]
    return ts_windows


def align_pair(x, y, sampling_rate_tuple=(1, 1), k=10, max_gamma=5, set_numbers="Z"):
    """
    :param x:
    :param y:
    :param sampling_rate_tuple:
    :param k:
    :param max_gamma:
    :param set_numbers: "Z", "N" or "-N"
    :return:
    """
    c1 = list()
    c2 = list()

    _, val = tmi(x, y, sampling_rate_tuple, k=k, p_value=False)
    c1.append(val)
    c2.append(val)

    if (set_numbers == "Z") or (set_numbers == "N"):
        for g in range(1, max_gamma):
            _, val = tmi(x, y, sampling_rate_tuple, k=k, gamma=g, p_value=False)
            c1.append(val)

    if (set_numbers == "Z") or (set_numbers == "-N"):
        for g in range(1, max_gamma):
            _, val = tmi(x, y, sampling_rate_tuple, k=k, gamma=-g, p_value=False)
            c2.append(val)

    if np.max(c1) >= np.max(c2):
        g = np.argmax(c1)
    else:
        g = -np.argmax(c2)
    print(c1, c2)
    return g


def align_one(x, y, sampling_rate_tuple=(1, 1), k=10, max_gamma=5):
    c1 = list()
    c2 = list()

    c = 0
    c1.append(c)
    c2.append(c)

    for g in range(1, max_gamma):
        _, val = tmi(x, y, sampling_rate_tuple, k=k, gamma=g, p_value=False)
        c1.append(val)

    g = np.argmax(c1)
    return g


# def update_gamma_for_d_ctmi(x, y, z, sampling_rate_tuple=(1, 1, 1), k=10, max_gamma=5):
#     c1 = list()
#     c2 = list()
#
#     _, val = tmi(x, y, z, sampling_rate_tuple, k=k, p_value=False)
#     c1.append(val)
#     c2.append(val)
#
#     for gxz in range(0, max_gamma):
#         for gyz in range(0, max_gamma):
#             _, val = tmi(x, y, sampling_rate_tuple, k=k, gamma=g, p_value=False)
#             c1.append(val)
#         for gyz in range(0, max_gamma):
#             _, val = tmi(x, y, sampling_rate_tuple, k=k, gamma=-g, p_value=False)
#             c1.append(val)
#     for gxz in range(0, max_gamma):
#         for gyz in range(0, max_gamma):
#             _, val = tmi(x, y, sampling_rate_tuple, k=k, gamma=g, p_value=False)
#             c1.append(val)
#         for gyz in range(0, max_gamma):
#             _, val = tmi(x, y, sampling_rate_tuple, k=k, gamma=-g, p_value=False)
#             c1.append(val)
#     g = np.argmax(c1)
#     return g

def align_matrix(data_dict, keys, sampling_rates, k=10, max_gamma=5):
    d = len(keys)
    g_matrix = np.zeros([d, d], dtype=int)
    for i in range(d):
        for j in range(i, d):
            if i != j:
                x = data_dict[keys[i]]
                y = data_dict[keys[j]]
                g = align_pair(x, y, (sampling_rates[keys[i]], sampling_rates[keys[j]]), k=k, max_gamma=max_gamma)
                g_matrix[i, j] = g
                g_matrix[j, i] = -g
            else:
                x = data_dict[keys[i]]
                y = data_dict[keys[j]]
                g = align_one(x, y, (sampling_rates[keys[i]], sampling_rates[keys[j]]), k=k, max_gamma=max_gamma)
                g_matrix[i, j] = g

    return pd.DataFrame(g_matrix, columns=keys, index=keys)


def tmi(x, y, sampling_rate_tuple=(1, 1), k=10, gamma=0, p_value=True, sig_samples=10000):
    sr1, sr2 = sampling_rate_tuple
    dsr = abs(sr1 - sr2)
    iter1 = (sr1 % sr2)*dsr
    iter2 = (sr2 % sr1)*dsr
    if gamma > 0:
        y = y[gamma:]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        x = x[x.index % (iter1+1) == 0]
        y = y[y.index % (iter2+1) == 0]

        x = x[:-gamma]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)

    elif gamma < 0:
        x = x[-gamma:]
        y = y.reset_index(drop=True)
        x = x.reset_index(drop=True)
        y = y[y.index % (iter2+1) == 0]
        x = x[x.index % (iter1+1) == 0]

        y = y[:gamma]
        y = y.reset_index(drop=True)
        x = x.reset_index(drop=True)

#        x = x[-gamma:]
#        y = y[:gamma]
    else:
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        x = x[x.index % (iter1 + 1) == 0]
        y = y[y.index % (iter2 + 1) == 0]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)

    m = min(x.shape[0], y.shape[0])
    x = x[:m]
    y = y[:m]

    if len(x.shape) == 1:
        x = x.to_frame()
    if len(y.shape) == 1:
        y = y.to_frame()

    mi_pval, mi_val = indep_test(x, y, sig_samples=sig_samples, p_value=p_value, measure="cmiknn", k=k)
    return mi_pval, mi_val


def get_index_of_aligned_dict(dict_data, keys):
    # check if dict is align
    legnth_list = []
    for k in keys:
        legnth_list.append(dict_data[k].shape[0])
    if legnth_list.count(legnth_list[0]) != len(legnth_list):
        print("Error: time series in dict are not aligned")
        exit(0)
    index_df = pd.DataFrame()
    for k in keys:
        index_df[k] = dict_data[k].index.map(str)
    index_df.insert(0, "concatenated", "id_")
    index_df["concatenated"] = index_df.sum(axis=1)
    index_df = index_df.set_index('concatenated')
    return index_df.index


def aligned_dict_to_df(dict_data):
    concat_df = pd.DataFrame()
    for k in dict_data.keys():
        if isinstance(dict_data[k], pd.Series):
            dict_data[k] = dict_data[k].to_frame()
        concat_df[dict_data[k].columns] = dict_data[k].reset_index(drop=True)
    return concat_df


def align_xy(x, y, gamma, sampling_rate_tuple):
    sr1, sr2 = sampling_rate_tuple
    dsr = abs(sr1 - sr2)
    iter1 = (sr1 % sr2) * dsr
    iter2 = (sr2 % sr1) * dsr
    idx_x = x.index
    idx_y = y.index
    if gamma > 0:
        y = y.loc[gamma:]
        idx_y = idx_y[gamma:]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        x = x[x.index % (iter1 + 1) == 0]
        y = y[y.index % (iter2 + 1) == 0]
        # idx_x = idx_x[x.index % (iter1 + 1) == 0]
        # idx_y = idx_y[y.index % (iter2 + 1) == 0]
        idx_x = idx_x[x.index]
        idx_y = idx_y[y.index]

        x = x[:-gamma]
        idx_x = idx_x[:-gamma]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)

    elif gamma < 0:
        x = x[-gamma:]
        idx_x = idx_x[-gamma:]
        y = y.reset_index(drop=True)
        x = x.reset_index(drop=True)
        y = y[y.index % (iter2 + 1) == 0]
        x = x[x.index % (iter1 + 1) == 0]
        # idx_y = idx_y[y.index % (iter2 + 1) == 0]
        # idx_x = idx_x[x.index % (iter1 + 1) == 0]
        idx_x = idx_x[x.index]
        idx_y = idx_y[y.index]

        y = y[:gamma]
        idx_y = idx_y[:gamma]
        y = y.reset_index(drop=True)
        x = x.reset_index(drop=True)
    else:
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        x = x[x.index % (iter1 + 1) == 0]
        y = y[y.index % (iter2 + 1) == 0]
        # idx_x = idx_x[x.index % (iter1 + 1) == 0]
        # idx_y = idx_y[y.index % (iter2 + 1) == 0]
        idx_x = idx_x[x.index]
        idx_y = idx_y[y.index]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)

    m = min(x.shape[0], y.shape[0])
    x = x[:m]
    y = y[:m]
    idx_x = idx_x[:m]
    idx_y = idx_y[:m]

    if len(x.shape) == 1:
        x = x.to_frame()
    if len(y.shape) == 1:
        y = y.to_frame()
    return x, y, idx_x, idx_y


def tmi_xy_z(x, y, z, sampling_rate_tuple, k, gamma, sig_samples=10000, measure="cmiknn"):
    sr1, sr2 = sampling_rate_tuple
    dsr = abs(sr1 - sr2)
    iter1 = (sr1 % sr2) * dsr
    iter2 = (sr2 % sr1) * dsr
    if gamma > 0:
        z = z[gamma:]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        z = z.reset_index(drop=True)
        x = x[x.index % (iter1 + 1) == 0]
        y = y[y.index % (iter1 + 1) == 0]
        z = z[z.index % (iter2 + 1) == 0]

        x = x[:-gamma]
        y = y[:-gamma]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        z = z.reset_index(drop=True)

    elif gamma < 0:
        x = x[-gamma:]
        y = y[-gamma:]
        z = z.reset_index(drop=True)
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        z = z[z.index % (iter2 + 1) == 0]
        x = x[x.index % (iter1 + 1) == 0]
        y = y[y.index % (iter1 + 1) == 0]

        z = z[:gamma]
        z = z.reset_index(drop=True)
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)

    else:
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        z = z.reset_index(drop=True)
        x = x[x.index % (iter1 + 1) == 0]
        y = y[y.index % (iter1 + 1) == 0]
        z = z[z.index % (iter2 + 1) == 0]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        z = z.reset_index(drop=True)

    m = min(x.shape[0], y.shape[0], z.shape[0])
    x = x[:m]
    y = y[:m]
    z = z[:m]

    if len(x.shape) == 1:
        x = x.to_frame()
    if len(y.shape) == 1:
        y = y.to_frame()
    if len(z.shape) == 1:
        z = z.to_frame()

    cmi_pval, cmi_val = indep_test(x, y, z=z, sig_samples=sig_samples, p_value=False, measure=measure, k=k)
    return cmi_pval, cmi_val


def align_cond(x, y, z, sampling_rate_tuple, k, max_gamma=5, measure="cmiknn", mission="dctmi"):
    z = z.loc[x.index[0]:]

    c1 = list()
    c2 = list()
    _, val = tmi_xy_z(x, y, z, sampling_rate_tuple, k=k, gamma=0, measure=measure)
    c1.append(val)
    c2.append(val)

    for g in range(1, max_gamma):
        _, val = tmi_xy_z(x, y, z, sampling_rate_tuple, k=k, gamma=g, measure=measure)
        c1.append(val)

    for g in range(1, max_gamma):
        _, val = tmi_xy_z(x, y, z, sampling_rate_tuple, k=k, gamma=-g, measure=measure)
        c2.append(val)
    if mission == "dctmi":
        if np.min(c1) <= np.min(c2):
            g = np.argmin(c1)
        else:
            g = -np.argmin(c2)
    else:
        if np.max(c1) >= np.max(c2):
            g = np.argmax(c1)
        else:
            g = -np.argmax(c2)
    return g


def align_xyz(xy, z, sampling_rate_dict, gamma_matrix, mission):
    names_xy = [*xy.keys()]
    name_x, name_y = names_xy[0], names_xy[1]
    sampling_rate_tuple = (sampling_rate_dict[name_x], sampling_rate_dict[name_y])

    g_xy = gamma_matrix[name_y].loc[name_x]
    v_x, v_y, idx_x, idx_y = align_xy(xy[name_x], xy[name_y], g_xy, sampling_rate_tuple)

    v_x2 = v_x.copy()
    v_y2 = v_y.copy()
    idx_x2 = idx_x.copy()
    idx_y2 = idx_y.copy()

    names_z = [*z.keys()]

    v_z = dict()
    v_z2 = dict()
    nz_visted = []

    k = 0
    for nz in names_z:
        g_xz = gamma_matrix[nz].loc[name_x]
        g_yz = gamma_matrix[nz].loc[name_y]

        # todo add constaint on gamma depending if it's ictmi  or dctmi
        # if mission == "dctmi":
        #     if (g_xz > abs(g_xy)) and (g_yz > abs(g_xy)):
        #         g_xz, g_yz = align_triple(xy[name_x], xy[name_y], z, (sampling_rate_dict[name_x], sampling_rate_dict[name_y],
        #                                               sampling_rate_dict[nz]), k=k, max_gamma=5)
        # elif mission == "ictmi":
        #     if g_xz < 0:
        #         g_xz = align_one(xy[name_x], z, (sampling_rate_dict[name_x], sampling_rate_dict[nz]), k=k, max_gamma=5)
        #     if g_yz < 0:
        #         g_yz = align_one(xy[name_y], z, (sampling_rate_dict[name_x], sampling_rate_dict[nz]), k=k, max_gamma=5)

        if (g_xz is not pd.NA) and (g_yz is not pd.NA):
            if idx_x[0] <= idx_y[0]:
                idx_xy = idx_x
                idx_xy2 = idx_x2
            else:
                idx_xy = idx_y
                idx_xy2 = idx_y2

            # g = min(g_xz, g_yz)
            # diff = max(g_xz, g_yz) - g
            # windows_size = diff + 1 # abs(g_xy)
            # nz_processed = concat_windows(z[nz], windows_size=windows_size, name=nz)

            if g_xz == g_yz:
                g = g_xz
            else:
                # v_xy = pd.concat([v_x.reindex(idx_xy), v_y.reindex(idx_xy)], axis=1)
                # print("v_xy: "+str(v_xy.shape))
                # g = align_pair(v_xy, z[nz], sampling_rate_tuple=(1, 1), k=10, max_gamma=5)
                sampling_rate_tuple_xyz = (max(sampling_rate_tuple), sampling_rate_dict[nz])
                v_x_new = v_x.copy()
                v_y_new = v_y.copy()
                v_x_new.index = idx_xy
                v_y_new.index = idx_xy
                g = align_cond(v_x_new, v_y_new, z[nz], sampling_rate_tuple_xyz,
                               k, max_gamma=5, measure="cmiknn", mission=mission)
            nz_processed = z[nz]

        elif g_xz is not pd.NA:
            idx_xy = idx_x
            idx_xy2 = idx_x2
            g = g_xz
            nz_processed = z[nz]
        elif g_yz is not pd.NA:
            idx_xy = idx_y
            idx_xy2 = idx_y2
            g = g_yz
            nz_processed = z[nz]
        else:
            g = np.nan
            nz_processed = None
            idx_xy = None
            idx_xy2 = None

        if g is not pd.NA:
            xyz_dict = {name_x: v_x, nz: nz_processed}
            xyz_dict[name_x].index = idx_xy

            xyz_dict2 = {name_y: v_y2, nz: nz_processed}
            xyz_dict2[name_y].index = idx_xy2

            sampling_rate_tuple = (sampling_rate_dict[name_x], sampling_rate_dict[nz])

            bool_idx = pd.DataFrame([False] * len(idx_xy), columns=['bool'])
            bool_idx.index = idx_xy
            bool_idx2 = pd.DataFrame([False] * len(idx_xy2), columns=['bool'])
            bool_idx2.index = idx_xy2

            v_x, z_processed, idx_x, _ = align_xy(xyz_dict[name_x], xyz_dict[nz], g, sampling_rate_tuple)
            bool_idx.loc[idx_x] = True
            bool_idx = bool_idx['bool'].values
            v_y = v_y[bool_idx]
            idx_y = idx_y[bool_idx]

            v_y2, z_processed2, idx_y2, _ = align_xy(xyz_dict2[name_y], xyz_dict[nz], g, sampling_rate_tuple)
            bool_idx2.loc[idx_y2] = True
            bool_idx2 = bool_idx2['bool'].values
            v_x2 = v_x2[bool_idx2]
            idx_x2 = idx_x2[bool_idx2]

            for nz_v in nz_visted:
                v_z[nz_v] = v_z[nz_v][bool_idx]
                v_z2[nz_v] = v_z2[nz_v][bool_idx2]
            v_z[nz] = z_processed
            v_z2[nz] = z_processed2
            nz_visted.append(nz)

    names_z = [*v_z.keys()]
    names_z2 = [*v_z.keys()]

    index_zx = get_index_of_aligned_dict(v_z, names_z)
    index_zy = get_index_of_aligned_dict(v_z2, names_z2)

    idx = index_zx.intersection(index_zy)
    if not (len(idx) == len(index_zx) == len(index_zy)):
        print("Something is wrong")
        input("enter")
    ##
    v_x = v_x.reset_index(drop=True)
    v_y = v_y.reset_index(drop=True)
    for nz_v in nz_visted:
        v_z[nz_v] = v_z[nz_v].reset_index(drop=True)

    return v_x, v_y, v_z


def ctmi(x, y, z, name_x, name_y, sampling_rate_dict, gamma_matrix, p_value=False, k=10, sig_samples=10000,
         mission="dctmi"):
    xy = {name_x: x, name_y: y}
    v_x, v_y, v_z = align_xyz(xy, z, sampling_rate_dict, gamma_matrix, mission)

    if len(v_x.shape) == 1:
        v_x = v_x.to_frame()
    if len(v_y.shape) == 1:
        v_y = v_y.to_frame()

    names_z = [*z.keys()]
    if len(names_z) > 0:
        v_z = aligned_dict_to_df(v_z)
    else:
        v_z = None
    cmi_pval, cmi_val = indep_test(v_x, v_y, z=v_z, sig_samples=sig_samples, p_value=p_value, measure="cmiknn", k=k)
    return cmi_pval, cmi_val


def get_sampling_rate(ts):
    # index of all non nan values in time series
    idx = np.argwhere(~np.isnan(ts).values)
    if len(idx) == len(ts):
        return True, 1
    # differentiate index, if there's no nan all values should be equal to 1
    diff = np.diff(idx, axis=0)
    udiff = np.unique(diff)
    if (len(udiff) == 1) and (udiff != 1):
        cd_bool = True
        cd_value = int(udiff)
    elif len(udiff) == 2:
        idx = np.argwhere(diff.reshape(-1) > 1)
        diff = diff[idx]
        udiff = np.unique(diff)
        if len(udiff) == 1:
            cd_bool = True
            cd_value = int(udiff)
        else:
            # ???
            cd_bool = False
            cd_value = np.nan
    else:
        # if there is no delay in the time series
        cd_bool = False
        cd_value = np.nan
    return cd_bool, cd_value


def i_ctmi(x, y, z, name_x, name_y, name_z, sampling_rate_dict, gamma_matrix, p_value=False, k=10, sig_samples=10000):
    # if u:
    #     xy = {name_x: x, name_y: y}
    #     sampling_rate_tuple = (sampling_rate_dict[name_x], sampling_rate_dict[name_y])
    #
    #     g_xy = gamma_matrix[name_y].loc[name_x]
    #     v_x, v_y, idx_x, idx_y = align_xy(xy[name_x], xy[name_y], g_xy, sampling_rate_tuple)
    #
    #     names_u = [*u.keys()]
    #
    #     v_u = dict()
    #     nu_visted = []
    #
    #     k = 0
    #     for nu in names_u:
    #         g_xu = gamma_matrix[nu].loc[name_x]
    #         g_yu = gamma_matrix[nu].loc[name_y]
    #
    #         if (g_xu is not pd.NA) and (g_yu is not pd.NA):
    #             if idx_x[0] <= idx_y[0]:
    #                 idx_xy = idx_x
    #             else:
    #                 idx_xy = idx_y
    #
    #             if g_xu == g_yu:
    #                 g = g_xu
    #             else:
    #                 sampling_rate_tuple_xyz = (max(sampling_rate_tuple), sampling_rate_dict[nu])
    #                 v_x_new = v_x.copy()
    #                 v_y_new = v_y.copy()
    #                 v_x_new.index = idx_xy
    #                 v_y_new.index = idx_xy
    #                 g = align_cond(v_x_new, v_y_new, u[nu], sampling_rate_tuple_xyz,
    #                                k, max_gamma=5, measure="cmiknn", mission="dctmi")
    #             nu_processed = u[nu]
    #
    #         elif g_xu is not pd.NA:
    #             idx_xy = idx_x
    #             g = g_xu
    #             nu_processed = u[nu]
    #         elif g_yu is not pd.NA:
    #             idx_xy = idx_y
    #             g = g_yu
    #             nu_processed = u[nu]
    #         else:
    #             g = np.nan
    #             nu_processed = None
    #             idx_xy = None
    #
    #         if g is not pd.NA:
    #             xyu_dict = {name_x: v_x, nu: nu_processed}
    #             xyu_dict[name_x].index = idx_xy
    #
    #             sampling_rate_tuple = (sampling_rate_dict[name_x], sampling_rate_dict[nu])
    #
    #             bool_idx = pd.DataFrame([False] * len(idx_xy), columns=['bool'])
    #             bool_idx.index = idx_xy
    #
    #             v_x, u_processed, idx_x, _ = align_xy(xyu_dict[name_x], xyu_dict[nu], g, sampling_rate_tuple)
    #             bool_idx.loc[idx_x] = True
    #             bool_idx = bool_idx['bool'].values
    #             v_y = v_y[bool_idx]
    #             idx_y = idx_y[bool_idx]
    #
    #             for nu_v in nu_visted:
    #                 v_u[nu_v] = v_u[nu_v][bool_idx]
    #             v_u[nu] = u_processed
    #             nu_visted.append(nu)
    #
    # else:
    v_x = x.copy()
    v_y = y.copy()
    v_z = z.copy()

    v_z_x = dict()
    v_z_y = dict()

    # name_z = z.names

    g_xz = gamma_matrix[name_z].loc[name_x]
    g_yz = gamma_matrix[name_z].loc[name_y]

    xz_dict = {name_x: v_x, name_z: v_z}
    yz_dict = {name_y: v_y, name_z: v_z}
    # xz_dict = {name_x: v_x, name_z: z[name_z]}
    # yz_dict = {name_y: v_y, name_z: z[name_z]}

    sampling_rate_tuple_x = (sampling_rate_dict[name_x], sampling_rate_dict[name_z])
    sampling_rate_tuple_y = (sampling_rate_dict[name_y], sampling_rate_dict[name_z])

    v_x, z_processed_x, idx_x, idx_z_x = align_xy(xz_dict[name_x], xz_dict[name_z], g_xz, sampling_rate_tuple_x)
    v_y, z_processed_y, idx_y, idx_z_y = align_xy(yz_dict[name_y], yz_dict[name_z], g_yz, sampling_rate_tuple_y)

    v_z_x[name_z] = z_processed_x
    v_z_y[name_z] = z_processed_y

    idx = idx_z_x.intersection(idx_z_y)

    bool_idx_x = pd.DataFrame([False] * len(v_x.index), columns=['bool'])
    bool_idx_x.index = idx_z_x
    bool_idx_x.loc[idx] = True
    bool_idx_x = bool_idx_x['bool'].values

    bool_idx_y = pd.DataFrame([False] * len(v_y.index), columns=['bool'])
    bool_idx_y.index = idx_z_y
    bool_idx_y.loc[idx] = True
    bool_idx_y = bool_idx_y['bool'].values

    v_x = v_x.loc[bool_idx_x]
    v_y = v_y.loc[bool_idx_y]
    v_x = v_x.reset_index(drop=True)
    v_y = v_y.reset_index(drop=True)
    v_z = v_z_x.copy()

    bool_idx_z = pd.DataFrame([False] * len(z_processed_x.index), columns=['bool'])
    bool_idx_z.index = idx_z_x
    bool_idx_z.loc[idx] = True
    bool_idx_z = bool_idx_z['bool'].values
    v_z = z_processed_x.loc[bool_idx_z]
    v_z = v_z.reset_index(drop=True)
    # if u:
    #     for nu in nu_visted:
    #         bool_idx_u = pd.DataFrame([False] * len(z_processed_x.index), columns=['bool'])
    #         bool_idx_u.index = idx_z_x
    #         bool_idx_u.loc[idx] = True
    #         bool_idx_u = bool_idx_u['bool'].values
    #         v_u[nu] = v_u[nu].loc[bool_idx_u]
    #     v_u = aligned_dict_to_df(v_u)
    #     v_zu = pd.concat([v_z, v_u], axis=1, sort=False)
    # else:
    v_zu = v_z

    if len(v_x.shape) == 1:
        v_x = v_x.to_frame()
    if len(v_y.shape) == 1:
        v_y = v_y.to_frame()
    if len(v_zu.shape) == 1:
        v_zu = v_zu.to_frame()

    cmi_pval, cmi_val = indep_test(v_x, v_y, z=v_zu, sig_samples=sig_samples, p_value=p_value, measure="cmiknn", k=k)
    return cmi_pval, cmi_val


class CTMI:
    def __init__(self, data, mission):
        self.data_dict = dict()
        self.sampling_rate_dict = dict()

        for col in range(data.shape[1]):
            _, s_r = get_sampling_rate(data[data.columns[col]])
            self.sampling_rate_dict[data.columns[col]] = s_r

        alpha = get_alpha(data)
        lags = []

        for col in range(data.shape[1]):
            lags.append(window_size(data[data.columns[col]], alpha=alpha))
            self.data_dict[data.columns[col]] = window_representation(data[data.columns[col]],
                                                                 windows_size=lags[col])

        self.am = align_matrix(self.data_dict, data.columns, self.sampling_rate_dict)

        self.mission = mission

    def fit(self):
        cond_dict = dict()
        cond_names = data.columns[2:]
        for name in cond_names:
            cond_dict[name] = self.data_dict[name]
        cmi = ctmi(self.data_dict[data.columns[0]], self.data_dict[data.columns[1]], cond_dict, data.columns[col1],
                   data.columns[col2], self.sampling_rate_dict, gamma_matrix=self.am, k=10, mission=self.mission)
        return cmi

# def mi(v1, v2, order1=1, order2=1):
#     v1 = ts_order(v1, order1)
#     v2 = ts_order(v2, order2)
#     if len(v1.shape) == 1:
#         v1 = v1.to_frame()
#     if len(v2.shape) == 1:
#         v2 = v2.to_frame()
#     res1 = v1
#     res2 = v2
#
#     cmi_knn = CMIknn(mask_type=None,
#                      significance='shuffle_test',
#                      fixed_thres=None,
#                      sig_samples=10000,
#                      sig_blocklength=3,
#                      knn=10,
#                      confidence='bootstrap',
#                      conf_lev=0.9,
#                      conf_samples=10000,
#                      conf_blocklength=1,
#                      verbosity=0)
#     X = np.concatenate((res1, res2), axis=1)
#     xyz = np.array([0] * order1 + [1] * order2)
#     mi_val = cmi_knn.get_dependence_measure(X.T, xyz)
#     return mi_val
#
#
# def cmi(v1, v2, v3, order1=1, order2=1, order3=1):
#     v1 = ts_order(v1, order1)
#     v2 = ts_order(v2, order2)
#     v3 = ts_order(v3, order3)
#     if len(v1.shape) == 1:
#         v1 = v1.to_frame()
#     if len(v2.shape) == 1:
#         v2 = v2.to_frame()
#     if len(v3.shape) == 1:
#         v3 = v3.to_frame()
#     res1 = v1
#     res2 = v2
#     res3 = v3
#
#     cmi_knn = CMIknn(mask_type=None,
#                      significance='shuffle_test',
#                      fixed_thres=None,
#                      sig_samples=10000,
#                      sig_blocklength=3,
#                      knn=10,
#                      confidence='bootstrap',
#                      conf_lev=0.9,
#                      conf_samples=10000,
#                      conf_blocklength=1,
#                      verbosity=0)
#     X = np.concatenate((res1, res2, res3), axis=1)
#     xyz = np.array([0] * order1 + [1] * order2 + [2] * order3)
#     cmi_val = cmi_knn.get_dependence_measure(X.T, xyz)
#     return cmi_val


if __name__ == "__main__":
    from data.sim_data import generate_fork, generate_v_structure, generate_mediator, generate_diamond
    get_data = {"fork": generate_fork, "v_structure": generate_v_structure, "mediator": generate_mediator,
                "diamond": generate_diamond}

    #######
    data_name = "v_structure"
    scale = False
    #######

    order = 1
    n_samples_list = [100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750, 4000,
                      4250, 4500, 4750, 5000, 10000]
    # n_samples_list = [100, 250, 500, 750, 1000, 1250]

    main_method = "i_ctmi"
    col1 = 0
    col2 = 1
    col3 = 2
    col4 = 2

    output = []
    for n_samples in n_samples_list:
        result = []
        for it in range(100):
            print("iteration: "+str(it))
            data = get_data[data_name](n_samples)
            if scale:
                data -= data.min()
                data /= data.max()

            if main_method == "tmi":
                data_dict1 = dict()
                s_rs1 = []
                s_rs_dict1 = dict()

                lags1 = []
                for col in range(data.shape[1]):
                    _, s_r = get_sampling_rate(data[data.columns[col]])
                    s_rs1.append(s_r)
                    s_rs_dict1[data.columns[col]] = s_r

                a = get_alpha(data)

                for col in range(data.shape[1]):
                    lags1.append(window_size(data[data.columns[col]], alpha=a))
                    data_dict1[data.columns[col]] = window_representation(data[data.columns[col]],
                                                                          windows_size=lags1[col])

                print("------------------------")
                am = align_matrix(data_dict1, data.columns, s_rs_dict1)
                print("alpha: "+str(a))
                print(am)

                data_col1 = data_dict1[data.columns[col1]]
                data_col2 = data_dict1[data.columns[col2]]
                dsr1 = abs(s_rs1[col1] - s_rs1[col2])

                res = tmi(data_col1, data_col2, sampling_rate_tuple=(s_rs1[col1], s_rs1[col2]),
                          gamma=am[data.columns[col2]].loc[data.columns[col1]])
                print("cti: " + str(res))
                result.append(res)
            elif main_method == "ctmi":
                data_dict1 = dict()
                lags1 = []
                s_rs1 = []
                s_rs_dict1 = dict()

                for col in range(data.shape[1]):
                    _, s_r1 = get_sampling_rate(data[data.columns[col]])
                    s_rs1.append(s_r1)
                    s_rs_dict1[data.columns[col]] = s_r1

                a = get_alpha(data)

                for col in range(data.shape[1]):
                    lags1.append(window_size(data[data.columns[col]], alpha=a))
                    data_dict1[data.columns[col]] = window_representation(data[data.columns[col]],
                                                                          windows_size=lags1[col])
                print("lags: "+str(lags1))
                print("sampling rates: "+str(s_rs1))
                print("sampling rates dict: "+str(s_rs_dict1))

                print("------------------------")
                am = align_matrix(data_dict1, data.columns, s_rs_dict1)
                print("gamam matrix: \n"+str(am))

                data_col1 = data_dict1[data.columns[col1]]
                data_col2 = data_dict1[data.columns[col2]]
                data_col3 = data_dict1[data.columns[col3]]
                data_col4 = data_dict1[data.columns[col4]]

                sampling_rate_dict1 = s_rs_dict1

                res = ctmi(data_col1, data_col2, {data.columns[col3]: data_col3}, data.columns[col2],
                           data.columns[col1],
                           sampling_rate_dict1, gamma_matrix=am, k=10)
                print("ccti: " + str(res))
                result.append(res)

            elif main_method == "i_ctmi":
                data_dict1 = dict()
                lags1 = []
                s_rs1 = []
                s_rs_dict1 = dict()

                for col in range(data.shape[1]):
                    _, s_r1 = get_sampling_rate(data[data.columns[col]])
                    s_rs1.append(s_r1)
                    s_rs_dict1[data.columns[col]] = s_r1

                a = get_alpha(data)

                for col in range(data.shape[1]):
                    lags1.append(window_size(data[data.columns[col]], alpha=a))
                    data_dict1[data.columns[col]] = window_representation(data[data.columns[col]],
                                                                          windows_size=lags1[col])
                print("lags: " + str(lags1))
                print("sampling rates: " + str(s_rs1))
                print("sampling rates dict: " + str(s_rs_dict1))

                print("------------------------")
                am = align_matrix(data_dict1, data.columns, s_rs_dict1)
                print("gamam matrix: \n" + str(am))

                data_col1 = data_dict1[data.columns[col1]]
                data_col2 = data_dict1[data.columns[col2]]
                data_col3 = data_dict1[data.columns[col3]]
                data_col4 = data_dict1[data.columns[col4]]

                sampling_rate_dict1 = s_rs_dict1

                res = i_ctmi(data_col1, data_col2, {data.columns[col3]: data_col3}, data.columns[col2],
                           data.columns[col1], sampling_rate_dict1, gamma_matrix=am, k=10)
                print("ccti: " + str(res))
                result.append(res)

        print(result)
        print("result:")
        print("("+str(n_samples)+", "+str(np.mean(result))+") +- ("+str(np.std(result))+", "+str(np.std(result))+")")
        output.append("("+str(n_samples)+", "+str(np.mean(result))+") +- ("+str(np.std(result))+", " +
                      str(np.std(result))+")")
