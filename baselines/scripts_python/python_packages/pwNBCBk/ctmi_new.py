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
import math

from baselines.scripts_python.python_packages.pwNBCBk.tigramite.tigramite.independence_tests import CMIknn
from baselines.scripts_python.python_packages.pwNBCBk.tigramite.tigramite.independence_tests import ParCorr


def indep_test(x, y, z=None, sig_samples=10000, p_value=True, measure="cmiknn", k=10, test_indep= True):
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

    # ws_xy = max(x.shape[1], x.shape[1])
    # ws_xy = 5
    ws_xy = 1
    if z is not None:
        # todo validate
        x_past = x[:-ws_xy]
        y_past = y[:-ws_xy]

        # x_past = x[:-1]
        # y_past = y[:-1]
        # x_past = window_representation(x_past[x_past.columns[0]], ws_xy)
        # y_past = window_representation(y_past[y_past.columns[0]], ws_xy)
        x = x[ws_xy:].reset_index(drop=True)
        y = y[ws_xy:].reset_index(drop=True)
        z = z[ws_xy:].reset_index(drop=True)

        x_past = x_past[x_past.columns[0]]
        y_past = y_past[y_past.columns[0]]

        z = pd.concat([z, x_past, y_past], axis=1)
        # z = pd.concat([z, y_past], axis=1)

        dim_z = z.shape[1]
        X = np.concatenate((x.values, y.values, z.values), axis=1)
        xyz = np.array([0] * dim_x + [1] * dim_y + [2] * dim_z)
    else:
        # todo validate
        x_past = x[:-ws_xy]
        y_past = y[:-ws_xy]

        # x_past = x[:-1]
        # y_past = y[:-1]
        # x_past = window_representation(x_past[x_past.columns[0]], ws_xy)
        # y_past = window_representation(y_past[y_past.columns[0]], ws_xy)
        x = x[ws_xy:].reset_index(drop=True)
        y = y[ws_xy:].reset_index(drop=True)

        x_past = x_past[x_past.columns[0]]
        y_past = y_past[y_past.columns[0]]

        z = pd.concat([x_past, y_past], axis=1)
        # z = y_past.to_frame()
        # print(x)
        # print(x_past)
        # print("ssssssssssssssssssss")
        dim_z = z.shape[1]
        X = np.concatenate((x.values, y.values, z.values), axis=1)
        xyz = np.array([0] * dim_x + [1] * dim_y + [2] * dim_z)

    value = cd.get_dependence_measure(X.T, xyz)
    if p_value:
        pvalue = cd.get_shuffle_significance(X.T, xyz, value)
    else:
        pvalue = value
    return pvalue, value
    # return value


def indep_test_past(x, z=None, sig_samples=10000, p_value=True, measure="cmiknn", k=10, test_indep= True):
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

    ws_xy = 5
    x_past = x[:-ws_xy]
    x = x[ws_xy:].reset_index(drop=True)
    if z is not None:
        z = z[ws_xy:].reset_index(drop=True)

        dim_z = z.shape[1]
        X = np.concatenate((x.values, x_past.values, z.values), axis=1)
        xyz = np.array([0] * dim_x + [1] * dim_y + [2] * dim_z)
    else:

        X = np.concatenate((x.values, y.values), axis=1)
        xyz = np.array([0] * dim_x + [1] * dim_y)

    value = cd.get_dependence_measure(X.T, xyz)
    if p_value:
        pvalue = cd.get_shuffle_significance(X.T, xyz, value)
    else:
        pvalue = value
    return pvalue, value



def indep_test_simple(x, y, z=None, sig_samples=10000, p_value=True, measure="cmiknn", k=10, test_indep= True):
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
        X = np.concatenate((x.values, y.values), axis=1)
        xyz = np.array([0] * dim_x + [1] * dim_y)

    value = cd.get_dependence_measure(X.T, xyz)
    if p_value:
        pvalue = cd.get_shuffle_significance(X.T, xyz, value)
    else:
        pvalue = value
    return pvalue, value
    # return value


def get_shuffle_significance_dependence(ci, array, xyz, value, return_null_dist=False):
    """Returns p-value for shuffle significance test.

    For residual-based test statistics only the residuals are shuffled.

    Parameters
    ----------
    array : array-like
        data array with X, Y, Z in rows and observations in columns

    xyz : array of ints
        XYZ identifier array of shape (dim,).

    value : number
        Value of test statistic for unshuffled estimate.

    Returns
    -------
    pval : float
        p-value
    """
    def get_single_residuals(array, target_var,
                              standardize=True,
                              return_means=False):
        """Returns residuals of linear multiple regression.

        Performs a OLS regression of the variable indexed by target_var on the
        conditions Z. Here array is assumed to contain X and Y as the first two
        rows with the remaining rows (if present) containing the conditions Z.
        Optionally returns the estimated regression line.

        Parameters
        ----------
        array : array-like
            data array with X, Y, Z in rows and observations in columns

        target_var : {0, 1}
            Variable to regress out conditions from.

        standardize : bool, optional (default: True)
            Whether to standardize the array beforehand. Must be used for
            partial correlation.

        return_means : bool, optional (default: False)
            Whether to return the estimated regression line.

        Returns
        -------
        resid [, mean] : array-like
            The residual of the regression and optionally the estimated line.
        """

        dim, T = array.shape
        dim_z = dim - 2

        # Standardize
        if standardize:
            array -= array.mean(axis=1).reshape(dim, 1)
            array /= array.std(axis=1).reshape(dim, 1)
            if np.isnan(array).sum() != 0:
                raise ValueError("nans after standardizing, "
                                 "possibly constant array!")

        y = array[target_var, :]

        if dim_z > 0:
            z = np.fastCopyAndTranspose(array[2:, :])
            beta_hat = np.linalg.lstsq(z, y, rcond=None)[0]
            mean = np.dot(z, beta_hat)
            resid = y - mean
        else:
            resid = y
            mean = None

        if return_means:
            return (resid, mean)
        return resid

    cd = CMIknn(mask_type=None, significance='shuffle_test', fixed_thres=None, sig_samples=10000,
                sig_blocklength=3, knn=10, confidence='bootstrap', conf_lev=0.9, conf_samples=10000,
                conf_blocklength=1, verbosity=0)

    x_vals = get_single_residuals(array, target_var=0)
    y_vals = get_single_residuals(array, target_var=1)
    array_resid = np.array([x_vals, y_vals])
    xyz_resid = np.array([0, 1])

    null_dist = cd._get_shuffle_dist(array_resid, xyz_resid, cd.get_dependence_measure,
                                       sig_samples=10000,
                                       sig_blocklength=3,
                                       verbosity=0)

    # test independence or test dependence
    pval = (null_dist < np.abs(value)).mean()

    # Adjust p-value for two-sided measures
    # if pval < 1.:
    #     pval *= 2.

    if return_null_dist:
        return pval, null_dist
    return pval


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


def align_cpair(x, y, sampling_rate_tuple=(1, 1), k=10, max_gamma=5, set_numbers="Z"):
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

    return g

def align_matrix(data_dict, keys, sampling_rates, k=10, max_gamma=5):
    d = len(keys)
    g_matrix = np.zeros([d, d], dtype=int)
    for i in range(d):
        for j in range(i, d):
            if i != j:
                x = data_dict[keys[i]]
                y = data_dict[keys[j]]
                g = align_cpair(x, y, (sampling_rates[keys[i]], sampling_rates[keys[j]]), k=k, max_gamma=max_gamma)
                g_matrix[i, j] = g
                g_matrix[j, i] = -g
            else:
                g = 1
                g_matrix[i, j] = g

    return pd.DataFrame(g_matrix, columns=keys, index=keys)


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
    mi_list = []
    for i in range(3, lag_max + 1):
        wts = window_representation(ts, windows_size=i)
        i_data = wts[wts.columns[0]].to_frame()
        j_data = wts[wts.columns[i - 1]].to_frame()
        c_data = wts[wts.columns[1:i - 1]]

        # mi_pval, mi_val = tmi(i_data, j_data, k=k, p_value=False)
        mi_pval, mi_val = indep_test_simple(i_data, j_data, z=c_data, sig_samples=10000, p_value=False)
        if mi_val == np.inf:
            mi_val = 1
        mi_list.append(mi_val)

    mi_array = np.array(mi_list)

    j = mi_array.argmax()
    if p_value:
        wts = window_representation(ts, windows_size=j + 3)
        i_data = wts[wts.columns[0]].to_frame()
        j_data = wts[wts.columns[j + 2]].to_frame()
        c_data = wts[wts.columns[1:j + 2]]

        # mi_pval, mi_val = tmi(i_data, j_data, k=k, p_value=True)
        mi_pval, mi_val = indep_test_simple(i_data, j_data, z=c_data, sig_samples=10000, p_value=True)

        if mi_pval <= alpha:
            window = j + 3
        else:
            window = 1
    else:
        if mi_array[j] > alpha:
            window = j + 3
        else:
            window = 1

    if window == 1:
        wts = window_representation(ts, windows_size=3)
        i_data = wts[wts.columns[0]].to_frame()
        j_data = wts[wts.columns[1]].to_frame()
        c_data = wts[wts.columns[2]].to_frame()
        _, mi_val1 = indep_test_simple(i_data, j_data, sig_samples=10000, p_value=False)
        _, mi_val2 = indep_test_simple(i_data, j_data, z=c_data, sig_samples=10000, p_value=False)
        # i_data = wts[wts.columns[2]].to_frame()
        # j_data = wts[wts.columns[1]].to_frame()
        # c_data = wts[wts.columns[0]].to_frame()
        # _, mi_val1 = indep_test(i_data, j_data, z=c_data, sig_samples=10000, p_value=False)

        print(mi_val2, mi_val1)
        if mi_val2 < mi_val1:
            window = 2

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


def get_index_of_aligned_dict(dict_data, keys):
    # check if dict is align
    legnth_list = []
    for name in keys:
        legnth_list.append(dict_data[name].shape[0])
    if legnth_list.count(legnth_list[0]) != len(legnth_list):
        print("Error: time series in dict are not aligned")
        exit(0)
    index_df = pd.DataFrame()
    for name in keys:
        index_df[name] = dict_data[name].index.map(str)
    index_df.insert(0, "concatenated", "id_")
    index_df["concatenated"] = index_df.sum(axis=1)
    index_df = index_df.set_index('concatenated')
    return index_df.index


def aligned_dict_to_df(dict_data):
    concat_df = pd.DataFrame()
    for name in dict_data.keys():
        if isinstance(dict_data[name], pd.Series):
            dict_data[name] = dict_data[name].to_frame()
        concat_df[dict_data[name].columns] = dict_data[name].reset_index(drop=True)
    return concat_df


def find_gamma_x_y(x, y, sampling_rate_tuple=(1, 1), k=10, max_gamma=5):
    c1 = list()
    c2 = list()

    _, val = tmi(x, y, sampling_rate_tuple, k=k, p_value=False)
    c1.append(val)
    c2.append(val)

    for g in range(1, max_gamma):
        _, val = tmi(x, y, sampling_rate_tuple, k=k, gamma=g, p_value=False)
        c1.append(val)

    for g in range(1, max_gamma):
        _, val = tmi(x, y, sampling_rate_tuple, k=k, gamma=-g, p_value=False)
        c2.append(val)
    if np.max(c1) >= np.max(c2):
        g = np.argmax(c1)
    else:
        g = -np.argmax(c2)
    return g


def find_gamma_x(x, sampling_rate_tuple=(1, 1), k=10, max_gamma=5):
    c1 = list()
    c2 = list()

    c = 0
    c1.append(c)
    c2.append(c)

    for g in range(1, max_gamma):
        _, val = tmi(x, x, sampling_rate_tuple, k=k, gamma=g, p_value=False)
        c1.append(val)

    g = np.argmax(c1)
    return g


def gamma_matrix(data_dict, keys, sampling_rates, k=10, max_gamma=5):
    d = len(keys)
    g_matrix = np.zeros([d, d], dtype=int)
    for i in range(d):
        for j in range(i, d):
            if i != j:
                x = data_dict[keys[i]]
                y = data_dict[keys[j]]
                g = find_gamma_x_y(x, y, (sampling_rates[keys[i]], sampling_rates[keys[j]]), k=k, max_gamma=max_gamma)
                g_matrix[i, j] = g
                g_matrix[j, i] = -g
            else:
                x = data_dict[keys[i]]
                g = find_gamma_x(x, (sampling_rates[keys[i]], sampling_rates[keys[j]]), k=k, max_gamma=max_gamma)
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


def find_gamma_lambda_x_y(x, y, sampling_rate_tuple=(1, 1), k=10, max_gamma=5):
    c = np.zeros([2*max_gamma - 1, max_gamma -3, max_gamma - 3])
    gamma_list = list(range(-max_gamma+1, max_gamma))

    ws_x_list = list(range(1, max_gamma-2))
    ws_y_list = list(range(1, max_gamma-2))
    for idx_g in range(len(gamma_list)):
        for idx_ws_x in range(len(ws_x_list)):
            x_w_rep = window_representation(x, windows_size=ws_x_list[idx_ws_x])
            for idx_ws_y in range(len(ws_y_list)):
                if ws_x_list[idx_ws_x] == ws_y_list[idx_ws_y] == 1:
                    y_w_rep = window_representation(y, windows_size=ws_y_list[idx_ws_y])
                    g = gamma_list[idx_g]
                    if g >= 0:
                        _, val = tmi(x_w_rep, y_w_rep, sampling_rate_tuple, k=k, gamma=g, p_value=False)
                    else:
                        val = 0
                    c[idx_g, idx_ws_x, idx_ws_y] = val
                else:
                    if ws_x_list[idx_ws_x] != ws_y_list[idx_ws_y]:
                        y_w_rep = window_representation(y, windows_size=ws_y_list[idx_ws_y])
                        g = gamma_list[idx_g]
                        if g >= 0:
                            _, val = tmi(x_w_rep, y_w_rep, sampling_rate_tuple, k=k, gamma=g, p_value=False)
                        else:
                            val = 0
                        c[idx_g, idx_ws_x, idx_ws_y] = val
                    else:
                        c[idx_g, idx_ws_x, idx_ws_y] = 0

    idx_g, idx_ws_x, idx_ws_y = np.where(c == np.max(c))
    idx_g= idx_g[0]
    idx_ws_x = idx_ws_x[0]
    idx_ws_y = idx_ws_y[0]
    g = gamma_list[idx_g]
    ws_x = ws_x_list[idx_ws_x]
    ws_y = ws_y_list[idx_ws_y]
    return g, ws_x, ws_y


def gamma_matrix_window_matrix(series, keys, sampling_rates, graph, k=10, max_gamma=5):
    d = len(keys)
    g_matrix = np.zeros([d, d], dtype=int)
    window_matrix = np.zeros([d, d], dtype=int)

    for i in range(d):
        for j in range(i, d):
            if i != j:
                x = series[keys[i]]
                y = series[keys[j]]
                if graph[i, j] == 2:
                    g, ws_x, ws_y = find_gamma_lambda_x_y(x, y, (sampling_rates[keys[i]], sampling_rates[keys[j]]), k=k, max_gamma=max_gamma)
                    g_matrix[i, j] = g
                    g_matrix[j, i] = -g
                    window_matrix[i, j] = ws_x
                    window_matrix[j, i] = ws_y
                else:
                    g, ws_y, ws_x = find_gamma_lambda_x_y(y, x, (sampling_rates[keys[j]], sampling_rates[keys[i]]), k=k,
                                                          max_gamma=max_gamma)
                    g_matrix[i, j] = -g
                    g_matrix[j, i] = g
                    window_matrix[i, j] = ws_x
                    window_matrix[j, i] = ws_y
            else:
                # x = series[keys[i]]
                # g = find_gamma_x(x, (sampling_rates[keys[i]], sampling_rates[keys[j]]), k=k, max_gamma=max_gamma)
                g_matrix[i, j] = 1
                window_matrix[i, j] = 1
                window_matrix[j, i] = 1
    return pd.DataFrame(g_matrix, columns=keys, index=keys), pd.DataFrame(window_matrix, columns=keys, index=keys)


def align_xy(x, y, gamma, sampling_rate_tuple):
    sr1, sr2 = sampling_rate_tuple
    dsr = abs(sr1 - sr2)
    iter1 = (sr1 % sr2) * dsr
    iter2 = (sr2 % sr1) * dsr
    idx_x = x.index
    idx_y = y.index
    if gamma > 0:
        y = y[gamma:]
        idx_y = idx_y[gamma:]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        x = x[x.index % (iter1 + 1) == 0]
        y = y[y.index % (iter2 + 1) == 0]
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


def find_gamma_xy_z_util(x, y, z, sampling_rate_tuple, k, gamma, sig_samples=10000, measure="cmiknn"):
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


def find_gamma_xy_z(name_x, name_y, x, y, z, sampling_rate_tuple, k, gamma_matrix, max_gamma=5, measure="cmiknn",
                    instantaneous=True):
    gamma_xy = gamma_matrix[name_y].loc[name_x]
    z = z.loc[x.index[0]:]

    c1 = list()
    c2 = list()

    if instantaneous:
        _, val = find_gamma_xy_z_util(x, y, z, sampling_rate_tuple, k=k, gamma=0, measure=measure)
        c1.append(val)
        c2.append(val)

        for g in range(1, abs(gamma_xy)):
            _, val = find_gamma_xy_z_util(x, y, z, sampling_rate_tuple, k=k, gamma=g, measure=measure)
            c1.append(val)

        for g in range(1, max_gamma):
            _, val = find_gamma_xy_z_util(x, y, z, sampling_rate_tuple, k=k, gamma=-g, measure=measure)
            c2.append(val)
    else:
        c1.append(1)
        c2.append(1)
        if abs(gamma_xy) >= 1:
            for g in range(1, 2):
                _, val = find_gamma_xy_z_util(x, y, z, sampling_rate_tuple, k=k, gamma=g, measure=measure)
                c1.append(val)

        for g in range(1, 2):
            _, val = find_gamma_xy_z_util(x, y, z, sampling_rate_tuple, k=k, gamma=-g, measure=measure)
            c2.append(val)

    if np.min(c1) <= np.min(c2):
        g = np.argmin(c1)
    else:
        g = -np.argmin(c2)

    return g


def align_xyz(name_x, name_y, v_x, v_y, idx_x, idx_y, z, sampling_rate_dict, k, gamma_matrix, instantaneous_dict, graph):
    # names_xy = [*xy.keys()]
    # name_x, name_y = names_xy[0], names_xy[1]
    sampling_rate_tuple = (sampling_rate_dict[name_x], sampling_rate_dict[name_y])

    # v_x, v_y, idx_x, idx_y = align_xy(xy[name_x], xy[name_y], g_xy, sampling_rate_tuple)

    v_x2 = v_x.copy()
    v_y2 = v_y.copy()
    idx_x2 = idx_x.copy()
    idx_y2 = idx_y.copy()

    names_z = [*z.keys()]

    v_z = dict()
    v_z2 = dict()
    nz_visted = []

    # k = 0
    for nz in names_z:
        graphical_optimization = False
        if graph is not None:
            g_xz = gamma_matrix[nz].loc[name_x]
            g_yz = gamma_matrix[nz].loc[name_y]

            graph = pd.DataFrame(graph, columns=gamma_matrix.columns, index=gamma_matrix.columns)
            is_xz_indep = graph[name_x].loc[nz] == 0
            is_yz_indep = graph[name_y].loc[nz] == 0

            if (not is_xz_indep) and (is_yz_indep):
                idx_xy = idx_x
                idx_xy2 = idx_x2
                g = g_xz
                # nz_processed = z[nz]
                graphical_optimization = True
            elif (not is_yz_indep) and (is_yz_indep):
                idx_xy = idx_y
                idx_xy2 = idx_y2
                g = g_yz
                # nz_processed = z[nz]
                graphical_optimization = True

        if not graphical_optimization:
            if idx_x[0] <= idx_y[0]:
                idx_xy = idx_x
                idx_xy2 = idx_x2
            else:
                idx_xy = idx_y
                idx_xy2 = idx_y2

            sampling_rate_tuple_xyz = (max(sampling_rate_tuple), sampling_rate_dict[nz])
            v_x_new = v_x.copy()
            v_y_new = v_y.copy()
            v_x_new.index = idx_xy
            v_y_new.index = idx_xy
            g = find_gamma_xy_z(name_x, name_y, v_x_new, v_y_new, z[nz], sampling_rate_tuple_xyz,
                           k, gamma_matrix, max_gamma=5, measure="cmiknn", instantaneous=instantaneous_dict[nz])
            print("gamma = "+str(g))
        nz_processed = z[nz]

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


def ctmi(x, y, z, name_x, name_y, sampling_rate_dict, gamma_matrix, instantaneous_dict, graph=None, p_value=False, k=10, sig_samples=10000):
    sampling_rate_tuple = (sampling_rate_dict[name_x], sampling_rate_dict[name_y])
    g_xy = gamma_matrix[name_y].loc[name_x]
    v_x, v_y, idx_x, idx_y = align_xy(x, y, g_xy, sampling_rate_tuple)

    if z:
        v_x, v_y, v_z = align_xyz(name_x, name_y, v_x, v_y, idx_x, idx_y, z, sampling_rate_dict, k, gamma_matrix, instantaneous_dict, graph)
        names_z = [*z.keys()]
        if len(names_z) > 0:
            v_z = aligned_dict_to_df(v_z)
        else:
            v_z = None
    else:
        v_z = None

    if len(v_x.shape) == 1:
        v_x = v_x.to_frame()
    if len(v_y.shape) == 1:
        v_y = v_y.to_frame()

    cmi_pval, cmi_val = indep_test(v_x, v_y, z=v_z, sig_samples=sig_samples, p_value=p_value, measure="cmiknn", k=k)
    return cmi_pval, cmi_val

def ctmi_past(x, z, name_x, name_y, sampling_rate_dict, gamma_matrix, instantaneous_dict, graph=None, p_value=False, k=10, sig_samples=10000):
    if z:
        v_x, v_y, v_z = align_xyz(name_x, name_y, v_x, v_y, idx_x, idx_y, z, sampling_rate_dict, k, gamma_matrix, instantaneous_dict, graph)
        names_z = [*z.keys()]
        if len(names_z) > 0:
            v_z = aligned_dict_to_df(v_z)
        else:
            v_z = None
    else:
        v_z = None

    if len(v_x.shape) == 1:
        v_x = v_x.to_frame()
    if len(v_y.shape) == 1:
        v_y = v_y.to_frame()

    cmi_pval, cmi_val = indep_test_past(v_x, v_y, z=v_z, sig_samples=sig_samples, p_value=p_value, measure="cmiknn", k=k)
    return cmi_pval, cmi_val

def find_gamma_xz_yz_util(x, y, z, sampling_rate_tuple, k, gamma_xz, gamma_yz, sig_samples=10000, measure="cmiknn"):
    sr1, sr2 = sampling_rate_tuple
    dsr = abs(sr1 - sr2)
    iter1 = (sr1 % sr2) * dsr
    iter2 = (sr2 % sr1) * dsr
    if (gamma_xz != 0) and (gamma_yz != 0):
        gamma = max(gamma_xz, gamma_yz)
        new_gamma_xz = gamma - gamma_xz
        new_gamma_yz = gamma - gamma_yz
        z = z[gamma:]
        if new_gamma_xz != 0:
            x = x[new_gamma_xz:]
        else:
            y = y[new_gamma_yz:]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        z = z.reset_index(drop=True)
        x = x[x.index % (iter1 + 1) == 0]
        y = y[y.index % (iter1 + 1) == 0]
        z = z[z.index % (iter2 + 1) == 0]

        x = x[:-gamma_xz]
        y = y[:-gamma_yz]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        z = z.reset_index(drop=True)
    else:
        if gamma_xz != 0:
            z = z[gamma_xz:]
            x = x.reset_index(drop=True)
            y = y.reset_index(drop=True)
            z = z.reset_index(drop=True)
            x = x[x.index % (iter1 + 1) == 0]
            y = y[y.index % (iter1 + 1) == 0]
            z = z[z.index % (iter2 + 1) == 0]

            x = x[:-gamma_xz]
            y = y[gamma_xz:]
            x = x.reset_index(drop=True)
            y = y.reset_index(drop=True)
            z = z.reset_index(drop=True)
        elif gamma_yz != 0:
            z = z[gamma_yz:]
            x = x.reset_index(drop=True)
            y = y.reset_index(drop=True)
            z = z.reset_index(drop=True)
            x = x[x.index % (iter1 + 1) == 0]
            y = y[y.index % (iter1 + 1) == 0]
            z = z[z.index % (iter2 + 1) == 0]

            x = x[gamma_yz:]
            y = y[:-gamma_yz]
            x = x.reset_index(drop=True)
            y = y.reset_index(drop=True)
            z = z.reset_index(drop=True)
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


def find_gamma_xz_yz(x, y, z, sampling_rate_tuple, k, max_gamma=5, measure="cmiknn"):
    z = z.loc[x.index[0]:]

    val_list = list()
    g_xz_list = list()
    g_yz_list = list()
    _, val = find_gamma_xz_yz_util(x, y, z, sampling_rate_tuple, k=k, gamma_xz=0, gamma_yz=0, measure=measure)

    for g_xz in range(max_gamma):
        for g_yz in range(max_gamma):
            _, val = find_gamma_xz_yz_util(x, y, z, sampling_rate_tuple, k=k, gamma_xz=g_xz, gamma_yz=g_yz, measure=measure)
            val_list.append(val)
            g_xz_list.append(g_xz)
            g_yz_list.append(g_yz)

    idx = np.argmax(val_list)
    g_xz = g_xz_list[idx]
    g_yz = g_yz_list[idx]
    return g_xz, g_yz


def i_ctmi(x, y, z, name_x, name_y, name_z, sampling_rate_dict, p_value=True, k=10, sig_samples=10000, test_indep=False):

    v_x = x.copy()
    v_y = y.copy()
    v_z = z.copy()

    sampling_rate_tuple_xy = (sampling_rate_dict[name_x], sampling_rate_dict[name_y])
    g_xz, g_yz = find_gamma_xz_yz(x, y, z, sampling_rate_tuple_xy, k, max_gamma=5, measure="cmiknn")

    xz_dict = {name_x: v_x, name_z: v_z}
    yz_dict = {name_y: v_y, name_z: v_z}

    sampling_rate_tuple_x = (sampling_rate_dict[name_x], sampling_rate_dict[name_z])
    sampling_rate_tuple_y = (sampling_rate_dict[name_y], sampling_rate_dict[name_z])

    v_x, z_processed_x, idx_x, idx_z_x = align_xy(xz_dict[name_x], xz_dict[name_z], g_xz, sampling_rate_tuple_x)
    v_y, z_processed_y, idx_y, idx_z_y = align_xy(yz_dict[name_y], yz_dict[name_z], g_yz, sampling_rate_tuple_y)

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

    bool_idx_z = pd.DataFrame([False] * len(z_processed_x.index), columns=['bool'])
    bool_idx_z.index = idx_z_x
    bool_idx_z.loc[idx] = True
    bool_idx_z = bool_idx_z['bool'].values
    v_z = z_processed_x.loc[bool_idx_z]
    v_z = v_z.reset_index(drop=True)

    if len(v_x.shape) == 1:
        v_x = v_x.to_frame()
    if len(v_y.shape) == 1:
        v_y = v_y.to_frame()
    if len(v_z.shape) == 1:
        v_z = v_z.to_frame()

    cmi_pval, cmi_val = indep_test(v_x, v_y, z=v_z, sig_samples=sig_samples, p_value=p_value, measure="cmiknn", k=k,
                                   test_indep=test_indep)
    return cmi_pval, cmi_val


if __name__ == "__main__":
    path = "../../../../data/simulated_ts_data/v_structure/data_"+str(0)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)
    data = data.loc[:1000]

    data_col0 = data[data.columns[0]]
    data_col1 = data[data.columns[1]]
    data_col2 = data[data.columns[2]]

    sampling_rate_dict1 = dict()
    for col in range(data.shape[1]):
        _, s_r1 = get_sampling_rate(data[data.columns[col]])
        sampling_rate_dict1[data.columns[col]] = s_r1

    Rando =  np.random.normal(0, 1, 1000)
    data_col3 = pd.DataFrame()
    N = 1000
    data_col3["Location"] = np.random.choice(Rando, size=N)

    res0 = tmi(data_col0, data_col1, sampling_rate_tuple=(1, 1), k=10, gamma=0, p_value=True, sig_samples=10000)
    print(res0)
    res = i_ctmi(data_col0, data_col1, data_col2, data.columns[0],
                 data.columns[1], data.columns[2], sampling_rate_dict1, k=10)
    print(11)
    print(res)

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
                am = gamma_matrix(data_dict1, data.columns, s_rs_dict1)
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
                am = gamma_matrix(data_dict1, data.columns, s_rs_dict1)
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
                am = gamma_matrix(data_dict1, data.columns, s_rs_dict1)
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
