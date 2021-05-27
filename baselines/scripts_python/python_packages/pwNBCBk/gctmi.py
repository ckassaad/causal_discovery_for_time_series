#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Graphical conditional mutual information between time series: script implementing
the GCTMI methods.

Date: Aug 2020
Author: Karim Assaad, karimassaad3@gmail.com, karim.assaad@univ.grenoble.alpes.fr, karim.assaad@coservit.com
"""

import numpy as np
import pandas as pd

from ctmi import align_xy, align_cond, get_index_of_aligned_dict, aligned_dict_to_df, indep_test


def g_align_xyz(xy, z, sampling_rate_dict, gamma_matrix, graph, mission):
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
    graph = pd.DataFrame(graph, columns=gamma_matrix.columns, index=gamma_matrix.columns)
    for nz in names_z:
        g_xz = gamma_matrix[nz].loc[name_x]
        g_yz = gamma_matrix[nz].loc[name_y]

        is_xz_indep = graph[name_x].loc[nz] == 0
        is_yz_indep = graph[name_y].loc[nz] == 0
        if (not is_xz_indep) and (not is_yz_indep):
            if idx_x[0] <= idx_y[0]:
                idx_xy = idx_x
                idx_xy2 = idx_x2
            else:
                idx_xy = idx_y
                idx_xy2 = idx_y2

            if g_xz == g_yz:
                g = g_xz
            else:
                sampling_rate_tuple_xyz = (max(sampling_rate_tuple), sampling_rate_dict[nz])
                v_x_new = v_x.copy()
                v_y_new = v_y.copy()
                v_x_new.index = idx_xy
                v_y_new.index = idx_xy
                g = align_cond(v_x_new, v_y_new, z[nz], sampling_rate_tuple_xyz,
                               k, max_gamma=5, measure="cmiknn", mission=mission)
            nz_processed = z[nz]

        elif not is_xz_indep:
            idx_xy = idx_x
            idx_xy2 = idx_x2
            g = g_xz
            nz_processed = z[nz]
        elif not is_yz_indep:
            idx_xy = idx_y
            idx_xy2 = idx_y2
            g = g_yz
            nz_processed = z[nz]
        else:
            g = None
            nz_processed = None
            idx_xy = None
            idx_xy2 = None

        if g is not None:
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


def gctmi(x, y, z, name_x, name_y, sampling_rate_dict, gamma_matrix, graph, p_value=False, k=10, sig_samples=10000,
         mission="ci"):
    xy = {name_x: x, name_y: y}
    v_x, v_y, v_z = g_align_xyz(xy, z, sampling_rate_dict, gamma_matrix, graph, mission)

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
