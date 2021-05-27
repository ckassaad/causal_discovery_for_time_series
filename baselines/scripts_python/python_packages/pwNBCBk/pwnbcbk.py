#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Discovery of summary causal graphs for time series: script implementing
the PCTMI and FCITMI methods.
Parallelization is done across variables for the skeleton construction step and for the rule origin of causality.

Date: Dec 2019
Author: Karim Assaad, karimassaad3@gmail.com, karim.assaad@univ.grenoble.alpes.fr, karim.assaad@coservit.com
paper: soon
"""

import numpy as np
import pandas as pd
import itertools
from joblib import Parallel, delayed

from datetime import datetime

from baselines.scripts_python.python_packages.pwNBCBk.ctmi import window_representation, get_sampling_rate, align_matrix, tmi, get_alpha, window_size, align_pair
from baselines.scripts_python.python_packages.pwNBCBk.ctmi_new import i_ctmi, ctmi, align_matrix, tmi, window_size, gamma_matrix_window_matrix
# from gctmi import gctmi


from baselines.scripts_python.python_packages.pwNBCBk.kitmi import nbcb_k



# todo: make ctmi use graph by adding None to cell in gamma matrix which is associated to two independent time series
# todo: make ranking based on statistics (not pvalue)


class Graph:
    """
    Graph structure
    0: no edge
    1: a tail -
    2: arrow head ->
    """
    def __init__(self, d):
        """
        :param d: number of nodes
        """
        self.d = d
        # self.edges = np.subtract(np.ones([n, n]), np.eye(n))
        self.edges = np.ones([d, d])
        self.sep = np.zeros([d, d, d])

    def del_edge(self, p, q):
        """
        :param p: index of a time series
        :param q: index of a time series
        """
        self.edges[p, q] = 0
        self.edges[q, p] = 0

    def add_sep(self, p, q, r):
        """
        :param p: index of a time series
        :param q: index of a time series
        :param r: index of seperation set
        """
        self.sep[p, q, r] = 1
        self.sep[q, p, r] = 1

    def search_adj(self, p):
        """
        :param p: index of a time series
        :return: list of adjacencies of time series p and the number of adjacencies
        """
        adj_1 = np.argwhere(self.edges[p, :] != 0)
        adj_2 = np.argwhere(self.edges[:, p] != 0)
        adj = np.intersect1d(adj_1, adj_2)
        if self.edges[p, p] == 1:
            adj = adj[adj != p]
        num_adj = len(adj)
        return adj, num_adj

    def search_adj_all(self):
        """
        :return: list of adjacencies of all time series and the number of adjacencies per time series
        """
        l_num_adj = []
        l_adj = []
        for p in range(self.d):
            adj, num_adj = self.search_adj(p)
            l_adj.append(adj.tolist())
            l_num_adj.append(num_adj)
        return l_adj, l_num_adj


    def search_par(self, p):
        """
        :param p: index of a time series
        :return: list of adjacencies of time series p and the number of adjacencies
        """
        # adj_1 = np.argwhere(self.edges[p, :] != 0)
        par = np.argwhere(self.edges[:, p] == 2)
        # adj = np.intersect1d(adj_1, adj_2)
        # if self.edges[p, p] == 1:
        #     adj = adj[adj != p]
        num_par = len(par)
        return par, num_par

    def search_par_all(self):
        """
        :return: list of adjacencies of all time series and the number of adjacencies per time series
        """
        l_num_par = []
        l_par = []
        for p in range(self.d):
            par, num_par = self.search_par(p)
            l_par.append(par.tolist())
            l_num_par.append(num_par)
        return l_par, l_num_par


class RankingList:
    def __init__(self):
        self.val = np.array([])
        self.elem_p = np.array([], dtype='int')
        self.elem_q = np.array([], dtype='int')
        self.elem_r = []

    def add(self, p, q, val, r):
        """
        :param p: index of a time series
        :param q: index of a time series
        :param val: value of mutual information
        :param r: index of set of conditionals
        """
        self.val = np.append(self.val, val)
        self.elem_p = np.append(self.elem_p, p)
        self.elem_q = np.append(self.elem_q, q)
        self.elem_r.append(r)

    def sort(self, descending=True):
        """
        :param descending: (bool) sort ascending vs. descending. By default True
        """
        idx = np.argsort(self.val)
        if descending:
            idx = np.flip(idx)
        # self.val = self.val[idx]
        # self.elem_p = self.elem_p[idx]
        # self.elem_q = self.elem_q[idx]
        # self.elem_r = self.elem_r[idx]
        self.val = np.take_along_axis(self.val, idx, axis=0)
        self.elem_p = np.take_along_axis(self.elem_p, idx, axis=0)
        self.elem_q = np.take_along_axis(self.elem_q, idx, axis=0)
        sorted_elem_r = []
        for i in idx:
            sorted_elem_r.append(self.elem_r[i])
        self.elem_r = sorted_elem_r


class pwNBCBk:
    def __init__(self, series, sig_lev=0.05, lag_max=5, p_value=True, rank_using_p_value=False, verbose=True, num_processor=-1,
                 graphical_optimization=False):
        """
        Causal inference (Wrapper) using TMI and CTMI (contain functions for skeleton construction)
        :param series: d-time series (with possibility of different sampling rate)
        :param sig_lev: significance level. By default 0.05
        :param p_value: Use p_value for decision making. By default True
        :param verbose: Print results. By default: True
        :param num_processor: number of processors for parallelization. By default -1 (all)
        """
        self.graph = Graph(series.shape[1])
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

        training_epoch = 1000
        noise = True  # d*(order-1)*2
        learning_rate = 0.01
        for i in range(series.shape[1]):
            for j in range(i+1, series.shape[1]):
                data_pair = series[[series.columns[i], series.columns[j]]]
                # res_order_pair = tskiko_mv(data_pair, lag_max, learning_rate, training_epoch, noise, sig_lev, "ParCorr", verbose)
                res_order_pair = nbcb_k(data_pair, lag_max, learning_rate, training_epoch, noise, sig_lev, "ParCorr", verbose)
                # res_order_pair = run_timino_pw_R([[data_pair, "data"], [0.00, "alpha"], [5, "nlags"]])
                # res_order_pair = pd.DataFrame(res_order_pair, columns=data_pair.columns, index=data_pair.columns)
                if res_order_pair[series.columns[j]].loc[series.columns[i]] == 2:
                    self.graph.edges[i, j] = 2
                if res_order_pair[series.columns[i]].loc[series.columns[j]] == 2:
                    self.graph.edges[j, i] = 2
        # order_kiko = tskiko_mv(series, lag_max, learning_rate, training_epoch, noise, sig_lev, "ParCorr", verbose)
        # print(order_kiko)
        # self.graph.edges = order_kiko.values

        if verbose:
            print("Order")
            print(self.graph.edges)

        self.adaptive_window = True
        self.series = series
        self.n = series.shape[0]
        self.d = series.shape[1]
        self.names = self.series.columns
        self.num_processor = num_processor
        self.p_value = p_value
        self.graphical_optimization = graphical_optimization
        if self.p_value == rank_using_p_value:
            self.rank_using_p_value = rank_using_p_value
        elif not rank_using_p_value:
            self.rank_using_p_value = rank_using_p_value
        else:
            print("Warning: rank_using_p_value can be True iff p_value is True. Using rank_using_p_value=False")
            self.rank_using_p_value = False
        self.verbose = verbose

        self.data_dict = dict()
        self.instantaneous_dict = dict()

        self.lags = []
        self.sampling_rate = dict()
        for col in range(series.shape[1]):
            _, s_r = get_sampling_rate(self.series[self.names[col]])
            self.sampling_rate[self.names[col]] = s_r

        self.sig_lev = sig_lev
        self.alpha = get_alpha(series)

        for col in range(series.shape[1]):
            # self.lags.append(window_size(series[series.columns[col]], alpha=self.alpha, lag_max=lag_max))
            if not self.adaptive_window:
                self.lags.append(1)
                self.data_dict[self.names[col]] = window_representation(self.series[self.names[col]],
                                                                        windows_size=self.lags[col])
            self.instantaneous_dict[self.names[col]] = True

        if self.adaptive_window:
            self.gamma_matrix, self.window_matrix = gamma_matrix_window_matrix(self.series, series.columns, self.sampling_rate, self.graph.edges)
        else:
            self.gamma_matrix = align_matrix(self.data_dict, series.columns, self.sampling_rate)

        self.cap_gamma_df = pd.DataFrame(columns=["p", "q", "r", "Grp", "Grq"])

        self.mi_array = np.ones([self.graph.d, self.graph.d])
        self.cmi_array = np.ones([self.graph.d, self.graph.d])

        self.biggamma = np.zeros([self.d, self.d, self.d])


        if self.verbose:
            print("n: "+str(self.n))
            print("d: "+str(self.d))
            print("names: "+str(self.names))
            print("sampling_rate: "+str(self.sampling_rate))
            print("significance level:"+str(self.sig_lev))
            print("alpha:"+str(self.alpha))
            print("window size:"+str(self.lags))
            print("gamma matrix:"+str(self.gamma_matrix))
            if self.adaptive_window:
                print("window matrix"+str(self.window_matrix))
            print("instantaneous dict :"+str(self.instantaneous_dict))


    def _mi_pq(self, p, q):
        """
        estimate tmi between two time series
        :param p: time series with index p
        :param q: time series with index q
        :return: p, q and the estimated value of tmi(p,q)
        """
        if self.adaptive_window:
            x = window_representation(self.series[self.names[p]], windows_size=self.window_matrix[self.names[q]].loc[self.names[p]])
            y = window_representation(self.series[self.names[q]], windows_size=self.window_matrix[self.names[p]].loc[self.names[q]])
            print("Nodes and windows:")
            print(self.names[p], self.window_matrix[self.names[q]].loc[self.names[p]])
            print(self.names[q], self.window_matrix[self.names[p]].loc[self.names[q]])
        else:
            x = self.data_dict[self.names[p]]
            y = self.data_dict[self.names[q]]

        mi_pval, mi_val = tmi(x, y, sampling_rate_tuple=(self.sampling_rate[self.names[p]],
                                                         self.sampling_rate[self.names[q]]),
                              gamma=self.gamma_matrix[self.names[q]].loc[self.names[p]], p_value=self.p_value)
        # mi_pval, mi_val = ctmi(x, y, None, self.names[p], self.names[q], self.sampling_rate,
        #                          gamma_matrix=self.gamma_matrix, p_value=self.rank_using_p_value)
        return p, q, mi_pval

    def skeleton_initialize(self):
        """
        initialize graph, remove all unconditional independencies and rank neighbors
        """
        if self.verbose:
            print("######################################")
            print("Skeletion Initialization")
            print("######################################")

        # p_list, q_list = np.where(np.triu(self.graph.edges) > 0)
        p_list, q_list = np.where((np.triu(self.graph.edges)-np.diag(np.diag(self.graph.edges))) > 0)
        res = Parallel(n_jobs=self.num_processor)(delayed(self._mi_pq)(p, q) for p, q in zip(p_list, q_list))

        for pq in range(len(res)):
            p, q, mi = res[pq][0], res[pq][1], res[pq][2]
            self.mi_array[p, q] = mi
            self.mi_array[q, p] = mi
            if self.verbose:
                print("p=" + str(p) + "; q=" + str(q) + "; I(p,q)=" + "{: 0.5f}".format(self.mi_array[p, q]), end=" ")
            if self.p_value:
                test = self.mi_array[p, q] > self.sig_lev
            else:
                test = self.mi_array[p, q] < self.alpha
            if test:
                if self.verbose:
                    print("=> Remove link between "+str(p)+" and "+str(q))
                self.graph.edges[p, q] = 0
                self.graph.edges[q, p] = 0
            else:
                if self.verbose:
                    print()

    def _cmi_sep_set_pq(self, p, q, set_size):
        """
        estimate ctmi between two time series conditioned on each set of neighbors with cardinality equal to set_size
        :param p: time series with index p
        :param q: time series with index q
        :param set_size: cardinality of the set of neighbors
        :return: p, q, list if estimated value of ctmi(p,q,r_set), and list of all r_sets
        """
        v_list = []
        # r_list = [r for r in range(self.graph.d) if (r != p) and (r != q) and ((
        #         (self.graph.edges[p, r] != 0) and (self.gamma_matrix[self.names[p]].loc[self.names[r]] >= 0)) or (
        #         (self.graph.edges[q, r] != 0) and (self.gamma_matrix[self.names[q]].loc[self.names[r]] >= 0)))]
        r_list = [r for r in range(self.graph.d) if (r != p) and (r != q) and (self.graph.edges[r, q] == 2)]

        r_list = [list(r) for r in itertools.combinations(r_list, set_size)]

        r_list_temp = r_list.copy()
        # if set_size == 1:
        for rs in r_list_temp:
            print(rs)
            print(all(elem >= self.d for elem in rs))
            if all(elem >= self.d for elem in rs):
                r_list.remove(rs)
        del r_list_temp

        if self.adaptive_window:
            x = window_representation(self.series[self.names[p]], windows_size=int(self.window_matrix[self.names[q]].loc[self.names[p]]))
            y = window_representation(self.series[self.names[q]], windows_size=int(self.window_matrix[self.names[p]].loc[self.names[q]]))
        else:
            x = self.data_dict[self.names[p]]
            y = self.data_dict[self.names[q]]

        for rs in r_list:
            z = dict()
            for r in rs:
                if self.adaptive_window:
                    # select and drop NA
                    z[self.names[r]] = self.series[self.names[r]].dropna()
                else:
                    z[self.names[r]] = self.data_dict[self.names[r]]
            if self.graphical_optimization:
                # cmi_pval, cmi_val = gctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
                #                           gamma_matrix=self.gamma_matrix, p_value=self.rank_using_p_value,
                #                           graph=self.graph.edges)
                cmi_pval, cmi_val = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
                                          gamma_matrix=self.gamma_matrix, graph=self.graph.edges,
                                         p_value=self.rank_using_p_value, instantaneous_dict=self.instantaneous_dict)
            else:
                cmi_pval, cmi_val = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
                                         gamma_matrix=self.gamma_matrix, p_value=self.rank_using_p_value,
                                         instantaneous_dict=self.instantaneous_dict)

            if self.rank_using_p_value:
                v_list.append(cmi_pval)
            else:
                v_list.append(cmi_val)
        if v_list:
            return p, q, v_list, r_list

    def rank_cmi_sep_set_parallel(self, set_size):
        """
        rank pairs of time series based on the estimation of ctmi between each pair of connected time series
        :param set_size: cardinality of the set of neighbors
        :return: ranking of each pair of connected time series based ctmi
        """
        list_par, list_num_par = self.graph.search_par_all()
        q_list = [q for q in range(len(list_num_par)) if list_num_par[q] > set_size]
        p_list = [list_par[q] for q in q_list]
        q_list = [q_list[q] for q in range(len(q_list)) for _ in p_list[q]]
        p_list = [p for sublist in p_list for p in sublist]
        p_list = [p for sublist in p_list for p in sublist]
        # pq_list = [(p, q) for p, q in zip(p_list, q_list)]
        # temp_pq = pq_list.copy()
        # temp_p = p_list.copy()
        # temp_q = q_list.copy()
        # for pq in range(len(temp_pq)):
        #     if (temp_pq[pq][1], temp_pq[pq][0]) in pq_list:
        #         pq_list.remove((temp_pq[pq][0], temp_pq[pq][1]))
        #         p_list.remove(temp_p[pq])
        #         q_list.remove(temp_q[pq])
        # del temp_pq, temp_p, temp_q
        print(list_par, list_num_par)
        print(p_list, q_list)
        print("set_size " +str(set_size))
        # res = Parallel(n_jobs=self.num_processor)(delayed(self._cmi_sep_set_pq)(p, q, set_size) for p, q in
        #                                           zip(p_list, q_list))
        res = []
        for p, q in zip(p_list, q_list):
            res.append(self._cmi_sep_set_pq(p, q, set_size))

        ranks = RankingList()
        for pq in range(len(res)):
            if res[pq] is not None:
                if isinstance(res[pq][2], list):
                    for r in range(len(res[pq][2])):
                        ranks.add(res[pq][0], res[pq][1], res[pq][2][r], res[pq][3][r])
                else:
                    ranks.add(res[pq][0], res[pq][1], res[pq][2], res[pq][3])
        if self.rank_using_p_value:
            ranks.sort(descending=True)
        else:
            ranks.sort(descending=False)
        return ranks

    def find_sep_set(self):
        """
        find the most contributing separation set (if it exists) between each pair of time series
        """
        if self.verbose:
            print("######################################")
            print("Skeletion Speperation")
            print("######################################")

        print("max set size = " + str(self.graph.d-1))
        for set_size in range(1, self.graph.d-1):
            ranks = self.rank_cmi_sep_set_parallel(set_size)
            if self.verbose:
                print("Ranking:")
                print("p: "+str(ranks.elem_p))
                print("p: " + str(ranks.elem_q))
                print("p: " + str(ranks.elem_r))
                print("p: " + str(ranks.val))
            for p, q, r_set, cmi in zip(ranks.elem_p, ranks.elem_q, ranks.elem_r, ranks.val):
                test = (self.graph.edges[p, q] != 0)
                for r in r_set:
                    if not test:
                        break
                    test = test and ((self.graph.edges[q, r] != 0) or (self.graph.edges[p, r] != 0))
                    # test = test and ((self.graph.sep[p, r, q] == 0) and (self.graph.sep[q, r, p] == 0))
                if test:
                    mi = self.mi_array[p, q]

                    if self.p_value != self.rank_using_p_value:
                        if self.adaptive_window:
                            x = window_representation(self.series[self.names[p]],
                                                      windows_size=self.window_matrix[self.names[q]].loc[self.names[p]])
                            y = window_representation(self.series[self.names[q]],
                                                      windows_size=self.window_matrix[self.names[p]].loc[self.names[q]])
                        else:
                            x = self.data_dict[self.names[p]]
                            y = self.data_dict[self.names[q]]

                        z = dict()
                        for r in r_set:
                            if self.adaptive_window:
                                # select and drop NA
                                z[self.names[r]] = self.series[self.names[r]].dropna()
                            else:
                                z[self.names[r]] = self.data_dict[self.names[r]]
                        if self.graphical_optimization:
                            # cmi, _ = gctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
                            #                gamma_matrix=self.gamma_matrix, p_value=self.p_value, graph=self.graph.edges)
                            cmi_pval, cmi_val = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
                                                     gamma_matrix=self.gamma_matrix, graph=self.graph.edges,
                                                     p_value=self.p_value,
                                                     instantaneous_dict=self.instantaneous_dict)
                        else:
                            cmi, _ = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
                                          gamma_matrix=self.gamma_matrix, p_value=self.p_value,
                                          instantaneous_dict=self.instantaneous_dict)
                    if self.verbose:
                        print("p=" + str(p) + "; q=" + str(q) + "; r=" + str(r_set) + "; I(p,q|r)=" + "{: 0.5f}".format(
                            cmi) + "; I(p,q)=" + "{: 0.5f}".format(mi), end=" ")

                    if self.p_value:
                        test = mi < self.sig_lev < cmi
                    else:
                        test = cmi < self.alpha
                    if test:
                        self.cmi_array[p, q] = cmi
                        self.cmi_array[q, p] = cmi
                        if self.verbose:
                            print("=> remove link between " + str(p) + " and " + str(q))
                        self.graph.edges[p, q] = 0
                        self.graph.edges[q, p] = 0

                        for r in r_set:
                            self.graph.add_sep(q, p, r)
                            self.biggamma[p,q,r] = self.gamma_matrix[self.names[p]].loc[self.names[r]]
                            self.biggamma[q,p,r] = self.gamma_matrix[self.names[q]].loc[self.names[r]]
                    else:
                        if self.verbose:
                            print()
        # self._exclude_past()


    def fit(self):
        """
        run ACITMI
        :return: graph (CPDAG)
        """
        if self.verbose:
            now = datetime.now()
            print("#######################################")
            print("########### Starting ACITMI ###########")
            print("########### " + now.strftime("%H:%M:%S" + " ###########"))
            print("#######################################")

        # Progressive Removal of Non-Causal Nodes
        self.skeleton_initialize()
        self.find_sep_set()

        if self.verbose:
            print("######################################")
            print("Final Results (KITMI)")
            print("######################################")
            print("Summary Graph:")
            print(self.graph.edges)
        return self.graph.edges


if __name__ == "__main__":
    import pandas as pd
    path = "../../../../data/simulated_ts_data/v_structure/data_"+str(0)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)
    data = data.loc[:1000]
    print(data)
    import random
    random.seed(2)

