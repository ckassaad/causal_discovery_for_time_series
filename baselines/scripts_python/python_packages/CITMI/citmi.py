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

from baselines.scripts_python.python_packages.CITMI.ctmi import window_representation, get_sampling_rate, align_matrix, tmi, get_alpha, window_size, align_pair
from baselines.scripts_python.python_packages.CITMI.ctmi_new import i_ctmi, ctmi, align_matrix, tmi, window_size, gamma_matrix_window_matrix, self_ctmi
# from gctmi import gctmi

from baselines.scripts_python.python_packages.CITMI.ctmi_new import indep_test


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


class CITMI:
    def __init__(self, series, sig_lev=0.05, lag_max=5, p_value=True, rank_using_p_value=False, verbose=True, num_processor=-1,
                 graphical_optimization=True):
        """
        Causal inference (Wrapper) using TMI and CTMI (contain functions for skeleton construction)
        :param series: d-time series (with possibility of different sampling rate)
        :param sig_lev: significance level. By default 0.05
        :param p_value: Use p_value for decision making. By default True
        :param verbose: Print results. By default: True
        :param num_processor: number of processors for parallelization. By default -1 (all)
        """
        self.adaptive_window = True
        self.series = series
        self.graph = Graph(series.shape[1])
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
            self.gamma_matrix, self.window_matrix = gamma_matrix_window_matrix(self.series, series.columns, self.sampling_rate)
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

        # self._include_past(series)

    def _include_past(self, series):
        self.graph.d = 2*self.graph.d
        self.graph.sep = np.zeros([2*self.d, 2*self.d, 2*self.d])

        self.graph.edges= np.append(self.graph.edges, np.zeros((self.d, self.d)), axis=1)
        self.graph.edges = np.append(self.graph.edges, np.zeros((self.d, self.d*2)), axis=0)

        self.mi_array = np.append(self.mi_array, np.zeros((self.d, self.d)), axis=1)
        self.mi_array = np.append(self.mi_array, np.zeros((self.d, self.d*2)), axis=0)
        self.cmi_array = np.append(self.cmi_array, np.zeros((self.d, self.d)), axis=1)
        self.cmi_array = np.append(self.cmi_array, np.zeros((self.d, self.d*2)), axis=0)

        counter = 0
        for i in range(self.d):
            self.instantaneous_dict[self.names[i] + "_past"] = False
            if self.gamma_matrix[self.names[i]].loc[self.names[i]] > 0:
                self.names = self.names.insert(self.d + counter, self.names[i] + "_past")

                counter = counter + 1
                self.gamma_matrix[self.names[i] + "_past"] = 0
                self.gamma_matrix = self.gamma_matrix.append(pd.Series([0] * (self.d + counter),
                                                                       name=self.names[i] + "_past",
                                                                       index=self.gamma_matrix.columns))

                self.gamma_matrix[self.names[i]].loc[self.names[i] + "_past"] = self.gamma_matrix[self.names[i]].loc[self.names[i]]
                self.gamma_matrix[self.names[i] + "_past"].loc[self.names[i]] = - self.gamma_matrix[self.names[i]].loc[self.names[i]]

                self.graph.edges[self.d + i, i] = 2
                self.graph.edges[i, self.d + i] = 1

                # self.data_dict[self.names[i] + "_past"] = self.data_dict[self.names[i]]
                self.data_dict[self.names[i] + "_past"] = series[self.names[i]]
                self.sampling_rate[self.names[i] + "_past"] = self.sampling_rate[self.names[i]]
        print(self.gamma_matrix)
        print(self.graph.edges)
        print(self.names)

    def _exclude_past(self):
        for i in range(self.d):
            if self.graph.edges[self.d + i, i] == 0:
                self.graph.edges[i, i] = 0
        self.graph.edges = self.graph.edges[:self.d, :self.d]
        self.graph.d = self.d

    def _mi_pq(self, p, q):
        """
        estimate tmi between two time series
        :param p: time series with index p
        :param q: time series with index q
        :return: p, q and the estimated value of tmi(p,q)
        """
        if self.adaptive_window:
            x = window_representation(self.series[self.names[p]], windows_size=self.window_matrix[self.names[p]].loc[self.names[p]])
            y = window_representation(self.series[self.names[q]], windows_size=self.window_matrix[self.names[q]].loc[self.names[q]])
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
        r_list = [r for r in range(self.graph.d) if (r != p) and (r != q) and ((
                (self.graph.edges[p, r] != 0) and (self.gamma_matrix[self.names[p]].loc[self.names[r]] >= 0)) or (
                (self.graph.edges[q, r] != 0) and (self.gamma_matrix[self.names[q]].loc[self.names[r]] >= 0)))]

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
            x = window_representation(self.series[self.names[p]], windows_size=self.window_matrix[self.names[p]].loc[self.names[p]])
            y = window_representation(self.series[self.names[q]], windows_size=self.window_matrix[self.names[q]].loc[self.names[q]])
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
        list_adj, list_num_adj = self.graph.search_adj_all()
        p_list = [p for p in range(len(list_num_adj)) if list_num_adj[p] > set_size]
        print(p_list)
        q_list = [list_adj[p] for p in p_list]
        p_list = [p_list[p] for p in range(len(p_list)) for _ in q_list[p]]
        q_list = [q for sublist in q_list for q in sublist]
        pq_list = [(p, q) for p, q in zip(p_list, q_list)]
        temp_pq = pq_list.copy()
        temp_p = p_list.copy()
        temp_q = q_list.copy()
        for pq in range(len(temp_pq)):
            if (temp_pq[pq][1], temp_pq[pq][0]) in pq_list:
                pq_list.remove((temp_pq[pq][0], temp_pq[pq][1]))
                p_list.remove(temp_p[pq])
                q_list.remove(temp_q[pq])
        del temp_pq, temp_p, temp_q
        print(list_adj, list_num_adj)
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
                                                      windows_size=self.window_matrix[self.names[p]].loc[self.names[p]])
                            y = window_representation(self.series[self.names[q]],
                                                      windows_size=self.window_matrix[self.names[q]].loc[self.names[q]])
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
                                                     p_value=self.rank_using_p_value,
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


class PCTMI(CITMI):
    def __init__(self, series, sig_lev=0.05, lag_max=5, p_value=True, rank_using_p_value=False, verbose=True, num_processor=-1,
                 graphical_optimization=False):
        """
        PC for time series using TMI and CTMI
        :param series: d-time series (with possibility of different sampling rate)
        :param sig_lev: significance level. By default 0.05
        :param p_value: Use p_value for decision making. By default True
        :param verbose: Print results. By default: True
        :param num_processor: number of processors for parallelization. By default -1 (all)
        """
        CITMI.__init__(self, series, sig_lev, lag_max, p_value, rank_using_p_value, verbose, num_processor,
                       graphical_optimization)

    def _find_shortest_directed_paths_util(self, p, q, visited, path, all_path):
        """
        sub function of _find_shortest_directed_paths
        :param p: index of time series
        :param q: index of time series
        :param visited: list of visited nodes
        :param path: current path
        :param all_path: list of all discovered paths
        """
        # Mark the current node as visited and store in path
        visited[p] = True
        path.append(p)

        # If current vertex is same as destination, then print
        # current path[]
        if p == q:
            if len(path) > 2:
                all_path.append(path.copy()[1:-1])
                return path
        else:
            # If current vertex is not destination
            # Recur for all the vertices child of this vertex
            child_p = np.where(self.graph.edges[p, :] == 2)[0]
            for k in child_p:
                if not visited[k]:
                    self._find_shortest_directed_paths_util(k, q, visited, path, all_path)

        # Remove current vertex from path[] and mark it as unvisited
        path.pop()
        visited[p] = False

    def _find_shortest_directed_paths(self, p, q):
        """
        find shortest directed path between time series of index p and time series of index q
        :param p: index of time series
        :param q: index of time series
        :return: all directed paths from p to q
        """
        # Mark all the vertices as not visited
        visited = [False] * self.d

        # Create an array to store paths
        path = []
        all_path = []

        # Call the recursive helper function to print all paths
        self._find_shortest_directed_paths_util(p, q, visited, path, all_path)
        return all_path

    def _oc_pq(self, p, q):
        """
        estimate ctmi between two time series conditioned of their sep set + each non oriented connected neighbor
        :param p: index of time series
        :param q: index of time series
        :return: p, q, the most contributing node and the MI conditioned on sep set and the most contributing node
        """
        v_list = []
        k_list = [k for k in range(self.d) if (k != p) and (k != q) and (self.graph.sep[p, q, k] == 0)
                  and (((self.graph.edges[q, k] == 1) and (self.graph.edges[p, k] == 1)
                       and (self.graph.edges[k, q] == 1) and (self.graph.edges[k, p] == 1)) or(
                    (self.graph.edges[q, k] == 2) and (self.graph.edges[p, k] == 1) and
                    (self.graph.edges[k, q] == 1) and (self.graph.edges[k, p] == 1)) or(
                    (self.graph.edges[q, k] == 1) and (self.graph.edges[p, k] == 2) and
                    (self.graph.edges[k, q] == 1) and (self.graph.edges[k, p] == 1)))]
        if len(k_list) > 0:
            if self.adaptive_window:
                # x = window_representation(self.series[self.names[p]],
                #                           windows_size=self.window_matrix[self.names[p]].loc[self.names[p]])
                # y = window_representation(self.series[self.names[q]],
                #                           windows_size=self.window_matrix[self.names[q]].loc[self.names[q]])
                x = self.series[self.names[p]]
                y = self.series[self.names[q]]
            else:
                x = self.data_dict[self.names[p]]
                y = self.data_dict[self.names[q]]
            sep = np.where(self.graph.sep[p, q, :] == 1)[0]
            for k in k_list:
                if k not in sep:
                    # sep_k = sep.tolist() + [k]
                    # z = dict()
                    # for name_k in self.names[sep_k]:
                    #     z[name_k] = self.data_dict[name_k]
                    # cmi_pval, cmi_val = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
                    #            gamma_matrix=self.gamma_matrix, p_value=self.p_value, mission="ictmi")
                    if self.adaptive_window:
                        # select and drop NA
                        z = self.series[self.names[k]].dropna()
                    else:
                        z = self.data_dict[self.names[k]]

                    cmi_pval, cmi_val = i_ctmi(x, y, z, self.names[p], self.names[q], self.names[k], self.sampling_rate,
                               p_value=self.p_value)
                    if self.verbose:
                        print("p=" + str(p) + "; q=" + str(q) + "; r=" + str(
                            k) + "; s=" + str(
                            sep) + "; I(p,q|r,s)=" + "{: 0.5f}".format(
                            cmi_pval) + "; I(p,q|s)=")
                    v_list.append(cmi_pval)
        if v_list:
            if self.p_value:
                idx = int(np.argmax(v_list))
            else:
                idx = int(np.argmin(v_list))
            return p, q, v_list[idx], k_list[idx]

    def rank_oc_parallel(self):
        """
        rank unsheilded triples based on the estimation of ctmi
        :return: ranking of each unsheilded triple based ctmi
        """
        p_list = []
        q_list = []
        for p in range(self.d):
            for q in range(p+1, self.d):
                if (self.graph.edges[p, q] == 0) and (self.graph.edges[q, p] == 0):
                    p_list.append(p)
                    q_list.append(q)
        res = Parallel(n_jobs=self.num_processor)(delayed(self._oc_pq)(p, q) for p, q in zip(p_list, q_list))
        ranks = RankingList()
        for pq in range(len(res)):
            if res[pq] is not None:
                ranks.add(res[pq][0], res[pq][1], res[pq][2], res[pq][3])
        if self.p_value:
            ranks.sort(descending=False)
        else:
            ranks.sort(descending=True)
        return ranks

    def rule_origin_causality(self):
        """
        rule 0 (origin of causality) from PC
        """
        if self.verbose:
            print("######################################")
            print("Rule Origin of Causality")
            print("######################################")

        ranks = self.rank_oc_parallel()
        for p, q, k, cmi in zip(ranks.elem_p, ranks.elem_q, ranks.elem_r, ranks.val):
            if ((self.graph.edges[q, k] == 1) and (self.graph.edges[p, k] == 1) \
                    and (self.graph.edges[k, q] == 1) and (self.graph.edges[k, p] == 1)) or(
                    (self.graph.edges[q, k] == 2) and (self.graph.edges[p, k] == 1) and
                    (self.graph.edges[k, q] == 1) and (self.graph.edges[k, p] == 1)) or(
                    (self.graph.edges[q, k] == 1) and (self.graph.edges[p, k] == 2) and
                    (self.graph.edges[k, q] == 1) and (self.graph.edges[k, p] == 1)):
                sep = np.where(self.graph.sep[p, q, :] == 1)[0]
                print("sep = " + str(sep))
                # if len(sep) > 0:
                #     mi = self.cmi_array[p, q]
                # else:
                mi = self.mi_array[p, q]
                if k not in sep:
                    if self.verbose:
                        print("p=" + str(p) + "; q=" + str(q) + "; r=" + str(
                            k) + "; s=" + str(
                            sep) + "; I(p,q|r,s)=" + "{: 0.5f}".format(
                            cmi) + "; I(p,q|s)=" + "{: 0.5f}".format(mi), end=" ")
                    if self.p_value:
                        test = cmi <= mi
                    else:
                        test = mi < cmi
                    if test:
                        if self.verbose:
                            print("=> orient " + str(p) + " -> " + str(k) + " <- " + str(q))
                        self.graph.edges[p, k] = 2
                        self.graph.edges[q, k] = 2
                    else:
                        if self.verbose:
                            print()

    def rule_propagation_causality(self):
        """
        rule 1 from PC
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule Propagation of Causality")
            print("######################################")

        test_find_orientation = False

        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 0) and (self.graph.edges[j, i] == 0):
                    k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and
                              (((self.graph.edges[j, k] == 2) and (self.graph.edges[k, i] == 1) and
                                (self.graph.edges[i, k] == 1)) or ((self.graph.edges[i, k] == 2) and
                                                                   (self.graph.edges[k, j] == 1) and
                                                                   (self.graph.edges[j, k] == 1)))]
                    if len(k_list) > 0:
                        test_find_orientation = True
                        for k in k_list:
                            if self.graph.edges[i, k] == 2:
                                if self.verbose:
                                    print(str(i) + "->" + str(k) + "-" + str(j), end=" ")
                                    print("=> orient " + str(i) + "-> " + str(k) + " -> " + str(j))
                                self.graph.edges[k, j] = 2
                            else:
                                if self.verbose:
                                    print(str(j) + "->" + str(k) + "-" + str(i), end=" ")
                                    print("=> orient " + str(j) + "-> " + str(k) + " -> " + str(i))
                                self.graph.edges[k, i] = 2
        return test_find_orientation

    def rule_2(self):
        """
        rule 2 from PC
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule 3")
            print("######################################")
        test_find_orientation = False

        for i in range(self.graph.d):
            j_list = np.where(self.graph.edges[i, :] == 1)[0].tolist()
            if i in j_list:
                j_list.remove(i)
            for j in j_list:
                if self.graph.edges[j, i] == 1:
                    shortest_directed_path = self._find_shortest_directed_paths(i, j)
                    if len(shortest_directed_path) > 0:
                        self.graph.edges[i, j] = 2
                        test_find_orientation = True
                        if self.verbose:
                            print_path = '->'.join(map(str, shortest_directed_path[0]))
                            print(str(i)+"-"+str(j)+" and "+str(i) + "->" + print_path + "->" + str(j), end=" ")
                            print("=> orient " + str(i) + "->" + str(j))
        return test_find_orientation

    def rule_3(self):
        """
        rule 3 from PC
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule 4")
            print("######################################")

        test_find_orientation = False

        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 0) and (self.graph.edges[j, i] == 0):
                    colliders = [k for k in range(self.graph.d) if (k != i) and (k != j) and (
                            (self.graph.edges[j, k] == 2) and (self.graph.edges[i, k] == 2))]
                    k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and (
                            (self.graph.edges[j, k] == 1) and (self.graph.edges[i, k] == 1))
                              and (self.graph.edges[k, j] == 1) and (self.graph.edges[k, i] == 1)]
                    if len(colliders) > 0 and len(k_list) > 0:
                        for c in colliders:
                            for k in k_list:
                                if (self.graph.edges[c, k] == 1) and (self.graph.edges[k, c] == 1):
                                    test_find_orientation = True
                                    self.graph.edges[k, c] = 2
                                    if self.verbose:
                                        print(str(i) + "->" + str(c) + "<-" + str(j) + " and " + str(i) + "-" +
                                              str(k) + "-" + str(j) + " and " + str(k) + "-" + str(c),
                                              end=" ")
                                        print("=> orient " + str(k) + "->" + str(c))
        return test_find_orientation

    def rule_confounder(self):
        if self.verbose:
            print("######################################")
            print("Rule confounder")
            print("######################################")

        print(np.where(self.graph.sep == 1))
        p_list, q_list, r_list = np.where(self.graph.sep == 1)
        print(p_list)
        print(q_list)
        print(r_list)
        for p,q,r in zip(p_list, q_list, r_list):
            if p > q:
                print(p,q,r)
                if self.graph.edges[p, r] != 0:
                    if self.biggamma[p,q,r] > 0:
                        self.graph.edges[r, p] = 2
                    elif self.biggamma[p,q,r] < 0:
                        self.graph.edges[p, r] = 2
                if self.graph.edges[q, r] != 0:
                    if self.biggamma[q,p,r] > 0:
                        self.graph.edges[r, q] = 2
                    elif self.biggamma[q,p,r] < 0:
                        self.graph.edges[q, r] = 2

    def rule_commun_confounder_and_causal_chain(self):
        """
        new rules (rule commun confounder (4) and rule causal_chain (5) from paper)
        """
        if self.verbose:
            print("######################################")
            print("Rule commun confounder and causal chain")
            print("######################################")

        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 0) and (self.graph.edges[j, i] == 0):
                    # k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and (
                    #     ((self.graph.edges[j, k] == 1) and (self.graph.edges[k, j] == 1) and
                    #         (self.graph.edges[i, k] == 1) and (self.graph.edges[k, i] == 1)))]
                    k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and (
                        ((self.graph.edges[j, k] == 1) and (self.graph.edges[k, j] == 1) and
                            (self.graph.edges[i, k] == 1) and (self.graph.edges[k, i] == 1)) or (
                        ((self.graph.edges[j, k] == 2) and (self.graph.edges[k, j] == 1) and
                         (self.graph.edges[i, k] == 1) and (self.graph.edges[k, i] == 1))) or (
                        ((self.graph.edges[j, k] == 1) and (self.graph.edges[k, j] == 2) and
                         (self.graph.edges[i, k] == 1) and (self.graph.edges[k, i] == 1))) or (
                        ((self.graph.edges[j, k] == 1) and (self.graph.edges[k, j] == 1) and
                         (self.graph.edges[i, k] == 2) and (self.graph.edges[k, i] == 1))) or (
                        ((self.graph.edges[j, k] == 1) and (self.graph.edges[k, j] == 1) and
                         (self.graph.edges[i, k] == 1) and (self.graph.edges[k, i] == 2))))]

                    if len(k_list) > 0:
                        for k in k_list:
                            gki = self.gamma_matrix[self.names[i]].loc[self.names[k]]
                            gkj = self.gamma_matrix[self.names[j]].loc[self.names[k]]
                            i_is_not_effet = (sum(self.graph.edges[:, i] == 2) == 0)
                            j_is_not_effet = (sum(self.graph.edges[:, j] == 2) == 0)
                            k_is_not_effet = (sum(self.graph.edges[:, k] == 2) == 0)

                            #Lagged common cause
                            if (gki > 0) and (gkj > 0):
                                if i_is_not_effet and j_is_not_effet:
                                    if self.verbose:
                                        print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)>0",
                                              end=" ")
                                        print("=> orient " + str(i) + "<- " + str(k) + " -> " + str(j))
                                    self.graph.edges[k, i] = 2
                                    self.graph.edges[k, j] = 2
                            #Lagged instantaneous confounder
                            elif (gki > 0) and (gkj == 0):
                                if i_is_not_effet:
                                    if j_is_not_effet and k_is_not_effet:
                                        if self.verbose:
                                            print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)==0",
                                                  end=" ")
                                            print("=> orient " + str(i) + "<- " + str(k) + " - " + str(j))
                                        self.graph.edges[k, i] = 2
                                    elif j_is_not_effet:
                                        if self.verbose:
                                            print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)==0",
                                                  end=" ")
                                            print("=> orient " + str(i) + "<- " + str(k) + " -> " + str(j))
                                        self.graph.edges[k, i] = 2
                                        self.graph.edges[k, j] = 2
                                    elif k_is_not_effet:
                                        if self.verbose:
                                            print(
                                                str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)==0",
                                                end=" ")
                                            print("=> orient " + str(i) + "<- " + str(k) + " <- " + str(j))
                                        self.graph.edges[k, i] = 2
                                        self.graph.edges[j, k] = 2
                            elif (gki == 0) and (gkj > 0):
                                if j_is_not_effet:
                                    if i_is_not_effet and k_is_not_effet:
                                        if self.verbose:
                                            print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)==0 and gamma(k,j)>0",
                                                  end=" ")
                                            print("=> orient " + str(i) + "- " + str(k) + " -> " + str(j))
                                        self.graph.edges[k, j] = 2
                                    elif i_is_not_effet:
                                        if self.verbose:
                                            print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)==0",
                                                  end=" ")
                                            print("=> orient " + str(i) + "<- " + str(k) + " -> " + str(j))
                                        self.graph.edges[k, i] = 2
                                        self.graph.edges[k, j] = 2
                                    elif k_is_not_effet:
                                        if self.verbose:
                                            print(
                                                str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,j)>0 and gamma(k,i)==0",
                                                end=" ")
                                            print("=> orient " + str(i) + "-> " + str(k) + " -> " + str(j))
                                        self.graph.edges[k, j] = 2
                                        self.graph.edges[i, k] = 2
                            # lagged instanteneous causal chain
                            elif (gki >= 0) and (gkj < 0):
                                if j_is_not_effet and k_is_not_effet:
                                    if self.verbose:
                                        print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)>0",
                                              end=" ")
                                        print("=> orient " + str(i) + "<- " + str(k) + " <- " + str(j))
                                    self.graph.edges[k, i] = 2
                                    self.graph.edges[j, k] = 2
                            elif (gki < 0) and (gkj >= 0):
                                if i_is_not_effet and k_is_not_effet:
                                    if self.verbose:
                                        print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)>0",
                                              end=" ")
                                        print("=> orient " + str(i) + "-> " + str(k) + " -> " + str(j))
                                    self.graph.edges[i, k] = 2
                                    self.graph.edges[k, j] = 2

    def rule_mediator(self):
        """
        new rules (rule mediator (6) from paper)
        """
        if self.verbose:
            print("######################################")
            print("Rule mediator")
            print("######################################")

        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] != 0) and (self.graph.edges[j, i] != 0):
                    k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and (
                        ((self.graph.edges[j, k] == 1) and (self.graph.edges[k, j] == 1) and
                            (self.graph.edges[i, k] == 1) and (self.graph.edges[k, i] == 1)))]
                    if len(k_list) > 0:
                        for k in k_list:
                            gij = self.gamma_matrix[self.names[j]].loc[self.names[i]]
                            gik = self.gamma_matrix[self.names[k]].loc[self.names[i]]
                            gjk = self.gamma_matrix[self.names[k]].loc[self.names[j]]
                            # g_list = [(gij, gik, gjk), (gij, gik, -gjk), (-gij, gjk, gik), (-gij, gjk, -gik),
                            #           (-gik, -gjk, gij), (-gik, -gjk, -gij)]
                            # i->j->k, i->k->j, j->i->k, j->k->->i, k->i->j, k->j->i
                            g_list = [(gij, gjk, gik), (gik, -gjk, gij), (-gij, gik, gjk), (gjk, -gik, -gij),
                                      (-gik, gij, -gjk), (-gjk, -gij, -gik)]
                            g_list_common = [(gij, gjk, gik), (gik, -gjk, gij), (gjk, -gik, -gij)]

                            msk = [(x[0] > 0) and (x[1] > 0) and (x[2] >= 0) for x in g_list]
                            msk_common = [(x[0] == 0) and (x[1] > 0) and (x[2] > 0) for x in g_list_common]
                            if any(msk):
                                print(g_list)
                                print(msk)
                                s = int(np.argwhere(msk)[0])
                                # g1, g2, g3 = g_list[s]
                                if s == 0:
                                    if (sum(self.graph.edges[:, j] == 2) == 0) and (
                                            sum(self.graph.edges[:, k] == 2) == 0):
                                        if (self.graph.edges[j, i] == 1) and (self.graph.edges[k, i] == 1) \
                                                and (self.graph.edges[k, j] == 1):
                                            if self.verbose:
                                                print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "-> " + str(j) + " -> " + str(k) + " <- "
                                                      + str(i))
                                            self.graph.edges[i, j] = 2
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[j, k] = 2
                                elif s == 1:
                                    if (sum(self.graph.edges[:, j] == 2) == 0) and (
                                            sum(self.graph.edges[:, k] == 2) == 0):
                                        if (self.graph.edges[j, i] == 1) and (self.graph.edges[k, i] == 1) \
                                                and (self.graph.edges[k, j] == 1):
                                            if self.verbose:
                                                print(str(i) + "-" + str(k) + "-" + str(j) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "-> " + str(k) + "-> " + str(j) + " <- "
                                                      + str(i))
                                            self.graph.edges[i, j] = 2
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[k, j] = 2
                                elif s == 2:
                                    if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                            sum(self.graph.edges[:, k] == 2) == 0):
                                        if (self.graph.edges[i, j] == 1) and (self.graph.edges[k, i] == 1) \
                                                and (self.graph.edges[k, j] == 1):
                                            if self.verbose:
                                                print(str(j) + "-" + str(i) + "-" + str(k) + "-" + str(j), end=" ")
                                                print("=> orient " + str(j) + "-> " + str(i) + " -> " + str(k) + " <- "
                                                      + str(j))
                                            self.graph.edges[j, i] = 2
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[j, k] = 2
                                if s == 3:
                                    if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                            sum(self.graph.edges[:, k] == 2) == 0):
                                        if (self.graph.edges[i, j] == 1) and (self.graph.edges[i, k] == 1) \
                                                and (self.graph.edges[k, j] == 1):
                                            if self.verbose:
                                                print(str(j) + "-" + str(k) + "-" + str(i) + "-" + str(j), end=" ")
                                                print("=> orient " + str(j) + "-> " + str(k) + "-> " + str(i) + " <- "
                                                      + str(j))
                                            self.graph.edges[j, i] = 2
                                            self.graph.edges[k, i] = 2
                                            self.graph.edges[j, k] = 2
                                elif s == 4:
                                    if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                            sum(self.graph.edges[:, j] == 2) == 0):
                                        if (self.graph.edges[j, i] == 1) and (self.graph.edges[i, k] == 1) \
                                                and (self.graph.edges[j, k] == 1):
                                            if self.verbose:
                                                print(str(k) + "-" + str(i) + "-" + str(j) + "-" + str(k), end=" ")
                                                print("=> orient " + str(k) + "-> " + str(i) + "-> " + str(j) + " <- "
                                                      + str(k))
                                            self.graph.edges[i, j] = 2
                                            self.graph.edges[k, i] = 2
                                            self.graph.edges[k, j] = 2
                                elif s == 5:
                                    if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                            sum(self.graph.edges[:, j] == 2) == 0):
                                        if (self.graph.edges[i, j] == 1) and (self.graph.edges[i, k] == 1) \
                                                and (self.graph.edges[j, k] == 1):
                                            if self.verbose:
                                                print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                      + str(k))
                                            self.graph.edges[j, i] = 2
                                            self.graph.edges[k, i] = 2
                                            self.graph.edges[k, j] = 2
                            elif any(msk_common):
                                s = int(np.argwhere(msk_common)[0])
                                if s == 0:
                                    if (self.graph.edges[j, i] == 1) and (self.graph.edges[k, i] == 1) \
                                            and (self.graph.edges[k, j] == 1):
                                        if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                sum(self.graph.edges[:, j] == 2) == 0) and (
                                                sum(self.graph.edges[:, k] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "- " + str(j) + " -> " + str(k) + " <- "
                                                      + str(i))
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[j, k] = 2
                                        elif (sum(self.graph.edges[:, j] == 2) == 0) and (
                                                sum(self.graph.edges[:, k] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "-> " + str(j) + " -> " + str(k) + " <- "
                                                      + str(i))
                                            self.graph.edges[i, j] = 2
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[j, k] = 2
                                        elif (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                sum(self.graph.edges[:, k] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
                                                print("=> orient " + str(j) + "-> " + str(i) + " -> " + str(k) + " <- "
                                                      + str(j))
                                            self.graph.edges[j, i] = 2
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[j, k] = 2
                                elif s == 1:
                                    if (self.graph.edges[j, i] == 1) and (self.graph.edges[k, i] == 1) \
                                            and (self.graph.edges[k, j] == 1):
                                        if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                sum(self.graph.edges[:, j] == 2) == 0) and (
                                                sum(self.graph.edges[:, k] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(k) + "-" + str(j) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "- " + str(k) + " -> " + str(j) + " <- "
                                                      + str(i))
                                            self.graph.edges[i, j] = 2
                                            self.graph.edges[k, j] = 2
                                        elif (sum(self.graph.edges[:, k] == 2) == 0) and (
                                                sum(self.graph.edges[:, j] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "-> " + str(j) + " -> " + str(k) + " <- "
                                                      + str(i))
                                            self.graph.edges[i, j] = 2
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[k, j] = 2
                                        elif (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                sum(self.graph.edges[:, j] == 2) == 0):
                                            if self.verbose:
                                                print(str(k) + "-" + str(i) + "-" + str(j) + "-" + str(k), end=" ")
                                                print("=> orient " + str(k) + "-> " + str(i) + " -> " + str(j) + " <- "
                                                      + str(k))
                                            self.graph.edges[k, i] = 2
                                            self.graph.edges[k, j] = 2
                                            self.graph.edges[i, j] = 2
                                elif s == 2:
                                    if (self.graph.edges[j, i] == 1) and (self.graph.edges[k, i] == 1) \
                                            and (self.graph.edges[k, j] == 1):
                                        if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                sum(self.graph.edges[:, j] == 2) == 0) and (
                                                sum(self.graph.edges[:, k] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(k) + "-" + str(j) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "- " + str(k) + " -> " + str(j) + " <- "
                                                      + str(i))
                                            self.graph.edges[j, i] = 2
                                            self.graph.edges[k, i] = 2
                                        elif (sum(self.graph.edges[:, k] == 2) == 0) and (
                                                sum(self.graph.edges[:, i] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "-> " + str(j) + " -> " + str(k) + " <- "
                                                      + str(i))
                                            self.graph.edges[j, k] = 2
                                            self.graph.edges[j, i] = 2
                                            self.graph.edges[k, i] = 2
                                        elif (sum(self.graph.edges[:, j] == 2) == 0) and (
                                                sum(self.graph.edges[:, i] == 2) == 0):
                                            if self.verbose:
                                                print(str(k) + "-" + str(i) + "-" + str(j) + "-" + str(k), end=" ")
                                                print("=> orient " + str(k) + "-> " + str(i) + " -> " + str(j) + " <- "
                                                      + str(k))
                                            self.graph.edges[k, j] = 2
                                            self.graph.edges[k, i] = 2
                                            self.graph.edges[j, i] = 2
                            else:
                                if (self.graph.edges[i, j] == 1) and (self.graph.edges[i, k] == 1) \
                                        and (self.graph.edges[k, j] == 1):
                                    if (gij!=0) and (gik==0) and (gjk==0):
                                        if gij>0:
                                            if (sum(self.graph.edges[:, j] == 2) == 0) and (
                                                    sum(self.graph.edges[:, k] == 2) == 0):
                                                if self.verbose:
                                                    print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                    print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                          + str(k))
                                                self.graph.edges[i, j] = 2
                                                self.graph.edges[i, k] = 2
                                                self.graph.edges[j, k] = 2
                                        elif gij<0:
                                            if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                    sum(self.graph.edges[:, k] == 2) == 0):
                                                if self.verbose:
                                                    print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                    print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                          + str(k))
                                                self.graph.edges[j, i] = 2
                                                self.graph.edges[i, k] = 2
                                                self.graph.edges[j, k] = 2
                                    elif (gij==0) and (gik!=0) and (gjk==0):
                                        if gik>0:
                                            if (sum(self.graph.edges[:, k] == 2) == 0) and (
                                                    sum(self.graph.edges[:, j] == 2) == 0):
                                                if self.verbose:
                                                    print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                    print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                          + str(k))
                                                self.graph.edges[i, k] = 2
                                                self.graph.edges[i, j] = 2
                                                self.graph.edges[k, j] = 2
                                        if gik<0:
                                            if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                    sum(self.graph.edges[:, j] == 2) == 0):
                                                if self.verbose:
                                                    print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                    print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                          + str(k))
                                                self.graph.edges[k, i] = 2
                                                self.graph.edges[i, j] = 2
                                                self.graph.edges[k, j] = 2
                                    elif (gij == 0) and (gik == 0) and (gjk != 0):
                                        if gjk>0:
                                            if (sum(self.graph.edges[:, k] == 2) == 0) and (
                                                    sum(self.graph.edges[:, i] == 2) == 0):
                                                if self.verbose:
                                                    print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                    print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                          + str(k))
                                                self.graph.edges[j, k] = 2
                                                self.graph.edges[j, i] = 2
                                                self.graph.edges[k, i] = 2
                                        if gjk<0:
                                            if (sum(self.graph.edges[:, j] == 2) == 0) and (
                                                    sum(self.graph.edges[:, i] == 2) == 0):
                                                if self.verbose:
                                                    print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                    print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                          + str(k))
                                                self.graph.edges[k, j] = 2
                                                self.graph.edges[j, i] = 2
                                                self.graph.edges[k, i] = 2

    def rule_proba_raising_principle(self):
        """
        new rules (rule prob raising principle from paper)
        """
        if self.verbose:
            print("######################################")
            print("Rule prob raising principle")
            print("######################################")
        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 1) and (self.graph.edges[j, i] == 1):

                    adjacent_i_is_1 = (sum(np.delete(self.graph.edges[:, i], i) != 0) == 1)
                    adjacent_j_is_1 = (sum(np.delete(self.graph.edges[:, j], j) != 0) == 1)
                    if adjacent_i_is_1 and adjacent_j_is_1:
                        gij = self.gamma_matrix[self.names[j]].loc[self.names[i]]
                        if gij > 0:
                            if self.verbose:
                                print(str(i) + "-" + str(j) + "g(i,j)>0", end=" ")
                                print("=> orient " + str(i) + "-> " + str(j))
                            self.graph.edges[i, j] = 2
                        elif gij < 0:
                            if self.verbose:
                                print(str(i) + "-" + str(j) + "g(i,j)<0", end=" ")
                                print("=> orient " + str(i) + "<- " + str(j))
                            self.graph.edges[j, i] = 2

    def gamma_zero_different_sampling_rate(self, i ,j):
        valid_idx_i = self.series[self.names[i]].first_valid_index()
        valid_idx_j = self.series[self.names[j]].first_valid_index()
        print(valid_idx_i, valid_idx_j)
        if valid_idx_i < valid_idx_j:
            self.graph.edges[i, j] = 2
            if self.verbose:
                print(str(i) + "-" + str(j) + "g(i,j)=0", end=" ")
                print("=> orient " + str(i) + "-> " + str(j))
        elif valid_idx_i < valid_idx_j:
            self.graph.edges[j, i] = 2
            if self.verbose:
                print(str(i) + "-" + str(j) + "g(i,j)=0", end=" ")
                print("=> orient " + str(i) + "<- " + str(j))

    def rule_entropy_reduction_gamma(self):
        """
        new rules (rule prob raising principle from paper)
        """
        if self.verbose:
            print("######################################")
            print("Rule entropy reduction principle Gamma")
            print("######################################")
        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 1) and (self.graph.edges[j, i] == 1):
                    gij = self.gamma_matrix[self.names[j]].loc[self.names[i]]
                    if gij > 0:
                        # count_possible_unsheilded_collider = sum(self.graph.edges[self.graph.edges[:, j] == 2, i] == 0)
                        # if count_possible_unsheilded_collider > 0:
                        if self.verbose:
                            print(str(i) + "-" + str(j) + "g(i,j)>0", end=" ")
                            print("=> orient " + str(i) + "-> " + str(j))
                        self.graph.edges[i, j] = 2
                    elif gij < 0:
                        # count_possible_unsheilded_collider = sum(self.graph.edges[self.graph.edges[:, i] == 2, j] == 0)
                        # if count_possible_unsheilded_collider > 0:
                        if self.verbose:
                            print(str(i) + "-" + str(j) + "g(i,j)<0", end=" ")
                            print("=> orient " + str(i) + "<- " + str(j))
                        self.graph.edges[j, i] = 2
                    else:
                        print("gamma")
                        print(gij)
                        print("sampling rate")
                        print(self.sampling_rate[self.names[i]], self.sampling_rate[self.names[j]])
                        if (self.sampling_rate[self.names[i]] != 1) or (self.sampling_rate[self.names[j]] != 1):
                            self.gamma_zero_different_sampling_rate(i, j)


    def rule_entropy_reduction_lambda(self):
        """
        new rules (rule prob raising principle from paper)
        """
        if self.verbose:
            print("######################################")
            print("Rule entropy reduction principle Lambda")
            print("######################################")
        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 1) and (self.graph.edges[j, i] == 1):
                    gij = self.gamma_matrix[self.names[j]].loc[self.names[i]]
                    if gij == 0:
                        test_parents = True
                        window_min = min(self.window_matrix[self.names[j]].loc[self.names[i]], self.window_matrix[self.names[i]].loc[
                            self.names[j]])
                        for k in range(self.graph.d):
                            if (self.graph.edges[k, j] == 2) and (self.graph.edges[k, i] == 2):
                                if (self.gamma_matrix[self.names[i]].loc[self.names[k]]>= window_min) or (
                                        self.gamma_matrix[self.names[j]].loc[self.names[k]]>= window_min):
                                    test_parents = False
                        if test_parents:
                            if self.window_matrix[self.names[j]].loc[self.names[i]] < self.window_matrix[self.names[i]].loc[self.names[j]]:
                                if self.verbose:
                                    print(str(i) + "-" + str(j) + "g(i,j)=0 and w(i,j)< w(j,i)", end=" ")
                                    print("=> orient " + str(i) + "-> " + str(j))
                                self.graph.edges[i, j] = 2
                            else:
                                if self.verbose:
                                    print(str(i) + "-" + str(j) + "g(i,j)=0 and w(i,j)> w(j,i)", end=" ")
                                    print("=> orient " + str(j) + "-> " + str(i))
                                self.graph.edges[j, i] = 2

    # def rule_mediator(self):
    #     """
    #     new rules (rule mediator (6) from paper)
    #     """
    #     if self.verbose:
    #         print("######################################")
    #         print("Rule mediator")
    #         print("######################################")
    #
    #     for i in range(self.graph.d):
    #         for j in range(i + 1, self.graph.d):
    #             if (self.graph.edges[i, j] != 0) and (self.graph.edges[j, i] != 0):
    #                 k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and (
    #                     ((self.graph.edges[j, k] == 1) and (self.graph.edges[k, j] == 1) and
    #                         (self.graph.edges[i, k] == 1) and (self.graph.edges[k, i] == 1)))]
    #                 if len(k_list) > 0:
    #                     for k in k_list:
    #                         gij = self.gamma_matrix[self.names[j]].loc[self.names[i]]
    #                         gik = self.gamma_matrix[self.names[k]].loc[self.names[i]]
    #                         gjk = self.gamma_matrix[self.names[k]].loc[self.names[j]]
    #                         g_list = [(gij, gik), (-gij, gjk), (-gik, -gjk)]
    #                         msk = [(x[0] > 0) and (x[1] > 0) for x in g_list]
    #                         if any(msk):
    #                             print(g_list)
    #                             print(msk)
    #                             s = int(np.argwhere(msk)[0])
    #                             g1, g2 = g_list[s]
    #                             if s == 0:
    #                                 #todo (done) check if j and k are the effect of any other node
    #                                 if (sum(self.graph.edges[:, j] == 2) == 0) and (
    #                                         sum(self.graph.edges[:, k] == 2) == 0):
    #                                     if g1 < g2:
    #                                         if (self.graph.edges[j, i] == 1) and (self.graph.edges[k, i] == 1) \
    #                                                 and (self.graph.edges[k, j] == 1):
    #                                             if self.verbose:
    #                                                 print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
    #                                                 print("=> orient " + str(i) + "-> " + str(j) + " -> " + str(k) + " <- "
    #                                                       + str(i))
    #                                             self.graph.edges[i, j] = 2
    #                                             self.graph.edges[i, k] = 2
    #                                             self.graph.edges[j, k] = 2
    #                                     elif g1 > g2:
    #                                         if (self.graph.edges[j, i] == 1) and (self.graph.edges[k, i] == 1) \
    #                                                 and (self.graph.edges[j, k] == 1):
    #                                             if self.verbose:
    #                                                 print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
    #                                                 print("=> orient " + str(i) + "-> " + str(j) + " <- " + str(k) + " <- "
    #                                                       + str(i))
    #                                             self.graph.edges[i, j] = 2
    #                                             self.graph.edges[i, k] = 2
    #                                             self.graph.edges[k, j] = 2
    #                             elif s == 1:
    #                                 if (sum(self.graph.edges[:, i] == 2) == 0) and (
    #                                         sum(self.graph.edges[:, k] == 2) == 0):
    #                                     if g1 < g2:
    #                                         if (self.graph.edges[i, j] == 1) and (self.graph.edges[k, i] == 1) \
    #                                                 and (self.graph.edges[k, j] == 1):
    #                                             if self.verbose:
    #                                                 print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
    #                                                 print("=> orient " + str(i) + "<- " + str(j) + " -> " + str(k) + " <- "
    #                                                       + str(i))
    #                                             self.graph.edges[j, i] = 2
    #                                             self.graph.edges[i, k] = 2
    #                                             self.graph.edges[j, k] = 2
    #                                     elif g1 > g2:
    #                                         if (self.graph.edges[i, j] == 1) and (self.graph.edges[i, k] == 1) \
    #                                                 and (self.graph.edges[k, j] == 1):
    #                                             if self.verbose:
    #                                                 print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
    #                                                 print("=> orient " + str(i) + "<- " + str(j) + " <- " + str(k) + " -> "
    #                                                       + str(i))
    #                                             self.graph.edges[j, i] = 2
    #                                             self.graph.edges[k, i] = 2
    #                                             self.graph.edges[j, k] = 2
    #                             elif s == 2:
    #                                 if (sum(self.graph.edges[:, i] == 2) == 0) and (
    #                                         sum(self.graph.edges[:, j] == 2) == 0):
    #                                     if g1 < g2:
    #                                         if (self.graph.edges[j, i] == 1) and (self.graph.edges[i, k] == 1) \
    #                                                 and (self.graph.edges[j, k] == 1):
    #                                             if self.verbose:
    #                                                 print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
    #                                                 print("=> orient " + str(i) + "-> " + str(j) + " <- " + str(k) + " -> "
    #                                                       + str(i))
    #                                             self.graph.edges[i, j] = 2
    #                                             self.graph.edges[k, i] = 2
    #                                             self.graph.edges[k, j] = 2
    #                                     elif g1 > g2:
    #                                         if (self.graph.edges[i, j] == 1) and (self.graph.edges[i, k] == 1) \
    #                                                 and (self.graph.edges[j, k] == 1):
    #                                             if self.verbose:
    #                                                 print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
    #                                                 print("=> orient " + str(i) + "<- " + str(j) + " <- " + str(k) + " -> "
    #                                                       + str(i))
    #                                             self.graph.edges[j, i] = 2
    #                                             self.graph.edges[k, i] = 2
    #                                             self.graph.edges[k, j] = 2

    def check_self_loops(self):
        if self.verbose:
            print("#######################################")
            print("########### Check past self causes ###########")
            print("#######################################")

        for p in range(self.d):
            r_list = [r for r in range(self.graph.d) if (r != p) and ((self.graph.edges[r, p] == 2))]

            # x = window_representation(self.series[self.names[p]], windows_size=5, overlap=True)
            # x_past = x[x.columns[:-1]]
            # x = x[x.columns[-1]]

            z = dict()
            for r in r_list:
                z[self.names[r]] = self.series[self.names[r]]
            cmi_pval, _ = self_ctmi(self.series[self.names[p]], z, self.names[p], self.gamma_matrix, self.sampling_rate, instantaneous_dict=self.instantaneous_dict, p_value=self.p_value)

            if cmi_pval > self.sig_lev:
                self.graph.edges[p, p] = 0
                if self.verbose:
                    print(str(p) + "is independent from its past conditionned on "+str(r_list) + ": pvalue = "+str(cmi_pval))
                    print("=> remove link between " + str(p) + " and " + str(p))

    # def check_past(self):
    #     if self.verbose:
    #         print("#######################################")
    #         print("########### Check past self causes ###########")
    #         print("#######################################")
    #
    #     data_dict_past, gamma_matrix_past = self._include_past_check_past()
    #
    #     for p in range(self.d):
    #         r_list = [r for r in range(self.graph.d) if (r != p) and ((self.graph.edges[r, p] == 2))]
    #
    #         x = data_dict_past[self.names[p]]
    #         x_past = data_dict_past[self.names[p]]
    #
    #         for rs in r_list:
    #             z = dict()
    #             for r in rs:
    #                 z[self.names[r]] = data_dict_past[self.names[r]]
    #             cmi_pval, cmi_val = ctmi(x, x_past, z, self.names[p], self.names[x_past], self.sampling_rate,
    #                                      gamma_matrix=gamma_matrix_past, p_value=self.rank_using_p_value,
    #                                      instantaneous_dict=self.instantaneous_dict)
    #
    #             if cmi_pval>=self.alpha:
    #                 self.graph.edges[p, p] = 0
    #     self._exclude_past()
    #
    # def _include_past_check_past(self):
    #     counter = 0
    #     lags_past = []
    #     data_dict_past = dict()
    #     for col in range(self.d):
    #         lags_past.append(1)
    #         data_dict_past[self.names[col]] = window_representation(self.series[self.names[col]],
    #                                                                     windows_size=self.lags[col])
    #     gamma_matrix_past = align_matrix(data_dict_past, self.names, self.sampling_rate)
    #
    #     for i in range(self.d):
    #         self.instantaneous_dict[self.names[i] + "_past"] = False
    #         if gamma_matrix_past[self.names[i]].loc[self.names[i]] > 0:
    #             self.names = self.names.insert(self.d + counter, self.names[i] + "_past")
    #
    #             counter = counter + 1
    #             gamma_matrix_past[self.names[i] + "_past"] = 0
    #             gamma_matrix_past = gamma_matrix_past.append(pd.Series([0] * (self.d + counter),
    #                                                                    name=self.names[i] + "_past",
    #                                                                    index=gamma_matrix_past.columns))
    #
    #             gamma_matrix_past[self.names[i]].loc[self.names[i] + "_past"] = \
    #                 gamma_matrix_past[self.names[i]].loc[self.names[i]]
    #             gamma_matrix_past[self.names[i] + "_past"].loc[self.names[i]] = - \
    #                 gamma_matrix_past[self.names[i]].loc[self.names[i]]
    #
    #             self.graph.edges[self.d + i, i] = 2
    #             self.graph.edges[i, self.d + i] = 1
    #
    #             data_dict_past[self.names[i] + "_past"] = data_dict_past[self.names[i]]
    #             self.sampling_rate[self.names[i] + "_past"] = self.sampling_rate[self.names[i]]
    #
    #     return data_dict_past, gamma_matrix_past

    def fit(self):
        """
        run PCTMI
        :return: graph (CPDAG)
        """
        if self.verbose:
            now = datetime.now()
            print("#######################################")
            print("########### Starting PCTMI ###########")
            print("########### " + now.strftime("%H:%M:%S" + " ###########"))
            print("#######################################")

        # initialize skeleton
        self.skeleton_initialize()

        # get separation sets
        self.find_sep_set()


        self.rule_entropy_reduction_gamma()
        self.rule_entropy_reduction_lambda()
        # self.rule_confounder()

        # orientation
        self.rule_origin_causality()

        test_rp = True
        test_r2 = True
        test_r3 = True
        while test_rp or test_r2 or test_r3:
            test_rp = self.rule_propagation_causality()
            test_r2 = self.rule_2()
            test_r3 = self.rule_3()

        # self.rule_commun_confounder_and_causal_chain()
        # self.rule_mediator()
        # self.rule_proba_raising_principle()

        # check self causes
        self.check_self_loops()

        if self.verbose:
            print("######################################")
            print("Final Results (PCTMI)")
            print("######################################")
            print("Summary Graph:")
            print(self.graph.edges)
        return self.graph.edges

    def rule_gap_orientation(self):
        """
        gamma heuristic rule from paper
        """
        if self.verbose:
            print("######################################")
            print("Rule gap orientation")
            print("######################################")

        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 1) and (self.graph.edges[j, i] == 1):
                    if self.gamma_matrix[self.names[j]].loc[self.names[i]] > 0:
                        if self.verbose:
                            print(str(i) + "-" + str(j) + "g(i,j)>0", end=" ")
                            print("=> orient " + str(i) + "-> " + str(j))
                        self.graph.edges[i, j] = 2
                    if self.gamma_matrix[self.names[j]].loc[self.names[i]] < 0:
                        if self.verbose:
                            print(str(i) + "-" + str(j) + "g(i,j)<0", end=" ")
                            print("=> orient " + str(i) + "<- " + str(j))
                        self.graph.edges[j, i] = 2

    def fit_gap_orientation(self):
        """
        run PCTMI-gamma (requirements: run PCTMI)
        :return: graph (CPDAG)
        """
        self.rule_gap_orientation()
        return self.graph.edges


class TPCTMI(PCTMI):
    def __init__(self, series, sig_lev=0.05, lag_max=5, p_value=True, rank_using_p_value=False, verbose=True, num_processor=-1,
                 graphical_optimization=False):
        """
        PC for time series using TMI and CTMI
        :param series: d-time series (with possibility of different sampling rate)
        :param sig_lev: significance level. By default 0.05
        :param p_value: Use p_value for decision making. By default True
        :param verbose: Print results. By default: True
        :param num_processor: number of processors for parallelization. By default -1 (all)
        """
        PCTMI.__init__(self, series, sig_lev, lag_max, p_value, rank_using_p_value, verbose, num_processor,
                       graphical_optimization)
        # self.ts_window_size = 0
        self.ts_data_dict = dict()
        self.ts_names = []
        self.ts_names_dict = dict()
        self.ts_data_df = pd.DataFrame()

        self.tgraph_dict = dict()
        for name in self.names:
            self.tgraph_dict[name] = []

    def check_cycles(self):
        if self.verbose:
            print("######################################")
            print("Check Cycles")
            print("######################################")
        print(self.graph.edges)
        for p in range(self.graph.d):
            for q in range(self.graph.d):
                if (self.graph.edges[p, q] == 2) and (self.graph.edges[q, p] == 1):
                    temp_matrix = self.gamma_matrix.copy()
                    temp_matrix[self.names[p] + "_future"] = None
                    temp_matrix = temp_matrix.append(pd.Series([None]*(self.graph.d+1), name=self.names[p] + "_future",
                                                               index=temp_matrix.columns))
                    for i in range(self.graph.d):
                        if i != p:
                            x = self.data_dict[self.names[i]]
                            y = self.data_dict[self.names[p]]
                            # todo lag max
                            g = align_pair(x, y, (self.sampling_rate[self.names[i]], self.sampling_rate[self.names[p]]),
                                           max_gamma=5, set_numbers="N")
                            temp_matrix[self.names[p] + "_future"].loc[self.names[i]] = g
                            temp_matrix[self.names[i]].loc[self.names[p] + "_future"] = -g
                        else:
                            temp_matrix[self.names[p] + "_future"].loc[self.names[p] + "_future"] = self.gamma_matrix[self.names[p]].loc[self.names[p]]
                            temp_matrix[self.names[p]].loc[self.names[p] + "_future"] = self.gamma_matrix[self.names[p]].loc[self.names[p]]
                            temp_matrix[self.names[p] + "_future"].loc[self.names[p]] = self.gamma_matrix[self.names[p]].loc[self.names[p]]

                    par_q = np.where(self.graph.edges[:, q] == 2)[0]
                    x = self.data_dict[self.names[q]]
                    y = self.data_dict[self.names[p]]
                    z = dict()
                    for name_z in self.names[par_q]:
                        z[name_z] = self.data_dict[name_z]

                    temp_sampling_rate = self.sampling_rate.copy()
                    temp_sampling_rate[self.names[p] + "_future"] = self.sampling_rate[self.names[p]]
                    cmi, _ = ctmi(x, y, z, self.names[q], self.names[p] + "_future", temp_sampling_rate,
                                  gamma_matrix=temp_matrix, p_value=self.p_value)
                    if self.verbose:
                        print("p=" + str(p) + "; q=" + str(q) + "; r=" + str(par_q) + "; I(p,q|r)=" + "{: 0.5f}".format(
                            cmi), end=" ")
                    if cmi > self.alpha:
                        if self.verbose:
                            print(str(p) + "->" + str(q), end=" ")
                            print("=> orient " + str(p) + "->" + str(q) + " and " + str(q) + "->" + str(p))
                        self.graph.edges[q, p] = 2

    def _ts_init(self):
        lag_gamma = []
        for i in range(self.d):
            for j in range(self.d):
                if (self.graph.edges[i, j] == 2) and (self.graph.edges[j, i] == 1):
                    lag_gamma.append(abs(self.gamma_matrix[self.names[j]].loc[self.names[i]]) + self.lags[j])
                elif (self.graph.edges[i, j] == 1) and (self.graph.edges[j, i] == 2):
                    lag_gamma.append(abs(self.gamma_matrix[self.names[i]].loc[self.names[j]]) + self.lags[i])
                elif (self.graph.edges[i, j] == 1) and (self.graph.edges[j, i] == 1):
                    lag_gamma.append(abs(self.gamma_matrix[self.names[j]].loc[self.names[i]]) + self.lags[j])
                    lag_gamma.append(abs(self.gamma_matrix[self.names[i]].loc[self.names[j]]) + self.lags[i])
        self.ts_window_size = max(lag_gamma)

        for col in range(self.graph.d):
            self.ts_data_dict[self.names[col]] = window_representation(self.series[self.names[col]],
                                                                       windows_size=self.ts_window_size)

        for col in range(self.d):
            self.ts_names_dict[self.names[col]] = self.ts_data_dict[self.names[col]].columns
            self.ts_names.extend(self.ts_data_dict[self.names[col]].columns.values)
            self.ts_data_df[self.ts_data_dict[self.names[col]].columns] = self.ts_data_dict[self.names[col]]
        self.ts_names_dict_inv = {v: k for k, v_list in self.ts_names_dict.items() for v in v_list}

        temporal_to_unit_idx = []
        ts_sampling_rate_temp = []
        j = 0
        for i in range(len(self.sampling_rate)):
            ts_sampling_rate_temp = ts_sampling_rate_temp + \
                                    [self.sampling_rate[self.names[i]]] * self.ts_window_size
            temporal_to_unit_idx = temporal_to_unit_idx + [j] * self.ts_window_size
            j = j + 1
        ts_sampling_rate_temp = np.array(ts_sampling_rate_temp)
        self.ts_sampling_rate = dict()
        for i in range(len(self.ts_names)):
            self.ts_sampling_rate[self.ts_names[i]] = ts_sampling_rate_temp[i]

        temp = []
        j = 0
        for i in range(self.d):
            for lag in range(self.ts_window_size):
                if lag == 0:
                    temp.append(lag)
                else:
                    temp.append(lag + self.ts_sampling_rate[self.ts_names[j]] - 1)
                j = j + 1
        self.time_table = pd.DataFrame([temp], columns=self.ts_names)

        # self.ts_mi_array = np.ones([self.graph.d, self.graph.d, self.ts_window_size])
        self.ts_mi_array = np.ones([self.ts_window_size * self.d, self.ts_window_size * self.d])

        # self.ts_cmi_array = np.ones([self.graph.d, self.graph.d, self.ts_window_size])
        self.ts_cmi_array = np.ones([self.ts_window_size * self.d, self.ts_window_size * self.d])

        self.tgraph = Graph(self.d * self.ts_window_size)

        if self.verbose:
            print("ts_names: " + str(self.ts_names))
            print("ts_sampling_rate: " + str(self.ts_sampling_rate))
            print("time_table: " + str(self.time_table))
            print("temporal window size:" + str(self.ts_window_size))

    def summary_to_temporal_array(self):
        """
        transfer knowledge from summary graph to temporal graph
        """

        self._ts_init()
        time_table = self.time_table.values.reshape(-1)
        counter = 0

        for i in range(len(self.names)):
            # lag of each time series i
            lag = len(self.ts_names_dict[self.names[i]])
            # use temporal priority property to add arrows from past to future in the same time series
            if self.graph.edges[i, i] == 1:
                a = self.temporal_priority(self.tgraph.edges[counter:counter + lag, counter:counter + lag])
                self.tgraph.edges[counter:counter + lag, counter:counter + lag] = a
            # iterate on each time step in time series i
            for ti in range(counter, counter + lag):
                counter_j = counter + lag
                for j in range(i + 1, len(self.names)):
                    sep = np.where(self.graph.sep[i, j, :] == 1)[0]
                    # lag of each time series j
                    lag_j = len(self.ts_names_dict[self.names[j]])
                    # iterate on each time step in time series j
                    for tj in range(counter_j, counter_j + lag_j):
                        # in case of independency between two time series, all sub variables will be independent
                        if self.graph.edges[i, j] == 0:
                            self.tgraph.edges[ti, tj] = 0
                            self.tgraph.edges[tj, ti] = 0

                            for s in sep:
                                counter_s = self.ts_window_size*s
                                # todo: take into account only gamma
                                g_is = self.gamma_matrix[self.names[i]].loc[self.names[s]]
                                g_js = self.gamma_matrix[self.names[j]].loc[self.names[s]]
                                for ts in range(counter_s, counter_s + self.ts_window_size):
                                    if (ts < ti) or (ts < tj):
                                        self.tgraph.sep[ti, tj, ts] = 1

                        elif self.graph.edges[i, j] == 2:
                            if time_table[tj]-time_table[ti] >= self.gamma_matrix[self.names[j]].loc[self.names[i]]:
                                self.tgraph.edges[ti, tj] = 2
                                self.tgraph.edges[tj, ti] = 1
                            else:
                                self.tgraph.edges[ti, tj] = 0
                                self.tgraph.edges[tj, ti] = 0
                        elif self.graph.edges[j, i] == 2:
                            if time_table[ti]-time_table[tj] >= self.gamma_matrix[self.names[i]].loc[self.names[j]]:
                                self.tgraph.edges[tj, ti] = 2
                                self.tgraph.edges[ti, tj] = 1
                            else:
                                self.tgraph.edges[ti, tj] = 0
                                self.tgraph.edges[tj, ti] = 0
                        else:
                            # in case of dependency between two time series, transform arrows to circles or dots
                            if time_table[ti] < time_table[tj]:
                                if time_table[tj] - time_table[ti] >= self.gamma_matrix[self.names[j]].loc[self.names[i]]:
                                    self.tgraph.edges[ti, tj] = 2
                                    self.tgraph.edges[tj, ti] = 1
                                else:
                                    self.tgraph.edges[ti, tj] = 0
                                    self.tgraph.edges[tj, ti] = 0
                            elif time_table[ti] > time_table[tj]:
                                if time_table[ti] - time_table[tj] >= self.gamma_matrix[self.names[i]].loc[self.names[j]]:
                                    self.tgraph.edges[tj, ti] = 2
                                    self.tgraph.edges[ti, tj] = 1
                                else:
                                    self.tgraph.edges[ti, tj] = 0
                                    self.tgraph.edges[tj, ti] = 0
                            else:
                                self.tgraph.edges[ti, tj] = 1
                                self.tgraph.edges[tj, ti] = 1
                    counter_j = counter_j + lag_j
            counter = counter + lag

    @staticmethod
    def temporal_priority(a):
        """
        :param a: square matrix that represent relations between sub variables in the same time series
        :param n_markov: if false consider one markov in the same time series
        :return: square matrix that takes into account temporal priority
        """
        for i in range(a.shape[0]-1):
            a[i, i] = 0
            for j in range(i+1, a.shape[0]):
                    a[i, j] = 2
                    a[j, i] = 1
        return a

    def temporal_to_temporal_dict(self):
        if self.verbose:
            print("######################################")
            print("Summary to temporal graph")
            print("######################################")

        for i in range(self.graph.d):
            for j in range(i, self.graph.d):
                if i == j:
                    # self.tgraph_dict[self.names[j]].append((self.names[i], -self.gamma_matrix[self.names[j]].loc[self.names[i]]))
                    for lag in range(self.gamma_matrix[self.names[j]].loc[self.names[i]],
                                     self.gamma_matrix[self.names[j]].loc[self.names[i]] + self.lags[j]):
                        self.tgraph_dict[self.names[j]].append((self.names[i], -lag))
                else:
                    if self.graph.edges[i, j] == 2:
                        for lag in range(self.gamma_matrix[self.names[j]].loc[self.names[i]], self.gamma_matrix[self.names[j]].loc[self.names[i]] + self.lags[j]):
                            self.tgraph_dict[self.names[j]].append((self.names[i], -lag))
                            for l in range(self.ts_window_size - lag):
                                self.tgraph[self.d*l + i, self.d*(l + lag) + j] = 2
                    elif self.graph.edges[j, i] == 2:
                        for lag in range(self.gamma_matrix[self.names[i]].loc[self.names[j]], self.gamma_matrix[self.names[i]].loc[self.names[j]] + self.lags[i]):
                            self.tgraph_dict[self.names[i]].append((self.names[j], -lag))
                            for l in range(self.ts_window_size - lag):
                                self.tgraph[self.d*l + j, self.d*(l + lag) + i] = 2
                    elif (self.graph.edges[i, j] == 1) and (self.graph.edges[j, i] == 1):
                        for lag in range(self.ts_window_size):
                            self.tgraph_dict[self.names[j]].append((self.names[i], -lag))
                            self.tgraph_dict[self.names[i]].append((self.names[j], -lag))
                        # for lag in range(self.gamma_matrix[j].loc[i], self.lags[j]):
                        #     self.tgraph_dict[self.names[i]].append((self.names[i], lag))
                        # for lag in range(self.gamma_matrix[i].loc[j], self.lags[i]):
                        #     self.tgraph_dict[self.names[j]].append((self.names[j], lag))

    #todo

    def _ts_cmi_sep_set_pq(self, p, q, set_size):
        """
        estimate ctmi between two time series conditioned on each set of neighbors with cardinality equal to set_size
        :param p: time series with index p
        :param q: time series with index q
        :param set_size: cardinality of the set of neighbors
        :return: p, q, list if estimated value of ctmi(p,q,r_set), and list of all r_sets
        """
        zeros_matrix = pd.DataFrame(np.zeros([self.tgraph.d, self.tgraph.d]), columns=self.ts_names,
                                    index=self.ts_names)
        v_list = []
        r_list = [r for r in range(self.tgraph.d) if (r != p) and (r != q) and ((
                self.tgraph.edges[p, r] == 2) or (self.tgraph.edges[q, r] == 2) or
                                                                         ((self.tgraph.edges[p, r] == 1)
                                                                          and (self.tgraph.edges[r, p] == 1))
                                                                         or ((self.tgraph.edges[q, r] == 1) and
                                                                             (self.tgraph.edges[r, q] == 1)))]
        r_list = [list(r) for r in itertools.combinations(r_list, set_size)]
        x = self.ts_data_df[self.ts_names[p]]
        y = self.ts_data_df[self.ts_names[q]]

        for rs in r_list:
            z = dict()
            for r in rs:
                z[self.ts_names[r]] = self.ts_data_df[self.ts_names[r]]

            cmi_pval, cmi_val = self.cmi_dsr(x, y, z, p_value=self.rank_using_p_value)
            if self.rank_using_p_value:
                v_list.append(cmi_pval)
            else:
                v_list.append(cmi_val)
        if v_list:
            return p, q, v_list, r_list

    def ts_rank_cmi_sep_set_parallel(self, set_size):
        """
        rank pairs of time series based on the estimation of ctmi between each pair of connected time series
        :param set_size: cardinality of the set of neighbors
        :return: ranking of each pair of connected time series based ctmi
        """
        list_adj, list_num_adj = self.tgraph.search_adj_all()
        p_list = [p for p in range(len(list_num_adj)) if list_num_adj[p] > set_size]
        q_list = [list_adj[p] for p in p_list]
        p_list = [p_list[p] for p in range(len(p_list)) for _ in q_list[p]]
        q_list = [q for sublist in q_list for q in sublist]
        pq_list = [(p, q) for p, q in zip(p_list, q_list)]
        temp_pq = pq_list.copy()
        temp_p = p_list.copy()
        temp_q = q_list.copy()
        for pq in range(len(temp_pq)):
            if (temp_pq[pq][1], temp_pq[pq][0]) in pq_list:
                pq_list.remove((temp_pq[pq][0], temp_pq[pq][1]))
                p_list.remove(temp_p[pq])
                q_list.remove(temp_q[pq])
        del temp_pq, temp_p, temp_q

        res = Parallel(n_jobs=self.num_processor)(delayed(self._ts_cmi_sep_set_pq)(p, q, set_size) for p, q in
                                                  zip(p_list, q_list))
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

    def ts_find_sep_set(self):
        """
        find the most contributing separation set (if it exists) between each pair of time series
        """
        if self.verbose:
            print("######################################")
            print("Temporal Skeletion Speperation")
            print("######################################")

        for set_size in range(1, self.tgraph.d-1):
            ranks = self.ts_rank_cmi_sep_set_parallel(set_size)
            if self.verbose:
                print("Ranking:")
                print("p: "+str(ranks.elem_p))
                print("q: " + str(ranks.elem_q))
                print("r: " + str(ranks.elem_r))
                print("val: " + str(ranks.val))
            for p, q, r_set, cmi in zip(ranks.elem_p, ranks.elem_q, ranks.elem_r, ranks.val):
                test = (self.tgraph.edges[p, q] != 0)
                for r in r_set:
                    if not test:
                        break
                    test = test and ((self.tgraph.edges[q, r] != 0) or (self.tgraph.edges[p, r] != 0))
                if test:
                    mi = self.ts_mi_array[p, q]

                    if self.p_value != self.rank_using_p_value:
                        x = self.ts_data_df[self.ts_names[p]]
                        y = self.ts_data_df[self.ts_names[q]]
                        z = dict()
                        for r in r_set:
                            z[self.ts_names[r]] = self.ts_data_df[self.ts_names[r]]
                        cmi, _ = self.cmi_dsr(x, y, z, p_value=self.p_value)
                    if self.verbose:
                        print("p=" + self.ts_names[p] + "; q=" + self.ts_names[q] + "; r=" + str(r_set) + "; I(p,q|r)=" + "{: 0.5f}".format(
                            cmi) + "; I(p,q)=" + "{: 0.5f}".format(mi), end=" ")

                    if self.p_value:
                        test = self.sig_lev < cmi
                    else:
                        test = cmi < self.alpha
                    if test:
                        self.ts_cmi_array[p, q] = cmi
                        self.ts_cmi_array[q, p] = cmi
                        self.ts_remove_hom_edges(p, q, self.ts_names_dict_inv[self.ts_names[p]],
                                           self.ts_names_dict_inv[self.ts_names[q]])
                        # self.tgraph.edges[p, q] = 0
                        # self.tgraph.edges[q, p] = 0

                        for r in r_set:
                            self.tgraph.add_sep(q, p, r)
                    else:
                        if self.verbose:
                            print()

    def ts_remove_hom_edges(self, p, q, p_name, q_name):
        if (p < self.tgraph.d) and (q < self.tgraph.d) and (self.ts_names_dict_inv[self.ts_names[p]] == p_name) and \
                (self.ts_names_dict_inv[self.ts_names[q]] == q_name):
            if self.verbose:
                print("=> remove link between " + str(p) + " and " + str(q))
            self.tgraph.edges[p, q] = 0
            self.tgraph.edges[q, p] = 0
            p = p+1
            q = q+1
            self.ts_remove_hom_edges(p, q, p_name, q_name)

    # def is_parent_of(self, p, q, lag):
    #     parent_list_of_q = self.tgraph_dict[self.names[q]]
    #     parent_list_of_p = self.tgraph_dict[self.names[p]]
    #     backward = True
    #     for l in range(self.ts_window_size):
    #         backward = backward and ((self.names[q], -l) not in parent_list_of_p)
    #     return ((self.names[p], -lag) in parent_list_of_q) and backward
    #
    # def is_neighbor_of(self, p, q, lag):
    #     parent_list_of_q = self.tgraph_dict[self.names[q]]
    #     parent_list_of_p = self.tgraph_dict[self.names[p]]
    #     return (p in parent_list_of_q) or (q in parent_list_of_p)

    def cmi_dsr(self, x, y, z_dict=dict(), p_value=False, k=10, sig_samples=10000):
        # todo adapt to different sampling rate
        if isinstance(x, pd.Series):
            x = x.to_frame()
        if isinstance(y, pd.Series):
            y = y.to_frame()
        # if z_dict:
        names_z = [*z_dict.keys()]
        if len(names_z) > 0:
            z = pd.DataFrame()
            for name in names_z:
                if isinstance(z_dict[name], pd.Series):
                    z_dict[name] = z_dict[name].to_frame()
                z[z_dict[name].columns] = z_dict[name].reset_index(drop=True)
        else:
            z = None
        return indep_test(x, y, z, sig_samples=sig_samples, p_value=p_value, measure="cmiknn", k=k)

    # todo
    def ts_dataframe_to_dict(self):
        names = self.names.tolist()
        print(names)
        nlags = self.ts_window_size
        df =pd.DataFrame(self.tgraph.edges, columns=self.ts_names, index=self.ts_names)
        g_dict = dict()
        for name_y in names:
            g_dict[name_y] = []
        for ty in range(nlags):
            for name_y in names:
                t_name_y = df.columns[names.index(name_y) * self.ts_window_size + ty]
                for tx in range(nlags):
                    for name_x in names:
                        t_name_x = df.columns[names.index(name_x) * self.ts_window_size + tx]
                        if df[t_name_y].loc[t_name_x] == 2:
                            if (name_x, tx - ty) not in g_dict[name_y]:
                                g_dict[name_y].append((name_x, tx - ty))
        return g_dict

    # def _ts_cmi_sep_set_pq(self, p, q, set_size):
    #     """
    #     estimate mi between two time series conditioned on each set of neighbors with cardinality equal to set_size
    #     :param p: time series with index p
    #     :param q: time series with index q
    #     :param set_size: cardinality of the set of neighbors
    #     :return: p, q, list if estimated value of ctmi(p,q,r_set), and list of all r_sets
    #     """
    #     v_list = []
    #     r_list = [r for r in range(self.d) if (self.is_parent_of(r, q)) or (self.is_neighbor_of(q, r) and
    #                                                                         not self.is_parent_of(r, q))]
    #     r_list = [list(r) for r in itertools.combinations(r_list, set_size)]
    #
    #     x = self.ts_data_dict[self.names[p]]
    #     y = self.ts_data_dict[self.names[q]]
    #
    #     print(x)
    #     for rs in r_list:
    #         z = dict()
    #         for r in rs:
    #             z[self.names[r]] = self.data_dict[self.names[r]]
    #
    #         cmi_pval, cmi_val = self.cmi(x, y, z)
    #         if self.rank_using_p_value:
    #             v_list.append(cmi_pval)
    #         else:
    #             v_list.append(cmi_val)
    #     if v_list:
    #         return p, q, v_list, r_list
    #
    # def ts_rank_cmi_sep_set_parallel(self, set_size):
    #     """
    #     rank pairs of time series based on the estimation of ctmi between each pair of connected time series
    #     :param set_size: cardinality of the set of neighbors
    #     :return: ranking of each pair of connected time series based ctmi
    #     """
    #     pq_list = [(p, q) for p in range(self.d) for q in range(p, self.d) if self.is_neighbor_of(p, q)]
    #     print(pq_list)
    #     res = Parallel(n_jobs=self.num_processor)(delayed(self._cmi_sep_set_pq)(p, q, set_size) for p, q in
    #                                               pq_list)
    #     ranks = RankingList()
    #     for pq in range(len(res)):
    #         if res[pq] is not None:
    #             if isinstance(res[pq][2], list):
    #                 for r in range(len(res[pq][2])):
    #                     ranks.add(res[pq][0], res[pq][1], res[pq][2][r], res[pq][3][r])
    #             else:
    #                 ranks.add(res[pq][0], res[pq][1], res[pq][2], res[pq][3])
    #     if self.rank_using_p_value:
    #         ranks.sort(descending=True)
    #     else:
    #         ranks.sort(descending=False)
    #     return ranks
    #
    # def ts_find_sep_set(self):
    #     """
    #     find the most contributing separation set (if it exists) between each pair of time series
    #     """
    #     if self.verbose:
    #         print("######################################")
    #         print("Ts Skeletion Speperation")
    #         print("######################################")
    #
    #     for set_size in range(1, (self.graph.d*self.ts_window_size)-1):
    #         ranks = self.rank_cmi_sep_set_parallel(set_size)
    #         if self.verbose:
    #             print("Ranking:")
    #             print("p: "+str(ranks.elem_p))
    #             print("p: " + str(ranks.elem_q))
    #             print("p: " + str(ranks.elem_r))
    #             print("p: " + str(ranks.val))
    #         for p, q, r_set, cmi in zip(ranks.elem_p, ranks.elem_q, ranks.elem_r, ranks.val):
    #             test = self.is_neighbor_of(p, q)
    #             for r in r_set:
    #                 if not test:
    #                     break
    #                 test = test and (self.is_neighbor_of(q, r) or self.is_neighbor_of(p, r))
    #             if test:
    #                 mi = self.ts_mi_array[p, q, lag]
    #
    #                 if self.p_value != self.rank_using_p_value:
    #                     x = self.data_dict[self.names[p]]
    #                     y = self.data_dict[self.names[q]]
    #                     z = dict()
    #                     for r in r_set:
    #                         z[self.names[r]] = self.data_dict[self.names[r]]
    #                     cmi, _ = self.cmi(x, y, z)
    #                 if self.verbose:
    #                     print("p=" + str(p) + "; q=" + str(q) + "; r=" + str(r_set) + "; I(p,q|r)=" + "{: 0.5f}".format(
    #                         cmi) + "; I(p,q)=" + "{: 0.5f}".format(mi), end=" ")
    #
    #                 if self.p_value:
    #                     test = mi < self.sig_lev < cmi
    #                 else:
    #                     test = cmi < self.alpha
    #                 if test:
    #                     self.cmi_array[p, q] = cmi
    #                     self.cmi_array[q, p] = cmi
    #                     if self.verbose:
    #                         print("=> remove link between " + str(p) + " and " + str(q))
    #                     if self.is_parent_of(p, q):
    #                         1
    #                     elif self.is_parent_of(q, p):
    #                         1
    #                     elif self.is_neighbor_of(p, q):
    #                         self.graph.edges[p, q] = 0
    #                         self.graph.edges[q, p] = 0
    #
    #                     for r in r_set:
    #                         self.graph.add_sep(q, p, r)
    #                 else:
    #                     if self.verbose:
    #                         print()

    def noise_based_hidden_counfounders(self):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel

        for p in range(self.tgraph.d):
            for q in range(self.tgraph.d):
                if self.tgraph.edges[p, q] == 2:
                    print(p,q)
                    X = self.ts_data_df[self.ts_names[p]].values.reshape(-1, 1)
                    y = self.ts_data_df[self.ts_names[q]].values.reshape(-1, 1)

                    kernel = DotProduct() + WhiteKernel()
                    gpr = GaussianProcessRegressor(kernel=kernel, random_state = 0).fit(X, y)
                    print(gpr.score(X, y))

    def fit(self):
        """
        run PCTMI
        :return: graph (CPDAG)
        """
        if self.verbose:
            print("#######################################")
            print("########### Starting TPCTMI ###########")
            print("#######################################")

        # initialize skeleton
        self.skeleton_initialize()

        # get separation sets
        self.find_sep_set()

        # orientation
        self.rule_origin_causality()

        test_rp = True
        test_r2 = True
        test_r3 = True
        while test_rp or test_r2 or test_r3:
            test_rp = self.rule_propagation_causality()
            test_r2 = self.rule_2()
            test_r3 = self.rule_3()

        self.rule_commun_confounder_and_causal_chain()
        self.rule_mediator()
        self.rule_proba_raising_principle()

        if self.verbose:
            print("######################################")
            print("Final Results (PCTMI)")
            print("######################################")
            print("Summary Graph:")
            print(self.graph.edges)

        self.check_cycles()

        self.summary_to_temporal_array()
        # self.summary_to_temporal()

        self.ts_find_sep_set()

        self.tgraph_dict = self.ts_dataframe_to_dict()

        if self.verbose:
            print("######################################")
            print("Final Results (TPCTMI)")
            print("######################################")
            print("Temporal Graph:")
        print(self.tgraph.edges)
        return self.tgraph_dict


class FCITMI(CITMI):
    # def __init__(self, series, sig_lev=0.05, p_value=True, rank_using_p_value=False, verbose=True, num_processor=-1,
    #              graphical_optimization=True):
    def __init__(self, series, sig_lev=0.05, lag_max=5, p_value=True, rank_using_p_value=False, verbose=True,
                     num_processor=-1, graphical_optimization=False):
        """
        FCI for time series using TMI and CTMI
        :param series: d-time series (with possibility of different sampling rate)
        :param sig_lev: significance level. By default 0.05
        :param p_value: Use p_value for decision making. By default True
        :param verbose: Print results. By default: True
        :param num_processor: number of processors for parallelization. By default -1 (all)
        """
        # CITMI.__init__(self, series, sig_lev, p_value, rank_using_p_value, verbose, num_processor,
        #                graphical_optimization)
        CITMI.__init__(self, series, sig_lev, lag_max, p_value, rank_using_p_value, verbose, num_processor,
                       graphical_optimization)

    def dag_to_pag(self):
        """
        transform dag to PAG (turn all tails to circle (undetermined))
        Graph structure
        0: no edge
        1: a tail -
        2: arrow head ->
        3: Undetermined -o
        """
        self.graph.edges[self.graph.edges == 1] = 3

    def _find_shortest_directed_paths_util(self, i, j, visited, path, all_path):
        """
        sub function of _find_shortest_directed_paths
        :param i: index of time series
        :param j: index of time series
        :param visited: list of visited nodes
        :param path: current path
        :param all_path: list of all discovered paths
        """
        # Mark the current node as visited and store in path
        visited[i] = True
        path.append(i)

        # If current vertex is same as destination, then print
        # current path[]
        if i == j:
            if len(path) > 2:
                all_path.append(path.copy()[1:-1])
                return path
        else:
            # If current vertex is not destination
            # Recur for all the vertices child of this vertex
            child_i = np.where(self.graph.edges[i, :] == 2)[0]
            for k in child_i:
                if not visited[k]:
                    self._find_shortest_directed_paths_util(k, j, visited, path, all_path)

        # Remove current vertex from path[] and mark it as unvisited
        path.pop()
        visited[i] = False

    def _find_shortest_directed_paths(self, i, j):
        """
        find shortest directed path between time series of index i and time series of index j
        :param i: index of time series
        :param j: index of time series
        :return: all directed paths from i to j
        """
        # Mark all the vertices as not visited
        visited = [False] * self.graph.d

        # Create an array to store paths
        path = []
        all_path = []

        # Call the recursive helper function to print all paths
        self._find_shortest_directed_paths_util(i, j, visited, path, all_path)
        return all_path


    # todo
    def _find_discriminating_paths_util(self, i, j, visited, path, all_path):
        """
        sub function of _find_shortest_directed_paths
        :param i: index of time series
        :param j: index of time series
        :param visited: list of visited nodes
        :param path: current path
        :param all_path: list of all discovered paths
        """
        # Mark the current node as visited and store in path
        visited[i] = True
        path.append(i)
        path.append(j)
        all_path.append(path.copy())
        path.pop()

        i_child = (self.graph.edges[i, :] == 2)
        i_parent = (self.graph.edges[:, i] == 2)
        j_adj1 = (self.graph.edges[:, j] != 0)
        j_adj2 = (self.graph.edges[j, :] != 0)
        next_i = np.where([a and b and c and d for a, b, c, d in zip(i_child, i_parent, j_adj1, j_adj2)])[0]

        for k in next_i:
            if not visited[k]:
                if (self.graph.edges[k, j] == 2) and (self.graph.edges[j, k] == 1):
                    self._find_shortest_directed_paths_util(k, j, visited, path, all_path)
                else:
                    visited[k] = True
                    path.append(k)
                    path.append(j)
                    all_path.append(path.copy())
                    path.pop()

        # Remove current vertex from path[] and mark it as unvisited
        path.pop()
        visited[i] = False

    def _find_discriminating_paths(self, i, j):
        """
        find discriminating  path between time series of index i and time series of index j
        :param i: index of time series
        :param j: index of time series
        :return: all discriminating paths from i to j
        """
        # Mark all the vertices as not visited
        visited = [False] * self.graph.d

        # Create an array to store paths
        path = [i]
        all_path = []

        i_child = (self.graph.edges[i, :] == 2)
        i_non_parent = (self.graph.edges[:, i] != 2)
        j_adj1 = (self.graph.edges[:, j] != 0)
        j_adj2 = (self.graph.edges[j, :] != 0)
        first_next_i = np.where([a and b and c and d for a, b, c, d in zip(i_child, i_non_parent, j_adj1, j_adj2)])[0]
        print(first_next_i)
        for f in first_next_i:
            # Call the recursive helper function to print all paths
            self._find_discriminating_paths_util(f, j, visited, path, all_path)

        return all_path

    def _find_ancestors_util(self, i, visited, path, all_path):
        """
        sub function of _find_ancestors
        :param i: index of time series
        :param visited: list of visited nodes
        :param path: current path
        :param all_path: list of all discovered paths
        :return:
        """
        # Mark the current node as visited and store in path
        visited[i] = True
        path.append(i)

        # If current vertex is same as destination, then print
        # current path[]
        parent_i = np.where(self.graph.edges[:, i] == 2)[0]
        if len(parent_i) == 0:
            if len(path) > 1:
                all_path.append(path.copy()[1:])
        else:
            # If current vertex is not destination
            # Recur for all the vertices child of this vertex
            for k in parent_i:
                if not visited[k]:
                    self._find_ancestors_util(k, visited, path, all_path)

        # Remove current vertex from path[] and mark it as unvisited
        path.pop()
        visited[i] = False

    def _find_ancestors(self, i):
        """
        find ancestors of time series of index i
        :param i: index if time series
        :return: list of ancestors
        """
        # Mark all the vertices as not visited
        visited = [False] * self.graph.d

        # Create an array to store paths
        path = []
        all_path = []

        # Call the recursive helper function to print all paths
        self._find_ancestors_util(i, visited, path, all_path)
        ancestors = [item for path in all_path for item in path]
        ancestors = list(set(ancestors))
        return ancestors

    def _find_possible_d_sep_ij_util(self, i, j, v, before_v, anc_i, anc_j, visited, path, possible_d_sep):
        # Mark the current node as visited and store in path
        visited[v] = True
        path.append(v)

        print(v, before_v, possible_d_sep)
        if (before_v != v) and (self.graph.edges[before_v, v] != 2):
            if len(path) > 1:
                print("one")
                possible_d_sep.append(v)
        else:
            # If current vertex is not destination
            # Recur for all the vertices child of this vertex
            adj_v = np.where(self.graph.edges[:, v] == 2)[0]
            adj_v = [k for k in adj_v if not visited[k]]
            print(visited)
            for k in adj_v:
                if not visited[k]:
                    if (v != i) and (v != j):
                        print("two")
                        possible_d_sep.append(v)
                    self._find_possible_d_sep_ij_util(i, j, k, v, anc_i, anc_j, visited, path, possible_d_sep)

    def _find_possible_d_sep_ij(self, i, j):
        """
        :param i: index of time series
        :param j: index of time series
        :return: all possible d-sep if i and j
        """
        anc_i = self._find_ancestors(i)
        anc_j = self._find_ancestors(j)

        # Mark all the vertices as not visited
        visited = [False] * self.graph.d

        # Create an array to store paths
        path = []
        possible_d_sep = []

        # Call the recursive helper function to print all paths
        self._find_possible_d_sep_ij_util(i, j, i, i, anc_i, anc_j, visited, path, possible_d_sep)
        return possible_d_sep

    def _cmi_possible_d_sep_ij(self, p, q, set_size):
        """
        estimate ctmi between two time series conditioned on each possible-d-set with cardinality equal to set_size
        :param i: time series with index i
        :param j: time series with index j
        :param set_size: cardinality of the set of neighbors
        :return: i, j, list of estimated values of ctmi(p,q,possible-d-set), and list of all possible-d-sets
        """
        v_list = []
        k_list = self._find_possible_d_sep_ij(p, q)
        k_list = [list(k) for k in itertools.combinations(k_list, set_size)]
        print(p, q, k_list)
        if self.adaptive_window:
            # x = self.series[self.names[p]]
            # y = self.series[self.names[p]]
            x = window_representation(self.series[self.names[p]],
                                      windows_size=self.window_matrix[self.names[p]].loc[self.names[p]])
            y = window_representation(self.series[self.names[q]],
                                      windows_size=self.window_matrix[self.names[q]].loc[self.names[q]])
        else:
            x = self.data_dict[self.names[p]]
            y = self.data_dict[self.names[q]]

        for ks in k_list:
            z = dict()
            for k in ks:
                if self.adaptive_window:
                    z[self.names[k]] = self.series[self.names[k]].dropna()
                else:
                    z[self.names[k]] = self.data_dict[self.names[k]]
            if self.graphical_optimization:
                cmi_pval, cmi_val = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
                                         gamma_matrix=self.gamma_matrix, graph=self.graph.edges,
                                         p_value=self.p_value, instantaneous_dict=self.instantaneous_dict)
            else:
                # cmi_pval, cmi_val = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
                #               gamma_matrix=self.gamma_matrix, p_value=self.p_value, instantaneous_dict=self.instantaneous_dict)
                # cmi_pval, cmi_val = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
                #                          gamma_matrix=self.gamma_matrix, p_value=self.p_value, instantaneous_dict=self.instantaneous_dict)
                cmi_pval, cmi_val = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
                                         gamma_matrix=self.gamma_matrix, p_value=self.rank_using_p_value,
                                         instantaneous_dict=self.instantaneous_dict)
            v_list.append(cmi_pval)
        if v_list:
            return p, q, v_list, k_list

    # def _cmi_possible_d_sep_ij(self, p, q, set_size):
    #     """
    #     estimate ctmi between two time series conditioned on each set of neighbors with cardinality equal to set_size
    #     :param p: time series with index p
    #     :param q: time series with index q
    #     :param set_size: cardinality of the set of neighbors
    #     :return: p, q, list if estimated value of ctmi(p,q,r_set), and list of all r_sets
    #     """
    #     v_list = []
    #     k_list = self._find_possible_d_sep_ij(p, q)
    #     r_list = [list(k) for k in itertools.combinations(k_list, set_size)]
    #
    #     r_list_temp = r_list.copy()
    #     # if set_size == 1:
    #     for rs in r_list_temp:
    #         print(rs)
    #         print(all(elem >= self.d for elem in rs))
    #         if all(elem >= self.d for elem in rs):
    #             r_list.remove(rs)
    #     del r_list_temp
    #
    #     if self.adaptive_window:
    #         x = window_representation(self.series[self.names[p]], windows_size=self.window_matrix[self.names[p]].loc[self.names[p]])
    #         y = window_representation(self.series[self.names[q]], windows_size=self.window_matrix[self.names[q]].loc[self.names[q]])
    #     else:
    #         x = self.data_dict[self.names[p]]
    #         y = self.data_dict[self.names[q]]
    #
    #     for rs in r_list:
    #         z = dict()
    #         for r in rs:
    #             if self.adaptive_window:
    #                 # select and drop NA
    #                 z[self.names[r]] = self.series[self.names[r]].dropna()
    #             else:
    #                 z[self.names[r]] = self.data_dict[self.names[r]]
    #         if self.graphical_optimization:
    #             # cmi_pval, cmi_val = gctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
    #             #                           gamma_matrix=self.gamma_matrix, p_value=self.rank_using_p_value,
    #             #                           graph=self.graph.edges)
    #             cmi_pval, cmi_val = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
    #                                       gamma_matrix=self.gamma_matrix, graph=self.graph.edges,
    #                                      p_value=self.rank_using_p_value, instantaneous_dict=self.instantaneous_dict)
    #         else:
    #             cmi_pval, cmi_val = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
    #                                      gamma_matrix=self.gamma_matrix, p_value=self.rank_using_p_value,
    #                                      instantaneous_dict=self.instantaneous_dict)
    #
    #         if self.rank_using_p_value:
    #             v_list.append(cmi_pval)
    #         else:
    #             v_list.append(cmi_val)
    #     if v_list:
    #         return p, q, v_list, r_list


    def rank_possible_d_sep_parallel(self, set_size):
        """
        rank pairs of connected time series conditioned of their possible-d-sep based on the estimation of ctmi
        :param set_size: cardinality of the possible-d-sep
        :return: ranking of each pair of connected time series based ctmi
        """
        list_adj, list_num_adj = self.graph.search_adj_all()
        i_list = [[i]*list_num_adj[i] for i in range(len(list_num_adj)) if list_num_adj[i] > 0]
        i_list = [i for sublist in i_list for i in sublist]
        j_list = [list_adj[j] for j in range(len(list_num_adj)) if list_num_adj[j] > 0]
        j_list = [j for sublist in j_list for j in sublist]

        # res = Parallel(n_jobs=self.num_processor)(delayed(self._cmi_possible_d_sep_ij)(i, j, set_size) for i, j in
        #                                           zip(i_list, j_list))
        res = []
        for p, q in zip(i_list, j_list):
            res.append(self._cmi_possible_d_sep_ij(p, q, set_size))

        ranks = RankingList()
        for ij in range(len(res)):
            if res[ij] is not None:
                if isinstance(res[ij][2], list):
                    for k in range(len(res[ij][2])):
                        ranks.add(res[ij][0], res[ij][1], res[ij][2][k], res[ij][3][k])
                else:
                    ranks.add(res[ij][0], res[ij][1], res[ij][2], res[ij][3])
        if self.rank_using_p_value:
            ranks.sort(descending=True)
        else:
            ranks.sort(descending=False)
        return ranks

    def find_d_sep(self):
        """
        find the most contributing d sep (if it exists) between each pair of time series
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("d-seperation")
            print("######################################")
        test_remove_links = False

        for set_size in range(1, self.graph.d-1):
            ranks = self.rank_possible_d_sep_parallel(set_size)
            for i, j, ks, cmi in zip(ranks.elem_p, ranks.elem_q, ranks.elem_r, ranks.val):
                test = (self.graph.edges[i, j] != 0)
                for k in ks:
                    if not test:
                        break
                    test = test and ((self.graph.edges[j, k] != 0) or (self.graph.edges[i, k] != 0))
                if test:
                    mi = self.mi_array[i, j]
                    if self.verbose:
                        print("i=" + str(i) + "; j=" + str(j) + "; z=" + str(ks) + "; I(i,j|z)=" + "{: 0.5f}".format(
                            cmi) + "; I(i,j)=" + "{: 0.5f}".format(mi), end=" ")

                    if self.p_value:
                        test = mi < self.sig_lev < cmi
                    else:
                        test = cmi < self.alpha
                    if test:
                        test_remove_links = True
                        self.cmi_array[i, j] = cmi
                        self.cmi_array[j, i] = cmi
                        if self.verbose:
                            print("=> remove link between " + str(i) + " and " + str(j))
                        self.graph.edges[i, j] = 0
                        self.graph.edges[j, i] = 0
                        for k in ks:
                            self.graph.add_sep(j, i, k)
                    else:
                        if self.verbose:
                            print()
        return test_remove_links


    # def find_d_sep(self):
    #     """
    #     find the most contributing separation set (if it exists) between each pair of time series
    #     """
    #     if self.verbose:
    #         print("######################################")
    #         print("d-seperation")
    #         print("######################################")
    #     test_remove_links = False
    #
    #
    #     print("max set size = " + str(self.graph.d-1))
    #     for set_size in range(1, self.graph.d-1):
    #         ranks = self.rank_possible_d_sep_parallel(set_size)
    #         if self.verbose:
    #             print("Ranking:")
    #             print("p: "+str(ranks.elem_p))
    #             print("p: " + str(ranks.elem_q))
    #             print("p: " + str(ranks.elem_r))
    #             print("p: " + str(ranks.val))
    #         for p, q, r_set, cmi in zip(ranks.elem_p, ranks.elem_q, ranks.elem_r, ranks.val):
    #             test = (self.graph.edges[p, q] != 0)
    #             for r in r_set:
    #                 if not test:
    #                     break
    #                 # test = test and ((self.graph.sep[p, r, q] == 0) and (self.graph.sep[q, r, p] == 0))
    #             if test:
    #                 mi = self.mi_array[p, q]
    #
    #                 if self.p_value != self.rank_using_p_value:
    #                     if self.adaptive_window:
    #                         x = window_representation(self.series[self.names[p]],
    #                                                   windows_size=self.window_matrix[self.names[p]].loc[self.names[p]])
    #                         y = window_representation(self.series[self.names[q]],
    #                                                   windows_size=self.window_matrix[self.names[q]].loc[self.names[q]])
    #                     else:
    #                         x = self.data_dict[self.names[p]]
    #                         y = self.data_dict[self.names[q]]
    #
    #                     z = dict()
    #                     for r in r_set:
    #                         if self.adaptive_window:
    #                             # select and drop NA
    #                             z[self.names[r]] = self.series[self.names[r]].dropna()
    #                         else:
    #                             z[self.names[r]] = self.data_dict[self.names[r]]
    #                     if self.graphical_optimization:
    #                         # cmi, _ = gctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
    #                         #                gamma_matrix=self.gamma_matrix, p_value=self.p_value, graph=self.graph.edges)
    #                         cmi_pval, cmi_val = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
    #                                                  gamma_matrix=self.gamma_matrix, graph=self.graph.edges,
    #                                                  p_value=self.rank_using_p_value,
    #                                                  instantaneous_dict=self.instantaneous_dict)
    #                     else:
    #                         cmi, _ = ctmi(x, y, z, self.names[p], self.names[q], self.sampling_rate,
    #                                       gamma_matrix=self.gamma_matrix, p_value=self.p_value,
    #                                       instantaneous_dict=self.instantaneous_dict)
    #                 if self.verbose:
    #                     print("p=" + str(p) + "; q=" + str(q) + "; r=" + str(r_set) + "; I(p,q|r)=" + "{: 0.5f}".format(
    #                         cmi) + "; I(p,q)=" + "{: 0.5f}".format(mi), end=" ")
    #
    #                 if self.p_value:
    #                     test = mi < self.sig_lev < cmi
    #                 else:
    #                     test = cmi < self.alpha
    #                 if test:
    #                     test_remove_links = True
    #                     self.cmi_array[p, q] = cmi
    #                     self.cmi_array[q, p] = cmi
    #                     if self.verbose:
    #                         print("=> remove link between " + str(p) + " and " + str(q))
    #                     self.graph.edges[p, q] = 0
    #                     self.graph.edges[q, p] = 0
    #
    #                     for r in r_set:
    #                         self.graph.add_sep(q, p, r)
    #                         self.biggamma[p,q,r] = self.gamma_matrix[self.names[p]].loc[self.names[r]]
    #                         self.biggamma[q,p,r] = self.gamma_matrix[self.names[q]].loc[self.names[r]]
    #                 else:
    #                     if self.verbose:
    #                         print()


    def remove_orientation(self):
        """
        turn all vertex into undetermined vertex
        """
        if self.verbose:
            print("######################################")
            print("Remove orientation")
            print("######################################")
        for i in range(self.graph.d):
            for j in range(self.graph.d):
                if i != j:
                    if self.graph.edges[i, j] != 0:
                        self.graph.edges[i, j] = 3

    def rule_origin_causality(self):
        """
        rule 0 (origin of causality) from FCI
        """
        if self.verbose:
            print("######################################")
            print("Rule Origin of Causality")
            print("######################################")
        for p in range(self.d):
            for q in range(p+1, self.d):
                if self.graph.edges[p, q] == 0:
                    sep = np.where(self.graph.sep[p, q, :] == 1)[0]
                    for r in range(self.d):
                        if (r != p) and (r != q) and (r not in sep):
                            if (self.graph.edges[q, r] == 3) and (self.graph.edges[p, r] == 3) and (
                                    self.graph.edges[r, q] == 3) and (self.graph.edges[r, p] == 3):
                                self.graph.edges[p, r] = 2
                                self.graph.edges[q, r] = 2

    # def _oc_pq(self, p, q):
    #     """
    #     estimate ctmi between two time series conditioned of their sep set + each non oriented connected neighbor
    #     :param p: index of time series
    #     :param q: index of time series
    #     :return: p, q, the most contributing node and the MI conditioned on sep set and the most contributing node
    #     """
    #     v_list = []
    #     k_list = [k for k in range(self.graph.d) if (k != p) and (k != q)
    #               and ((self.graph.edges[q, k] == 3) and (self.graph.edges[p, k] == 3)
    #                    and (self.graph.edges[k, q] == 3) and (self.graph.edges[k, p] == 3)
    #                    and (self.graph.sep[p, q, k] == 0))]
    #     if len(k_list) > 0:
    #         x = self.data_dict[self.names[p]]
    #         y = self.data_dict[self.names[q]]
    #         sep = np.where(self.graph.sep[p, q, :] == 1)[0]
    #         for k in k_list:
    #             if k not in sep:
    #                 z = self.data_dict[self.names[k]]
    #                 print(z.shape)
    #                 cmi_pval, cmi_val = i_ctmi(x, y, z, self.names[p], self.names[q], self.names[k], self.sampling_rate,
    #                            p_value=self.p_value)
    #                 v_list.append(cmi_pval)
    #     if v_list:
    #         if self.p_value:
    #             idx = int(np.argmax(v_list))
    #         else:
    #             idx = int(np.argmin(v_list))
    #         return p, q, v_list[idx], k_list[idx]
    #
    # def rank_oc_parallel(self):
    #     """
    #     rank unsheilded triples based on the estimation of ctmi
    #     :return: ranking of each unsheilded triple based ctmi
    #     """
    #     p_list = []
    #     q_list = []
    #     for p in range(self.d):
    #         for q in range(p+1, self.d):
    #             if (self.graph.edges[p, q] == 0) and (self.graph.edges[q, p] == 0):
    #                 p_list.append(p)
    #                 q_list.append(q)
    #     res = Parallel(n_jobs=self.num_processor)(delayed(self._oc_pq)(p, q) for p, q in zip(p_list, q_list))
    #     ranks = RankingList()
    #     for pq in range(len(res)):
    #         if res[pq] is not None:
    #             ranks.add(res[pq][0], res[pq][1], res[pq][2], res[pq][3])
    #     if self.p_value:
    #         ranks.sort(descending=False)
    #     else:
    #         ranks.sort(descending=True)
    #     return ranks
    #
    # def rule_origin_causality(self):
    #     """
    #     rule 0 (origin of causality) from FCI
    #     """
    #     if self.verbose:
    #         print("######################################")
    #         print("Rule Origin of Causality")
    #         print("######################################")
    #
    #     ranks = self.rank_oc_parallel()
    #     for p, q, k, cmi in zip(ranks.elem_p, ranks.elem_q, ranks.elem_r, ranks.val):
    #         if (self.graph.edges[q, k] == 3) and (self.graph.edges[p, k] == 3) \
    #                 and (self.graph.edges[k, q] == 3) and (self.graph.edges[k, p] == 3):
    #             sep = np.where(self.graph.sep[p, q, :] == 1)[0]
    #             print("sep = " + str(sep))
    #             # if len(sep) > 0:
    #             #     mi = self.cmi_array[p, q]
    #             # else:
    #             mi = self.mi_array[p, q]
    #             if k not in sep:
    #                 if self.verbose:
    #                     print("p=" + str(p) + "; q=" + str(q) + "; r=" + str(
    #                         k) + "; s=" + str(
    #                         sep) + "; I(p,q|r,s)=" + "{: 0.5f}".format(
    #                         cmi) + "; I(p,q|s)=" + "{: 0.5f}".format(mi), end=" ")
    #                 if self.p_value:
    #                     test = cmi < mi
    #                 else:
    #                     test = mi - cmi < 0
    #                 if test:
    #                     if self.verbose:
    #                         print("=> orient " + str(p) + " -> " + str(k) + " <- " + str(q))
    #                     self.graph.edges[p, k] = 2
    #                     self.graph.edges[q, k] = 2
    #                 else:
    #                     if self.verbose:
    #                         print()

    # def rule_origin_causality(self):
    #     """
    #     rule 0 (origin of causality) from FCI
    #     """
    #     # todo parallelization
    #     if self.verbose:
    #         print("######################################")
    #         print("Rule Origin of Causality")
    #         print("######################################")
    #
    #     for i in range(self.graph.d):
    #         for j in range(i + 1, self.graph.d):
    #             if (self.graph.edges[i, j] == 0) and (self.graph.edges[j, i] == 0):
    #                 k_list = [k for k in range(self.graph.d) if (k != i) and (k != j)
    #                           and ((self.graph.edges[j, k] == 3) and (self.graph.edges[i, k] == 3))]
    #                 if len(k_list) > 0:
    #                     x = self.data_dict[self.names[i]]
    #                     y = self.data_dict[self.names[j]]
    #                     sep = np.where(self.graph.sep[i, j, :] == 1)[0]
    #                     print(str(i), str(j) + "sep = " + str(sep))
    #                     if len(sep) > 0:
    #                         mi = self.cmi_array[i, j]
    #                     else:
    #                         mi = self.mi_array[i, j]
    #                     for k in k_list:
    #                         if k not in sep:
    #                             sep_k = sep.tolist() + [k]
    #                             z = dict()
    #                             for name_k in self.names[sep_k]:
    #                                 z[name_k] = self.data_dict[name_k]
    #                             # cmi_pval, cmi_val = ctmi(x, y, z, self.names[i], self.names[j], self.sampling_rate,
    #                             #            gamma_matrix=self.gamma_matrix, p_value=self.p_value, mission="ictmi")
    #                             cmi_pval, cmi_val = i_ctmi(x, y, z, self.names[i], self.names[j], self.names[k],
    #                                                        self.sampling_rate,
    #                                                        p_value=self.p_value)
    #                             if self.verbose:
    #                                 print("i=" + str(i) + "; j=" + str(j) + "; z=" + str(
    #                                     k) + "; u=" + str(
    #                                     sep) + "; I(i,j|u,z)=" + "{: 0.5f}".format(
    #                                     cmi_pval) + "; I(i,j|u)=" + "{: 0.5f}".format(mi), end=" ")
    #                             if self.p_value:
    #                                 test = cmi_pval < mi
    #                             else:
    #                                 test = mi - cmi_pval < 0
    #                             if test:
    #                                 if self.verbose:
    #                                     print("=> orient " + str(i) + " -> " + str(k) + " <- " + str(j))
    #                                 self.graph.edges[i, k] = 2
    #                                 self.graph.edges[j, k] = 2
    #                             else:
    #                                 if self.verbose:
    #                                     print()

    def rule_temporal_priority_within_time_series(self):
        for i in range(self.graph.d):
            if self.graph.edges[i, i] == 3:
                self.graph.edges[i, i] = 1

    def rule_propagation_causality(self):
        """
        rule 1 from FCI
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule Propagation of Causality")
            print("######################################")

        test_find_orientation = False

        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 0) and (self.graph.edges[j, i] == 0):
                    k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and
                              (((self.graph.edges[j, k] == 2) and
                                (self.graph.edges[k, j] != 0) and (self.graph.edges[k, i] != 0) and
                                (self.graph.edges[i, k] == 3)) or ((self.graph.edges[i, k] == 2) and
                                                                   (self.graph.edges[k, i] != 0) and
                                                                   (self.graph.edges[k, j] != 0) and
                                                                   (self.graph.edges[j, k] == 3)))]
                    if len(k_list) > 0:
                        test_find_orientation = True
                        for k in k_list:
                            if self.graph.edges[i, k] == 2:
                                if self.verbose:
                                    print(str(i) + "*->" + str(k) + "°-*" + str(j), end=" ")
                                    print("=> orient " + str(i) + "*-> " + str(k) + " -> " + str(j))
                                self.graph.edges[k, j] = 2
                                self.graph.edges[j, k] = 1
                            else:
                                if self.verbose:
                                    print(str(j) + "*->" + str(k) + "°-*" + str(i), end=" ")
                                    print("=> orient " + str(j) + "*-> " + str(k) + " -> " + str(i))
                                self.graph.edges[k, i] = 2
                                self.graph.edges[i, k] = 1
        return test_find_orientation

    def rule_2(self):
        """
        rule 2 from FCI
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule 3")
            print("######################################")
        test_find_orientation = False

        for i in range(self.graph.d):
            j_list = np.where(self.graph.edges[i, :] == 3)[0].tolist()
            if i in j_list:
                j_list.remove(i)
            for j in j_list:
                shortest_directed_path = self._find_shortest_directed_paths(i, j)
                if len(shortest_directed_path) > 0:
                    self.graph.edges[i, j] = 2
                    test_find_orientation = True
                    if self.verbose:
                        print_path = '*->'.join(map(str, shortest_directed_path[0]))
                        print(str(i)+"*-0"+str(j)+" and "+str(i) + "*->" + print_path + "*->" + str(j), end=" ")
                        print("=> orient " + str(i) + "*->" + str(j))
        return test_find_orientation

    def rule_3(self):
        """
        rule 3 from FCI
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule 3")
            print("######################################")

        test_find_orientation = False

        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 0) and (self.graph.edges[j, i] == 0):
                    colliders = [k for k in range(self.graph.d) if (k != i) and (k != j) and (
                            (self.graph.edges[j, k] == 2) and (self.graph.edges[i, k] == 2))]
                    k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and (
                            (self.graph.edges[j, k] == 3) and (self.graph.edges[i, k] == 3))]
                    if len(colliders) > 0 and len(k_list) > 0:
                        for c in colliders:
                            for k in k_list:
                                if self.graph.edges[k, c] == 3:
                                    test_find_orientation = True
                                    self.graph.edges[k, c] = 2
                                    if self.verbose:
                                        print(str(i) + "*->" + str(c) + "<-*" + str(j) + " and " + str(i) + "*-0" +
                                              str(k) + "0-*" + str(j) + " and " + str(k) + "*-0" + str(c),
                                              end=" ")
                                        print("=> orient " + str(k) + "*->" + str(c))
        return test_find_orientation

    def rule_4(self):
        """
        rule 4 from FCI
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule 4")
            print("######################################")

        test_find_orientation = False

        for i in range(self.graph.d):
            for j in range(self.graph.d):
                if (i != j and self.graph.edges[i, j] == 0) and (self.graph.edges[j, i] == 0):
                    discriminating_paths = self._find_discriminating_paths(i, j)
                    for dp in discriminating_paths:
                        k = dp[-2]
                        if self.graph.edges[j, k] == 3:
                            self.graph.edges[j, k] = 1
                            self.graph.edges[k, j] = 2
                        else:
                            self.graph.edges[j, k] = 2
                            self.graph.edges[k, j] = 2
                            s = dp[-3]
                            self.graph.edges[s, k] = 2
                            self.graph.edges[k, s] = 2

        return test_find_orientation

    def _find_uncovered_path_util(self, i, j, i_2, i_1, visited, path, all_path):
        """
        sub function of _find_uncovered_path
        :param i: index of time series
        :param j: index of time series
        :param i_2: index of time series at before the previous iteration
        :param i_1: index of time series at the previous iteration
        :param visited: list of visited nodes
        :param path: current path
        :param all_path: list of all discovered paths
        :return:
        """
        # Mark the current node as visited and store in path
        visited[i] = True
        path.append(i)

        # If current vertex is same as destination, then print
        # current path[]
        if i == j:
            if len(path) > 2:
                if len(path) == 3:
                    print(i, i_2)
                    print(path)
                    print(self.graph.edges[i, i_2])
                    if self.graph.edges[i, i_2] == 0:
                        all_path.append(path.copy())
                else:
                    all_path.append(path.copy())
        else:
            # If current vertex is not destination
            # Recur for all the vertices child of this vertex
            child_i = np.where(self.graph.edges[i, :] != 0)[0]
            for k in child_i:
                if not visited[k]:
                    if len(path) > 2:
                        if self.graph.edges[i, i_2] == 0:
                            self._find_uncovered_path_util(k, j, i_1, i, visited, path, all_path)
                    elif len(path) == 2:
                        self._find_uncovered_path_util(k, j, i_2, i, visited, path, all_path)
                    else:
                        self._find_uncovered_path_util(k, j, i_2, i_1, visited, path, all_path)

        # Remove current vertex from path[] and mark it as unvisited
        path.pop()
        visited[i] = False

    def _find_uncovered_path(self, i, j):
        """
        find uncovered path between time series of index i and time series of index j
        :param i: index of time series
        :param j: index of time series
        :return: all uncovered paths from i to j
        """
        # Mark all the vertices as not visited
        visited = [False] * self.graph.d

        # Create an array to store paths
        path = []
        all_uncovered_path = []

        # Call the recursive helper function to print all paths
        self._find_uncovered_path_util(i, j, i, i, visited, path, all_uncovered_path)
        return all_uncovered_path

    def _is_circle_path(self, path):
        """
        check if path is a circle path
        :param path: any path in the graph
        :return: bool
        """
        test_list = []
        for p in range(len(path)-1):
            test_list.append((self.graph.edges[path[p], path[p+1]] == 3) and
                             (self.graph.edges[path[p+1], path[p]] == 3))
        return all(test_list)

    def rule_5(self):
        """
        rule 5 from FCI
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule 5")
            print("######################################")

        test_find_orientation = False
        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 3) and (self.graph.edges[j, i] == 3):
                    uncovered_path_list = self._find_uncovered_path(i, j)
                    if len(uncovered_path_list) > 0:
                        for ucp in uncovered_path_list:
                            if self._is_circle_path(ucp[1:-1]):
                                if (self.graph.edges[ucp[0], ucp[-2]] == 0) and \
                                        (self.graph.edges[ucp[-1], ucp[1]] == 0):
                                    test_find_orientation = True
                                    if self.verbose:
                                        print(str(i) + "0-0" + str(j) + " and found an uncovered path", end=" ")
                                        print("=> orient " + str(i) + "- " + str(j))
                                    self.graph.edges[i, j] = 1
                                    self.graph.edges[j, i] = 1
                                    for p in range(len(ucp)-1):
                                        if self.verbose:
                                            print(str(ucp[p]) + "0-0" + str(ucp[p+1]), end=" ")
                                            print("=> orient " + str(ucp[p]) + "- " + str(ucp[p+1]))
                                        self.graph.edges[ucp[p], ucp[p + 1]] = 1
                                        self.graph.edges[ucp[p + 1], ucp[p]] = 1
        return test_find_orientation

    def rule_6(self):
        """
        rule 6 from FCI
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule 6")
            print("######################################")

        test_find_orientation = False

        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and
                          (((self.graph.edges[j, k] == 3) and (self.graph.edges[k, j] != 0) and
                            (self.graph.edges[k, i] == 1) and (self.graph.edges[i, k] == 1)) or
                           ((self.graph.edges[i, k] == 3) and (self.graph.edges[k, i] != 0) and
                            (self.graph.edges[k, j] == 1) and (self.graph.edges[j, k] == 1)))]
                if len(k_list) > 0:
                    test_find_orientation = True
                    for k in k_list:
                        if self.graph.edges[j, k] == 3:
                            if self.verbose:
                                print(str(i) + "-" + str(k) + "0-*" + str(j), end=" ")
                                print("=> orient " + str(i) + "- " + str(k) + " -* " + str(j))
                            self.graph.edges[j, k] = 1
                        else:
                            if self.verbose:
                                print(str(j) + "-" + str(k) + "0-*" + str(i), end=" ")
                                print("=> orient " + str(j) + "- " + str(k) + " -* " + str(i))
                            self.graph.edges[i, k] = 1
        return test_find_orientation

    def rule_7(self):
        """
        rule 7 from FCI
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule 7")
            print("######################################")

        test_find_orientation = False

        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 0) and (self.graph.edges[j, i] == 0):
                    k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and
                              (((self.graph.edges[j, k] == 3) and (self.graph.edges[k, j] != 0) and
                                (self.graph.edges[k, i] == 1) and (self.graph.edges[i, k] == 3)) or
                               ((self.graph.edges[i, k] == 3) and (self.graph.edges[k, i] != 0) and
                                (self.graph.edges[k, j] == 1) and (self.graph.edges[j, k] == 3)))]
                    if len(k_list) > 0:
                        test_find_orientation = True
                        for k in k_list:
                            if self.graph.edges[k, i] == 1:
                                if self.verbose:
                                    print(str(i) + "-0" + str(k) + "0-*" + str(j), end=" ")
                                    print("=> orient " + str(i) + "-0 " + str(k) + " -* " + str(j))
                                self.graph.edges[j, k] = 1
                            else:
                                if self.verbose:
                                    print(str(j) + "-0" + str(k) + "0-*" + str(i), end=" ")
                                    print("=> orient " + str(j) + "-0 " + str(k) + " -* " + str(i))
                                self.graph.edges[i, k] = 1
        return test_find_orientation

    def _is_potentially_directed(self, path):
        """
        check if path is a potentially directed path
        :param path: any path in the graph
        :return: bool
        """
        test_list1 = []
        for p in range(len(path)-1):
            test_list1.append((self.graph.edges[path[p+1], path[p]] != 2))
        return all(test_list1)

    def rule_8(self):
        """
        rule 8 from FCI
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule 8")
            print("######################################")

        test_find_orientation = False

        for i in range(self.graph.d):
            for j in range(self.graph.d):
                if (self.graph.edges[i, j] == 2) and (self.graph.edges[j, i] == 3):
                    k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and
                              (((self.graph.edges[i, k] == 2) and (self.graph.edges[k, i] == 1) and
                                (self.graph.edges[k, j] == 2) and (self.graph.edges[j, k] == 1)) or
                               ((self.graph.edges[i, k] == 3) and (self.graph.edges[k, i] == 1) and
                                (self.graph.edges[k, j] == 2) and (self.graph.edges[j, k] == 1)))]
                    if len(k_list) > 0:
                        test_find_orientation = True
                        for k in k_list:
                            if self.verbose:
                                if self.graph.edges[i, k] == 3:
                                    print(str(i) + "-0" + str(k) + "->" + str(j) + " and "+str(i) + "0->" + str(j),
                                          end=" ")
                                    print("=> orient " + str(i) + " -> " + str(j))
                                else:
                                    print(str(i) + "->" + str(k) + "->" + str(j) + " and "+str(i) + "0->" + str(j),
                                          end=" ")
                                    print("=> orient " + str(i) + " -> " + str(j))
                            self.graph.edges[j, i] = 1
        return test_find_orientation

    def rule_9(self):
        """
        rule 9 from FCI
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule 9")
            print("######################################")
        test_find_orientation = False
        for i in range(self.graph.d):
            for j in range(self.graph.d):
                if (self.graph.edges[i, j] == 2) and (self.graph.edges[j, i] == 3):
                    uncovered_path_list = self._find_uncovered_path(i, j)
                    if len(uncovered_path_list) > 0:
                        for p_d in uncovered_path_list:
                            if self._is_potentially_directed(p_d):
                                if self.graph.edges[p_d[-1], p_d[1]] == 0:
                                    test_find_orientation = True
                                    if self.verbose:
                                        print(str(i) + "0->" + str(j) + " and found a potential directed path", end=" ")
                                        print("=> orient " + str(i) + "->" + str(j))
                                    self.graph.edges[j, i] = 1
        return test_find_orientation

    def rule_10(self):
        """
        rule 10 from FCI
        :return: (bool) True if the rule made a change in the graph and False otherwise
        """
        if self.verbose:
            print("######################################")
            print("Rule 10")
            print("######################################")
        test_find_orientation = False
        for i in range(self.graph.d):
            for j in range(self.graph.d):
                if (self.graph.edges[i, j] == 2) and (self.graph.edges[j, i] == 3):
                    colliders_tails = [k for k in range(self.graph.d) if (k != j) and (
                            (self.graph.edges[k, j] == 2) and (self.graph.edges[j, k] == 1))]
                    colliders_tails = [list(k) for k in itertools.combinations(colliders_tails, 2)]
                    for ks in colliders_tails:
                        beta = ks[0]
                        theta = ks[1]
                        uncovered_path_list1 = self._find_uncovered_path(i, beta)
                        uncovered_path_list2 = self._find_uncovered_path(i, theta)
                        if (len(uncovered_path_list1) > 0) and (len(uncovered_path_list2) > 0):
                            for p1 in uncovered_path_list1:
                                if self._is_potentially_directed(p1):
                                    for p2 in uncovered_path_list2:
                                        if self._is_potentially_directed(p2):
                                            mu = p1[1]
                                            w = p2[1]
                                            if (mu != w) and (self.graph.edges[mu, w] == 0):
                                                test_find_orientation = True
                                                if self.verbose:
                                                    print(str(i) + "0->" + str(j) + " and ...", end=" ")
                                                    print("=> orient " + str(i) + "->" + str(j))
                                                self.graph.edges[j, i] = 1
        return test_find_orientation

    def rule_commun_confounder_and_causal_chain(self):
        """
        new rules (rule commun confounder (4) and rule causal_chain (5) from paper)
        """
        if self.verbose:
            print("######################################")
            print("Rule commun confounder and causal chain")
            print("######################################")

        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 0) and (self.graph.edges[j, i] == 0):
                    k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and (
                        ((self.graph.edges[j, k] == 3) and (self.graph.edges[k, j] == 3) and
                            (self.graph.edges[i, k] == 3) and (self.graph.edges[k, i] == 3)))]
                    if len(k_list) > 0:
                        for k in k_list:
                            gki = self.gamma_matrix[self.names[i]].loc[self.names[k]]
                            gkj = self.gamma_matrix[self.names[j]].loc[self.names[k]]
                            i_is_not_effet = (sum(self.graph.edges[:, i] == 2) == 0)
                            j_is_not_effet = (sum(self.graph.edges[:, j] == 2) == 0)
                            k_is_not_effet = (sum(self.graph.edges[:, k] == 2) == 0)

                            #Lagged common cause
                            if (gki > 0) and (gkj > 0):
                                if i_is_not_effet and j_is_not_effet:
                                    if self.verbose:
                                        print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)>0",
                                              end=" ")
                                        print("=> orient " + str(i) + "<- " + str(k) + " -> " + str(j))
                                    self.graph.edges[k, i] = 2
                                    self.graph.edges[k, j] = 2
                            #Lagged instantaneous confounder
                            elif (gki > 0) and (gkj == 0):
                                if i_is_not_effet:
                                    if j_is_not_effet and k_is_not_effet:
                                        if self.verbose:
                                            print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)==0",
                                                  end=" ")
                                            print("=> orient " + str(i) + "<- " + str(k) + " - " + str(j))
                                        self.graph.edges[k, i] = 2
                                    elif j_is_not_effet:
                                        if self.verbose:
                                            print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)==0",
                                                  end=" ")
                                            print("=> orient " + str(i) + "<- " + str(k) + " -> " + str(j))
                                        self.graph.edges[k, i] = 2
                                        self.graph.edges[k, j] = 2
                                    elif k_is_not_effet:
                                        if self.verbose:
                                            print(
                                                str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)==0",
                                                end=" ")
                                            print("=> orient " + str(i) + "<- " + str(k) + " <- " + str(j))
                                        self.graph.edges[k, i] = 2
                                        self.graph.edges[j, k] = 2
                            elif (gki == 0) and (gkj > 0):
                                if j_is_not_effet:
                                    if i_is_not_effet and k_is_not_effet:
                                        if self.verbose:
                                            print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)==0 and gamma(k,j)>0",
                                                  end=" ")
                                            print("=> orient " + str(i) + "- " + str(k) + " -> " + str(j))
                                        self.graph.edges[k, j] = 2
                                    elif i_is_not_effet:
                                        if self.verbose:
                                            print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)==0",
                                                  end=" ")
                                            print("=> orient " + str(i) + "<- " + str(k) + " -> " + str(j))
                                        self.graph.edges[k, i] = 2
                                        self.graph.edges[k, j] = 2
                                    elif k_is_not_effet:
                                        if self.verbose:
                                            print(
                                                str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,j)>0 and gamma(k,i)==0",
                                                end=" ")
                                            print("=> orient " + str(i) + "-> " + str(k) + " -> " + str(j))
                                        self.graph.edges[k, j] = 2
                                        self.graph.edges[i, k] = 2
                            # lagged instanteneous causal chain
                            elif (gki >= 0) and (gkj < 0):
                                if j_is_not_effet and k_is_not_effet:
                                    if self.verbose:
                                        print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)>0",
                                              end=" ")
                                        print("=> orient " + str(i) + "<- " + str(k) + " <- " + str(j))
                                    self.graph.edges[k, i] = 2
                                    self.graph.edges[j, k] = 2
                            elif (gki < 0) and (gkj >= 0):
                                if i_is_not_effet and k_is_not_effet:
                                    if self.verbose:
                                        print(str(i) + "-" + str(k) + "-" + str(j) + "and gamma(k,i)>0 and gamma(k,j)>0",
                                              end=" ")
                                        print("=> orient " + str(i) + "-> " + str(k) + " -> " + str(j))
                                    self.graph.edges[i, k] = 2
                                    self.graph.edges[k, j] = 2

    def rule_mediator(self):
        """
        new rules (rule mediator (6) from paper)
        """
        if self.verbose:
            print("######################################")
            print("Rule mediator")
            print("######################################")

        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] != 0) and (self.graph.edges[j, i] != 0):
                    k_list = [k for k in range(self.graph.d) if (k != i) and (k != j) and (
                        ((self.graph.edges[j, k] == 3) and (self.graph.edges[k, j] == 3) and
                            (self.graph.edges[i, k] == 3) and (self.graph.edges[k, i] == 3)))]
                    if len(k_list) > 0:
                        for k in k_list:
                            gij = self.gamma_matrix[self.names[j]].loc[self.names[i]]
                            gik = self.gamma_matrix[self.names[k]].loc[self.names[i]]
                            gjk = self.gamma_matrix[self.names[k]].loc[self.names[j]]
                            # g_list = [(gij, gik, gjk), (gij, gik, -gjk), (-gij, gjk, gik), (-gij, gjk, -gik),
                            #           (-gik, -gjk, gij), (-gik, -gjk, -gij)]
                            # i->j->k, i->k->j, j->i->k, j->k->->i, k->i->j, k->j->i
                            g_list = [(gij, gjk, gik), (gik, -gjk, gij), (-gij, gik, gjk), (gjk, -gik, -gij),
                                      (-gik, gij, -gjk), (-gjk, -gij, -gik)]
                            g_list_common = [(gij, gjk, gik), (gik, -gjk, gij), (gjk, -gik, -gij)]

                            msk = [(x[0] > 0) and (x[1] > 0) and (x[2] >= 0) for x in g_list]
                            msk_common = [(x[0] == 0) and (x[1] > 0) and (x[2] > 0) for x in g_list_common]
                            if any(msk):
                                print(g_list)
                                print(msk)
                                s = int(np.argwhere(msk)[0])
                                # g1, g2, g3 = g_list[s]
                                if s == 0:
                                    if (sum(self.graph.edges[:, j] == 2) == 0) and (
                                            sum(self.graph.edges[:, k] == 2) == 0):
                                        if (self.graph.edges[j, i] == 3) and (self.graph.edges[k, i] == 3) \
                                                and (self.graph.edges[k, j] == 3):
                                            if self.verbose:
                                                print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "-> " + str(j) + " -> " + str(k) + " <- "
                                                      + str(i))
                                            self.graph.edges[i, j] = 2
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[j, k] = 2
                                elif s == 1:
                                    if (sum(self.graph.edges[:, j] == 2) == 0) and (
                                            sum(self.graph.edges[:, k] == 2) == 0):
                                        if (self.graph.edges[j, i] == 3) and (self.graph.edges[k, i] == 3) \
                                                and (self.graph.edges[k, j] == 3):
                                            if self.verbose:
                                                print(str(i) + "-" + str(k) + "-" + str(j) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "-> " + str(k) + "-> " + str(j) + " <- "
                                                      + str(i))
                                            self.graph.edges[i, j] = 2
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[k, j] = 2
                                elif s == 2:
                                    if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                            sum(self.graph.edges[:, k] == 2) == 0):
                                        if (self.graph.edges[i, j] == 3) and (self.graph.edges[k, i] == 3) \
                                                and (self.graph.edges[k, j] == 3):
                                            if self.verbose:
                                                print(str(j) + "-" + str(i) + "-" + str(k) + "-" + str(j), end=" ")
                                                print("=> orient " + str(j) + "-> " + str(i) + " -> " + str(k) + " <- "
                                                      + str(j))
                                            self.graph.edges[j, i] = 2
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[j, k] = 2
                                if s == 3:
                                    if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                            sum(self.graph.edges[:, k] == 2) == 0):
                                        if (self.graph.edges[i, j] == 3) and (self.graph.edges[i, k] == 3) \
                                                and (self.graph.edges[k, j] == 3):
                                            if self.verbose:
                                                print(str(j) + "-" + str(k) + "-" + str(i) + "-" + str(j), end=" ")
                                                print("=> orient " + str(j) + "-> " + str(k) + "-> " + str(i) + " <- "
                                                      + str(j))
                                            self.graph.edges[j, i] = 2
                                            self.graph.edges[k, i] = 2
                                            self.graph.edges[j, k] = 2
                                elif s == 4:
                                    if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                            sum(self.graph.edges[:, j] == 2) == 0):
                                        if (self.graph.edges[j, i] == 3) and (self.graph.edges[i, k] == 3) \
                                                and (self.graph.edges[j, k] == 3):
                                            if self.verbose:
                                                print(str(k) + "-" + str(i) + "-" + str(j) + "-" + str(k), end=" ")
                                                print("=> orient " + str(k) + "-> " + str(i) + "-> " + str(j) + " <- "
                                                      + str(k))
                                            self.graph.edges[i, j] = 2
                                            self.graph.edges[k, i] = 2
                                            self.graph.edges[k, j] = 2
                                elif s == 5:
                                    if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                            sum(self.graph.edges[:, j] == 2) == 0):
                                        if (self.graph.edges[i, j] == 3) and (self.graph.edges[i, k] == 3) \
                                                and (self.graph.edges[j, k] == 3):
                                            if self.verbose:
                                                print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                      + str(k))
                                            self.graph.edges[j, i] = 2
                                            self.graph.edges[k, i] = 2
                                            self.graph.edges[k, j] = 2
                            elif any(msk_common):
                                s = int(np.argwhere(msk_common)[0])
                                if s == 0:
                                    if (self.graph.edges[j, i] == 3) and (self.graph.edges[k, i] == 3) \
                                            and (self.graph.edges[k, j] == 3):
                                        if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                sum(self.graph.edges[:, j] == 2) == 0) and (
                                                sum(self.graph.edges[:, k] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "- " + str(j) + " -> " + str(k) + " <- "
                                                      + str(i))
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[j, k] = 2
                                        elif (sum(self.graph.edges[:, j] == 2) == 0) and (
                                                sum(self.graph.edges[:, k] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "-> " + str(j) + " -> " + str(k) + " <- "
                                                      + str(i))
                                            self.graph.edges[i, j] = 2
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[j, k] = 2
                                        elif (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                sum(self.graph.edges[:, k] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
                                                print("=> orient " + str(j) + "-> " + str(i) + " -> " + str(k) + " <- "
                                                      + str(j))
                                            self.graph.edges[j, i] = 2
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[j, k] = 2
                                elif s == 1:
                                    if (self.graph.edges[j, i] == 3) and (self.graph.edges[k, i] == 3) \
                                            and (self.graph.edges[k, j] == 3):
                                        if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                sum(self.graph.edges[:, j] == 2) == 0) and (
                                                sum(self.graph.edges[:, k] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(k) + "-" + str(j) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "- " + str(k) + " -> " + str(j) + " <- "
                                                      + str(i))
                                            self.graph.edges[i, j] = 2
                                            self.graph.edges[k, j] = 2
                                        elif (sum(self.graph.edges[:, k] == 2) == 0) and (
                                                sum(self.graph.edges[:, j] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "-> " + str(j) + " -> " + str(k) + " <- "
                                                      + str(i))
                                            self.graph.edges[i, j] = 2
                                            self.graph.edges[i, k] = 2
                                            self.graph.edges[k, j] = 2
                                        elif (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                sum(self.graph.edges[:, j] == 2) == 0):
                                            if self.verbose:
                                                print(str(k) + "-" + str(i) + "-" + str(j) + "-" + str(k), end=" ")
                                                print("=> orient " + str(k) + "-> " + str(i) + " -> " + str(j) + " <- "
                                                      + str(k))
                                            self.graph.edges[k, i] = 2
                                            self.graph.edges[k, j] = 2
                                            self.graph.edges[i, j] = 2
                                elif s == 2:
                                    if (self.graph.edges[j, i] == 3) and (self.graph.edges[k, i] == 3) \
                                            and (self.graph.edges[k, j] == 3):
                                        if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                sum(self.graph.edges[:, j] == 2) == 0) and (
                                                sum(self.graph.edges[:, k] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(k) + "-" + str(j) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "- " + str(k) + " -> " + str(j) + " <- "
                                                      + str(i))
                                            self.graph.edges[j, i] = 2
                                            self.graph.edges[k, i] = 2
                                        elif (sum(self.graph.edges[:, k] == 2) == 0) and (
                                                sum(self.graph.edges[:, i] == 2) == 0):
                                            if self.verbose:
                                                print(str(i) + "-" + str(j) + "-" + str(k) + "-" + str(i), end=" ")
                                                print("=> orient " + str(i) + "-> " + str(j) + " -> " + str(k) + " <- "
                                                      + str(i))
                                            self.graph.edges[j, k] = 2
                                            self.graph.edges[j, i] = 2
                                            self.graph.edges[k, i] = 2
                                        elif (sum(self.graph.edges[:, j] == 2) == 0) and (
                                                sum(self.graph.edges[:, i] == 2) == 0):
                                            if self.verbose:
                                                print(str(k) + "-" + str(i) + "-" + str(j) + "-" + str(k), end=" ")
                                                print("=> orient " + str(k) + "-> " + str(i) + " -> " + str(j) + " <- "
                                                      + str(k))
                                            self.graph.edges[k, j] = 2
                                            self.graph.edges[k, i] = 2
                                            self.graph.edges[j, i] = 2
                            else:
                                if (self.graph.edges[i, j] == 3) and (self.graph.edges[i, k] == 3) \
                                        and (self.graph.edges[k, j] == 3):
                                    if (gij!=0) and (gik==0) and (gjk==0):
                                        if gij>0:
                                            if (sum(self.graph.edges[:, j] == 2) == 0) and (
                                                    sum(self.graph.edges[:, k] == 2) == 0):
                                                if self.verbose:
                                                    print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                    print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                          + str(k))
                                                self.graph.edges[i, j] = 2
                                                self.graph.edges[i, k] = 2
                                                self.graph.edges[j, k] = 2
                                        elif gij<0:
                                            if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                    sum(self.graph.edges[:, k] == 2) == 0):
                                                if self.verbose:
                                                    print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                    print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                          + str(k))
                                                self.graph.edges[j, i] = 2
                                                self.graph.edges[i, k] = 2
                                                self.graph.edges[j, k] = 2
                                    elif (gij==0) and (gik!=0) and (gjk==0):
                                        if gik>0:
                                            if (sum(self.graph.edges[:, k] == 2) == 0) and (
                                                    sum(self.graph.edges[:, j] == 2) == 0):
                                                if self.verbose:
                                                    print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                    print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                          + str(k))
                                                self.graph.edges[i, k] = 2
                                                self.graph.edges[i, j] = 2
                                                self.graph.edges[k, j] = 2
                                        if gik<0:
                                            if (sum(self.graph.edges[:, i] == 2) == 0) and (
                                                    sum(self.graph.edges[:, j] == 2) == 0):
                                                if self.verbose:
                                                    print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                    print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                          + str(k))
                                                self.graph.edges[k, i] = 2
                                                self.graph.edges[i, j] = 2
                                                self.graph.edges[k, j] = 2
                                    elif (gij == 0) and (gik == 0) and (gjk != 0):
                                        if gjk>0:
                                            if (sum(self.graph.edges[:, k] == 2) == 0) and (
                                                    sum(self.graph.edges[:, i] == 2) == 0):
                                                if self.verbose:
                                                    print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                    print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                          + str(k))
                                                self.graph.edges[j, k] = 2
                                                self.graph.edges[j, i] = 2
                                                self.graph.edges[k, i] = 2
                                        if gjk<0:
                                            if (sum(self.graph.edges[:, j] == 2) == 0) and (
                                                    sum(self.graph.edges[:, i] == 2) == 0):
                                                if self.verbose:
                                                    print(str(k) + "-" + str(j) + "-" + str(i) + "-" + str(k), end=" ")
                                                    print("=> orient " + str(k) + "-> " + str(j) + " -> " + str(i) + " <- "
                                                          + str(k))
                                                self.graph.edges[k, j] = 2
                                                self.graph.edges[j, i] = 2
                                                self.graph.edges[k, i] = 2

    def rule_proba_raising_principle(self):
        """
        new rules (rule prob raising principle from paper)
        """
        if self.verbose:
            print("######################################")
            print("Rule prob raising principle")
            print("######################################")
        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 3) and (self.graph.edges[j, i] == 3):

                    adjacent_i_is_1 = (sum(np.delete(self.graph.edges[:, i], i) != 0) == 1)
                    adjacent_j_is_1 = (sum(np.delete(self.graph.edges[:, j], j) != 0) == 1)
                    if adjacent_i_is_1 and adjacent_j_is_1:
                        gij = self.gamma_matrix[self.names[j]].loc[self.names[i]]
                        if gij > 0:
                            if self.verbose:
                                print(str(i) + "-" + str(j) + "g(i,j)>0", end=" ")
                                print("=> orient " + str(i) + "-> " + str(j))
                            self.graph.edges[i, j] = 2
                        elif gij < 0:
                            if self.verbose:
                                print(str(i) + "-" + str(j) + "g(i,j)<0", end=" ")
                                print("=> orient " + str(i) + "<- " + str(j))
                            self.graph.edges[j, i] = 2

    def rule_gap_orientation(self):
        """
        gamma heuristic rule from paper
        """
        if self.verbose:
            print("######################################")
            print("Rule gap orientation")
            print("######################################")

        for i in range(self.graph.d):
            for j in range(i + 1, self.graph.d):
                if (self.graph.edges[i, j] == 3) and (self.graph.edges[j, i] == 3):
                    if self.gamma_matrix[self.names[j]].loc[self.names[i]] > 0:
                        if self.verbose:
                            print(str(i) + "-" + str(j) + "g(i,j)>0", end=" ")
                            print("=> orient " + str(i) + "-> " + str(j))
                        self.graph.edges[i, j] = 2
                    if self.gamma_matrix[self.names[j]].loc[self.names[i]] < 0:
                        if self.verbose:
                            print(str(i) + "-" + str(j) + "g(i,j)<0", end=" ")
                            print("=> orient " + str(i) + "<- " + str(j))
                        self.graph.edges[j, i] = 2

    def fit(self):
        """
        run FCITMI
        :return: graph (PAG)
        """
        if self.verbose:
            print("#######################################")
            print("########### Starting FCITMI ###########")
            print("#######################################")

        # initialize skeleton
        self.skeleton_initialize()

        # get separation sets
        self.find_sep_set()

        # include circle in the skeleton
        self.dag_to_pag()

        # orientation
        self.rule_temporal_priority_within_time_series()
        self.rule_origin_causality()


        # find possible d sep
        # test_remove_links = self.find_d_sep()
        test_remove_links = False

        # remove orientation
        if test_remove_links:
            self.remove_orientation()
            # orientation
            self.rule_origin_causality()
        test_rp = True
        test_r2 = True
        test_r3 = True
        test_r4 = True
        while test_rp or test_r2 or test_r3 or test_r4:
            test_rp = self.rule_propagation_causality()
            test_r2 = self.rule_2()
            test_r3 = self.rule_3()
            test_r4 = self.rule_4()

        self.rule_commun_confounder_and_causal_chain()
        self.rule_mediator()
        self.rule_proba_raising_principle()

        test_r8 = True
        test_r9 = True
        test_r10 = True
        while test_r8 or test_r9 or test_r10:
            test_r8 = self.rule_8()
            test_r9 = self.rule_9()
            test_r10 = self.rule_10()

        if self.verbose:
            print("######################################")
            print("Final Results (FCITMI)")
            print("######################################")
            print("Summary Graph:")
            print(self.graph.edges)
        return self.graph.edges

    def fit_gap_orientation(self):
        """
        run FCITMI-gamma (requirements: run FCITMI)
        :return: graph (PAG)
        """
        self.rule_gap_orientation()
        return self.graph.edges


if __name__ == "__main__":
    import pandas as pd
    path = "../../../../data/simulated_ts_data/v_structure/data_"+str(0)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)
    data = data.loc[:1000]
    print(data)
    import random
    random.seed(2)

    ci = TPCTMI(data, verbose=True, p_value=False, num_processor=1)

    gamma_matrix = ci.gamma_matrix
    print(gamma_matrix)
    print(data.shape)

    x = data[data.columns[0]].loc[1:].reset_index(drop=True)
    y = data[data.columns[2]].loc[1:].reset_index(drop=True)
    z = data[data.columns[0]].loc[:999].reset_index(drop=True)
    z.name = "V1_past"
    z = {'V1_past': z}
    print(x.shape, y.shape)
    sampling_rate = {'V1': 1, 'V2': 1, 'V3': 1, 'V1_past': 1, 'V2_past':1, 'V3_past':1}
    cmi, _ = ctmi(x, y, z, data.columns[0], data.columns[2], sampling_rate,
                  gamma_matrix=gamma_matrix, p_value=True)
    print(cmi)
    print("done")

    # ci = FCITMI(data, verbose=True, p_value=True, num_processor=1)
    ci = TPCTMI(data, verbose=True, p_value=False, num_processor=1)
    ci.fit()

    print(ci.tgraph_dict)

    ci.noise_based_hidden_counfounders()

    # ar = np.zeros([4, 4])
    # ar[0,1] = 2
    # ar[1,0] = 1
    #
    # ar[1,2] = 2
    # ar[2,1] = 2
    # ar[1,3] = 2
    # ar[3,1] = 1
    #
    # ar[2,3] = 2
    # ar[3,2] = 2
    # ci.graph.edges = ar
    # print(ar)
    # print(ci._find_discriminating_paths(0, 3))