import random
import torch
import torch.nn as nn
from torch import optim

import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
import time

from tools.tigramite.tigramite.independence_tests import CMIknn, ParCorr

import itertools
from joblib import Parallel, delayed

# from ctmi import window_representation, get_sampling_rate, align_matrix, tmi, get_alpha, window_size, align_pair
# from ctmi_new import i_ctmi, ctmi
# from gctmi import gctmi

from ctmi import window_representation, get_sampling_rate, align_matrix, tmi, get_alpha
from ctmi_new import ctmi, align_matrix, tmi, gamma_matrix_window_matrix
# from gctmi import gctmi

from datetime import datetime

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class TestMI:
    def __init__(self, p_value= True):
        self.cd = CMIknn(mask_type=None, significance='shuffle_test', fixed_thres=None, sig_samples=10000,
                sig_blocklength=3, knn=10, confidence='bootstrap', conf_lev=0.9, conf_samples=10000,
                conf_blocklength=1, verbosity=0)
        self.p_value = p_value

    def fit(self, x, y, z=None):
        if len(x.shape) == 1:
            x = x.reshape(-1, 1)
        if len(y.shape) == 1:
            y = y.reshape(-1, 1)
        dim_x = x.shape[1]
        dim_y = y.shape[1]

        ws_xy = 1
        # y_past = y[:-ws_xy]
        # x = x[ws_xy:]#.reset_index(drop=True)
        # y = y[ws_xy:]#.reset_index(drop=True)

        if z is not None:
            # z = z[ws_xy:]  # .reset_index(drop=True)
            # z = np.concatenate((z, y_past), axis=1)

            dim_z = z.shape[1]
            X = np.concatenate((x, y, z), axis=1)
            xyz = np.array([0] * dim_x + [1] * dim_y+ [2] * dim_z)
        else:
            # X = np.concatenate((x, y, y_past), axis=1)
            X = np.concatenate((x, y), axis=1)
            # xyz = np.array([0] * dim_x + [1] * dim_y + [2] * ws_xy)
            xyz = np.array([0] * dim_x + [1] * dim_y)
        value = self.cd.get_dependence_measure(X.T, xyz)
        if self.p_value:
            pvalue = self.cd.get_shuffle_significance(X.T, xyz, value)
            return pvalue, value
        else:
            return 0, value


class TestParCorr:
    def __init__(self):
        self.cd = ParCorr(mask_type=None, significance='shuffle_test', fixed_thres=None, sig_samples=10000,
                sig_blocklength=3, confidence='bootstrap', conf_lev=0.9, conf_samples=10000,
                conf_blocklength=1, verbosity=0)

    def fit(self, x, y, z=None):
        if len(x.shape) == 1:
            x = x.reshape(-1, 1)
        if len(y.shape) == 1:
            y = y.reshape(-1, 1)
        dim_x = x.shape[1]
        dim_y = y.shape[1]
        if z is not None:
            dim_z = z.shape[1]
            X = np.concatenate((x, y, z), axis=1)
            xyz = np.array([0] * dim_x + [1] * dim_y+ [2] * dim_z)
        else:
            X = np.concatenate((x, y), axis=1)
            xyz = np.array([0] * dim_x + [1] * dim_y)
        value = self.cd.get_dependence_measure(X.T, xyz)
        # pvalue = self.cd.get_shuffle_significance(X.T, xyz, value)
        pvalue = 0
        return pvalue, value


class tsCNN(nn.Module):
    def __init__(self, input_size, output_size, input_lag):
        super(tsCNN, self).__init__()
        self.input_size = input_size
        self.output_size = output_size
        self.input_lag = input_lag
        self.compact_ts = nn.Linear(input_lag, 1)
        self.compact_in = nn.Linear(int(input_size/input_lag), 1)
        self.conv1 = nn.Sequential(
            nn.Conv1d(
                in_channels=1,
                out_channels=input_size,  # Some random number
                kernel_size=5,
                stride=1,
                padding=2,
            ),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2),  # size after pooling?
        )
        # self.conv2 = nn.Sequential(
        #     nn.Conv1d(
        #         in_channels=16,
        #         out_channels=8,
        #         kernel_size=5,
        #         stride=1,
        #         padding=2,
        #     ),
        #     nn.ReLU(),
        #     nn.MaxPool1d(kernel_size=2),
        #
        # )
        # self.compact_out = nn.Linear(8, 1)
        self.out = nn.Linear(input_size, output_size)

    def forward(self, x_dict):
        # print(x.size())
        compact_ts_dict = dict()
        names = list(x_dict.keys())
        for name in names:
            x = x_dict[name].view(-1, self.input_lag)
            compact_ts_i = self.compact_ts(x)
            compact_ts_dict[name] = compact_ts_i
            if name == names[0]:
                compact_ts = compact_ts_i
            else:
                compact_ts = torch.cat((compact_ts, compact_ts_i), 1)
        compact_in = self.compact_in(compact_ts)
        x = compact_in.view(-1, 1, 1)
        x = self.conv1(x)
        # print(x.size())
        # x = self.conv2(x)
        x = x.view(-1, self.input_size)
        # compact_out = self.compact_out(x)
        # output = self.out(compact_out)
        output = self.out(x)
        return output, compact_in, compact_ts_dict


def train(input_dict, target_tensor, model, optimizer, criterion):

    optimizer.zero_grad()
    # model.zero_grad()

    output, _, _ = model(input_dict)
    loss = criterion(output, target_tensor)

    loss.backward()
    optimizer.step()
    return loss.item()


def predict(input_dict, model):
    output, compact_out, compact_ts_dict = model(input_dict)
    return output, compact_out, compact_ts_dict


# Function to produce noise
def add_noise(x, d, order, beta=0.5):
    x = x.copy()
    rand = np.random.randint(0, high=d, size=x.shape[0])
    for j in range(d):
        proba = np.random.random(size=1)
        if proba > beta:
            for o in range(order-1):
                i = j + o*d
                x[i, rand[j]] = 0
    return x


def mts_order(mts, order=4):
    new_mts = pd.DataFrame()
    for i in range(order):
        if i == order:
            i_data = mts[i:]
        else:
            i_data = mts[i:(-order + i)]
        if isinstance(mts, pd.DataFrame):
            names_col = mts.columns.values+ "_" + str(i + 1)
        elif isinstance(mts, pd.Series):
            names_col = mts.name + "_" + str(i + 1)
        else:
            print('error!')
            exit(0)
        for j in range(len(names_col)):
            new_mts[names_col[j]] = i_data[mts.columns.values[j]].values
    return new_mts


def tskiko_mv(data, max_lag, learning_rate, training_epoch, noise=True, alpha=0.05, cond_ind_test="ParCorr",
           verbose=True):
    """
    :param data: input
    :param max_lag: max_lag
    :param learning_rate: learning rate of the autoencoder
    :param training_epoch: number of training epochs
    :param num_neurons: number of neurones in the hidden layer
    :param noise: boolean value, if true a denoising autoencoder should be use
    :param alpha:
    :param cond_ind_test: CMI or ParCorr
    :param verbose:
    :return: dict
    """
    # cond_ind_test = "CMI"
    option = 1
    # Start Causal ordering
    start = time.time()

    # scaler = MinMaxScaler(feature_range=(-1, 1))
    # data = pd.DataFrame(scaler.fit_transform(data.values), columns=data.columns)
    # data.columns = data.columns.values.astype(str)
    d = data.shape[1]

    x = mts_order(data, order=max_lag)  # [:-order]
    names_x = x.columns[:-d]
    names_y = x.columns[-d:]
    y = x[names_y]
    x = x[names_x]

    summary_names = list(data.columns)
    temporal_names = dict()
    for s in range(d):
        temporal_names[summary_names[s]] = []
        for o in range(max_lag - 1):
            i = s + o * d
            temporal_names[summary_names[s]].append(names_x[i])

    cost_history = []
    indep_history = []
    test_indep_history = []

    x_train = x.copy()

    S = list(data.columns)
    sig = []
    pa = dict()
    for name in summary_names:
        pa[name] = []

    # todo compare each time series with its past
    for j in range(d-1):
        # if j != 0:
        #     x_train.drop(temporal_names[selected], axis=1, inplace=True)
        #     y.drop(y.columns[selected_loc], axis=1, inplace=True)
        #     del S[selected_loc]
        criterion = nn.MSELoss()
        model = tsCNN(x_train.shape[1], y.shape[1], max_lag-1).to(device)
        optimizer = optim.Adam(model.parameters(), lr=learning_rate)

        n_epochs_stop = 100
        epochs_no_improve = 0
        min_loss = np.inf
        for iter in range(training_epoch + 1):
            # mini_batch_size = 5
            # N = x_train.shape[0] - 2
            # n_batch = N // mini_batch_size + (N % mini_batch_size != 0)
            # i_batch = (iter % N)
            if noise:
                x_train_n = add_noise(x_train.values, d=len(S), order=max_lag)
            else:
                x_train_n = x_train.values.copy()
            x_train_n = pd.DataFrame(x_train_n, columns=x_train.columns)

            input_dict = dict()
            for i in range(int(x_train_n.shape[1]/(max_lag-1))):
                input_tensor = torch.tensor(
                    x_train_n[temporal_names[S[i]]].values.reshape(-1, max_lag-1, 1), dtype=torch.float,
                    device=device)
                input_dict[S[i]] = input_tensor
            # input_tensor = torch.tensor(
            #     x_train_n.reshape(-1, x_train_n.shape[1], 1), dtype=torch.float,
            #     device=device)
            target_tensor = torch.tensor(
                y.values.reshape(-1, y.shape[1]), dtype=torch.float,
                device=device)

            loss = train(input_dict, target_tensor, model, optimizer, criterion)
            if loss < min_loss:
                min_loss = loss
                epochs_no_improve = 0
            else:
                epochs_no_improve = epochs_no_improve + 1
            if iter > 100 and epochs_no_improve == n_epochs_stop:
                if verbose:
                    print('Early stopping!')
                    print("Epoch ", iter, "MSE: ", "{:.4f}".format(loss, 4))
                break

            if verbose:
                if iter % 250 == 0:
                    print("Epoch ", iter, "MSE: ", "{:.4f}".format(loss, 4))

            # loss = train(input_dict, target_tensor, model, optimizer, criterion)
            # if verbose:
            #     if iter % 100 == 0:
            #         print("Epoch ", iter, "MSE: ", "{:.4f}".format(loss, 4))
        cost_history.append(loss)

        # k = d - 1 - j
        test_indep_values = []
        indep_values = []
        c = x_train.copy()
        input_dict = dict()
        for i in range(int(c.shape[1] / (max_lag - 1))):
            input_tensor = torch.tensor(
                c[temporal_names[S[i]]].values.reshape(-1, max_lag - 1, 1), dtype=torch.float,
                device=device)
            input_dict[S[i]] = input_tensor
        # input_tensor = torch.tensor(c.values.reshape(c.shape[0], c.shape[1], 1), dtype=torch.float, device=device)
        res, compact_res, compact_ts_dict = predict(input_dict, model)
        res = res.detach().numpy()
        compact_res = compact_res.detach().numpy()
        for s in range(len(S)):
            # for o in range(order-1):
            #     i = s + o*len(S)
            #     c[:, i] = np.zeros((x_train.shape[0]))
            e = y.values[:, s] - res[:, s]
            # hs = TestHSIC(kernel='rbf')
            if cond_ind_test == "ParCorr":
                hs = TestParCorr()
            elif cond_ind_test == "CMI":
                hs = TestMI()
            # pval, val = hs.fit(compact_res, e, c[temporal_names[S[s]]])
            cond = compact_ts_dict[S[s]].detach().numpy()
            pval, val = hs.fit(compact_res, e, cond)
            # pval, val = hs.fit(compact_res, e)
            test_indep_values.append(pval)
            indep_values.append(abs(val))
            # indep_values.append(hs.fit(compact_res.detach().numpy(), e))
        indep_history.append(indep_values)
        test_indep_history.append(test_indep_values)
        # test_indep_array = np.array(test_indep_values).reshape(-1, len(S))
        test_indep_array = pd.DataFrame(np.array(test_indep_values).reshape(-1, len(S)), columns=S, index=[0])
        # indep_array = np.array(indep_values).reshape(-1, len(S))
        indep_array = pd.DataFrame(np.array(indep_values).reshape(-1, len(S)), columns=S, index=[0])
        if test_indep_values.count(test_indep_values[0]) == len(test_indep_values):
            selected = indep_array.idxmin(axis=1).loc[0]
            if verbose:
                print("since all p-values are the same, we are looking at the statistics...")
                print('indeps :' + str(indep_values))
        else:
            if verbose:
                print('test indeps :' + str(test_indep_values))
            selected = test_indep_array.idxmax(axis=1).loc[0]
        pval_init = test_indep_array[selected].loc[0]
        sig.insert(0, selected)

        pa[selected] = summary_names.copy()
        # pa[S[idp_init]].remove(S[idp_init])
        for name in sig:
            pa[selected].remove(name)
        selected_loc = test_indep_array.columns.get_loc(selected)

        c = x_train.copy()

        print("selected:" +str(selected))
        print("candidate parents" +str(pa[selected]))

        x_train.drop(temporal_names[selected], axis=1, inplace=True)
        y.drop(y.columns[selected_loc], axis=1, inplace=True)
        del S[selected_loc]

    if len(S) == 1:
        sig[0] = S[0]
    print(sig)

    end = time.time()
    discovery_time = end - start
    print('time causal discovery: '+str(discovery_time))

    print(pa)

    res_unit_array = pd.DataFrame(np.zeros([d, d]), columns=summary_names, index=summary_names, dtype=int)
    for k in pa.keys():
        res_unit_array[k].loc[k] = 1
        temp = pa[k]
        for i in temp:
            # if k == i:
            # res_unit_array[i].loc[i] = 1
            # else:
            if res_unit_array[i].loc[k] == 0:
                res_unit_array[i].loc[k] = 1
            res_unit_array[k].loc[i] = 2

    return res_unit_array


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


class KITMI:
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
        self.series = series
        self.graph = Graph(series.shape[1])

        training_epoch = 1000
        noise = True  # d*(order-1)*2
        learning_rate = 0.01
        for i in range(series.shape[1]):
            for j in range(i+1, series.shape[1]):
                data_pair = series[[series.columns[i], series.columns[j]]]
                res_order_pair = tskiko_mv(data_pair, lag_max, learning_rate, training_epoch, noise, sig_lev, "ParCorr", verbose)
                if res_order_pair[series.columns[j]].loc[series.columns[i]] == 2:
                    self.graph.edges[i, j] = 2
                if res_order_pair[series.columns[i]].loc[series.columns[j]] == 2:
                    self.graph.edges[j, i] = 2
        # self.graph.edges = tskiko_mv(self.series[[series.columns[0], series.columns[1]]], lag_max, learning_rate, training_epoch, noise, sig_lev, "ParCorr", verbose)

        if verbose:
            print("Order")
            print(self.graph.edges)


        self.adaptive_window = True
        self.series = series
        # self.graph = Graph(series.shape[1])
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
            self.gamma_matrix, self.window_matrix = self.gamma_matrix_window_matrix(self.series, series.columns)
        else:
            self.gamma_matrix = align_matrix(self.data_dict, series.columns, self.sampling_rate)

        self.cap_gamma_df = pd.DataFrame(columns=["p", "q", "r", "Grp", "Grq"])

        self.mi_array = np.ones([self.graph.d, self.graph.d])
        self.cmi_array = np.ones([self.graph.d, self.graph.d])

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
            print("Orderrrr")
            print(self.graph.edges)

    def find_gamma_lambda_x_y(self, x, y, k=10, max_gamma=5):
        gamma_list = list(range(1, max_gamma))
        # todo add windows
        # ws_x_list = list(range(1, max_gamma - 2))
        # ws_y_list = list(range(1, max_gamma - 2))
        ws_x_list = [1]
        ws_y_list = [1]

        c = np.zeros([len(gamma_list), len(ws_x_list), len(ws_y_list)])

        for idx_g in range(len(gamma_list)):
            for idx_ws_x in range(len(ws_x_list)):
                x_w_rep = window_representation(x, windows_size=ws_x_list[idx_ws_x])
                for idx_ws_y in range(len(ws_y_list)):
                    # if ws_x_list[idx_ws_x] == ws_y_list[idx_ws_y] == 1:
                    y_w_rep = window_representation(y, windows_size=ws_y_list[idx_ws_y])
                    g = gamma_list[idx_g]

                    if g > 0:
                        y_w_rep = y_w_rep[g:]
                        x_w_rep = x_w_rep.reset_index(drop=True)
                        y_w_rep = y_w_rep.reset_index(drop=True)

                        x_w_rep = x_w_rep[:-g]
                        x_w_rep = x_w_rep.reset_index(drop=True)
                        y_w_rep = y_w_rep.reset_index(drop=True)
                    m = min(x_w_rep.shape[0], y_w_rep.shape[0])
                    x_w_rep = x_w_rep[:m]
                    y_w_rep = y_w_rep[:m]
                    if len(x_w_rep.shape) == 1:
                        x_w_rep = x_w_rep.to_frame()
                    if len(y_w_rep.shape) == 1:
                        y_w_rep = y_w_rep.to_frame()
                    cmi = TestMI(p_value=False)
                    _, val = cmi.fit(x_w_rep, y_w_rep)

                    c[idx_g, idx_ws_x, idx_ws_y] = val
                    # else:
                    #     if ws_x_list[idx_ws_x] != ws_y_list[idx_ws_y]:
                    #         y_w_rep = window_representation(y, windows_size=ws_y_list[idx_ws_y])
                    #         g = gamma_list[idx_g]
                    #         _, val = tmi(x_w_rep, y_w_rep, sampling_rate_tuple, k=k, gamma=g, p_value=False)
                    #         c[idx_g, idx_ws_x, idx_ws_y] = val
                    #     else:
                    #         c[idx_g, idx_ws_x, idx_ws_y] = 0

        idx_g, idx_ws_x, idx_ws_y = np.where(c == np.max(c))
        idx_g = idx_g[0]
        idx_ws_x = idx_ws_x[0]
        idx_ws_y = idx_ws_y[0]
        g = gamma_list[idx_g]
        ws_x = ws_x_list[idx_ws_x]
        ws_y = ws_y_list[idx_ws_y]
        return g, ws_x, ws_y

    def gamma_matrix_window_matrix(self, series, keys, k=10, max_gamma=5):
        d = len(keys)
        g_matrix = np.zeros([d, d], dtype=int)
        window_matrix = np.zeros([d, d], dtype=list)

        for i in range(d):
            for j in range(d):
                if i != j:
                    x = series[keys[i]]
                    y = series[keys[j]]
                    g, ws_x, ws_y = self.find_gamma_lambda_x_y(x, y, k=k, max_gamma=max_gamma)
                    g_matrix[i, j] = g
                    window_matrix[i, j] = [ws_x, ws_y]
                    # window_matrix[j, i] = ws_y
                else:
                    g_matrix[i, j] = 1
                    window_matrix[i, j] = [1, 1]
                    # window_matrix[j, i] = 1
        return pd.DataFrame(g_matrix, columns=keys, index=keys), pd.DataFrame(window_matrix, columns=keys, index=keys)

    def align_pq(self, x, y, gamma):
        x = x.loc[y.index[0]:]

        idx_x = x.index
        idx_y = y.index
        if gamma > 0:
            y = y[gamma:]
            idx_y = idx_y[gamma:]
            x = x.reset_index(drop=True)
            y = y.reset_index(drop=True)
            idx_x = idx_x[x.index]
            idx_y = idx_y[y.index]

            x = x[:-gamma]
            idx_x = idx_x[:-gamma]
            x = x.reset_index(drop=True)
            y = y.reset_index(drop=True)
        else:
            print("Error: gamma <= 0")
            exit(0)

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

    def find_gamma_z_xy_util(self, x, y, z, k, Gamma, sig_samples=10000, measure="cmiknn"):
        if Gamma > 0:
            x = x[Gamma:]
            y = y[Gamma:]
            x = x.reset_index(drop=True)
            y = y.reset_index(drop=True)
            z = z.reset_index(drop=True)

            z = z[:-Gamma]
            x = x.reset_index(drop=True)
            y = y.reset_index(drop=True)
            z = z.reset_index(drop=True)
        else:
            print("Error: Gamma <= 0")
            exit(0)

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

        cmi = TestMI(p_value=False)
        _, cmi_val = cmi.fit(x, y, z)
        return cmi_val

    def find_gamma_z_xy(self, x, y, z, k, max_gamma=5, measure="cmiknn"):
        z = z.loc[y.index[0]:]

        c1 = list()

        c1.append(1)
        for G in range(1, max_gamma):
            val = self.find_gamma_z_xy_util(x, y, z, k=k, Gamma=G, measure=measure)
            c1.append(val)

        G = np.argmin(c1) + 1
        return G

    def align_pqr(self, v_p, v_q, idx_q, r, k):

        names_r = [*r.keys()]
        v_r = dict()
        nr_visted = []

        for nr in names_r:
            # idx_pq = idx_q

            v_p_new = v_p.copy()
            v_q_new = v_q.copy()
            v_p_new.index = idx_q
            v_q_new.index = idx_q
            g = self.find_gamma_z_xy(v_p_new, v_q_new, r[nr], k, max_gamma=5, measure="cmiknn")
            print("Gamma = " + str(g))
            # nr_processed = r[nr]

            # xyz_dict = {name_q: v_p_new, nr: nr_processed}
            # xyz_dict[name_q].index = idx_pq

            bool_idx = pd.DataFrame([False] * len(idx_q), columns=['bool'])
            bool_idx.index = idx_q

            v_q, r_processed, idx_q, _ = self.align_pq(v_q_new, r[nr], g)
            bool_idx.loc[idx_q] = True
            bool_idx = bool_idx['bool'].values
            v_p = v_p[bool_idx]
            # idx_p = idx_p[bool_idx]

            for nr_v in nr_visted:
                v_r[nr_v] = v_r[nr_v][bool_idx]
            v_r[nr] = r_processed
            nr_visted.append(nr)

        v_p = v_p.reset_index(drop=True)
        v_q = v_q.reset_index(drop=True)
        for nr_v in nr_visted:
            v_r[nr_v] = v_r[nr_v].reset_index(drop=True)

        return v_p, v_q, v_r

    def causation_entropy(self, p, q, r_list=[], k=10):
        gamma = self.gamma_matrix[self.names[q]].loc[self.names[p]]
        pt_1 = self.series[self.names[p]].copy()
        qt = self.series[self.names[q]].copy()
        pt_1, qt, idx_p, idx_q = self.align_pq(pt_1, qt, gamma)

        # qt = q.iloc[1:].values
        # pt_1 = p.iloc[:-1].values
        if len(r_list) > 0:
            # rt_1 = r.iloc[:-1].values
            r_1_dict = dict()
            for r_i in r_list:
                r_1_dict[self.names[r_i]] = self.series[self.names[r_i]].copy()
            pt_1, qt, r_1_dict = self.align_pqr(pt_1, qt, idx_q, r_1_dict, k)

            # Dict to df
            rt_1 = pd.DataFrame()
            for name in r_1_dict.keys():
                if isinstance(r_1_dict[name], pd.Series):
                    r_1_dict[name] = r_1_dict[name].to_frame()
                rt_1[r_1_dict[name].columns] = r_1_dict[name].reset_index(drop=True)
            rt_1 = rt_1.values
        else:
            rt_1 = None

        qt = qt.values
        pt_1 = pt_1.values

        cmi = TestMI()
        pval, val = cmi.fit(qt, pt_1, rt_1)
        return pval

    def causation_entropy_simple(self, p, q, r_list=[], k=10):
        pt_1 = self.series[self.names[p]].copy()
        qt = self.series[self.names[q]].copy()

        qt = qt.iloc[1:].values
        pt_1 = pt_1.iloc[:-1].values
        if len(r_list) > 0:
            rt_1 = self.series[self.names[r_list]].copy()
            rt_1 = rt_1.iloc[:-1].values

        else:
            rt_1 = None

        cmi = TestMI()
        pval, val = cmi.fit(qt, pt_1, rt_1)
        return pval

    def progressive_removal_of_non_causal_nodes(self):
        if self.verbose:
            print("######################################")
            print("Progressive Removal of Non-Causal Nodes")
            print("######################################")

        parents = dict()
        for q in range(self.d):
            parents[self.names[q]] = []
            for p in range(self.d):
                if p != q:
                    if self.graph.edges[p, q] == 2:
                        parents[self.names[q]].append(self.names[p])
                else:
                    parents[self.names[q]].append(self.names[q])
        print(parents)

        for q in range(self.d):
            name_q = self.series.columns[q]
            # series_q = self.series[name_q]
            parents_q = parents[name_q].copy()
            for name_p in parents_q:
                p = self.names.tolist().index(name_p)
                parents_q_without_p = list(set(parents[self.series.columns[q]]) - {name_p})
                r_list = []
                for par_name in parents_q_without_p:
                    r_list.append(self.names.tolist().index(par_name))
                # series_p = self.series[name_p]
                # series_cond = self.series[parents_q_without_p]
                print(name_p, name_q)
                pval = self.causation_entropy_simple(p, q, r_list)
                if self.verbose:
                    print('CE('+name_p+'->'+name_q+'|'+str(parents_q_without_p)+') = '+str(pval))
                if pval > self.sig_lev:
                    if self.verbose:
                        print('Remove '+name_p+' from parents of '+name_q)
                    parents[self.series.columns[q]].remove(name_p)
                    self.graph.edges[p, q] = 0
                    self.graph.edges[q, p] = 0


    def fit(self):
        """
        run KITMI
        :return: graph (CPDAG)
        """
        if self.verbose:
            now = datetime.now()
            print("#######################################")
            print("########### Starting KITMI ###########")
            print("########### " + now.strftime("%H:%M:%S" + " ###########"))
            print("#######################################")

        # Progressive Removal of Non-Causal Nodes
        self.progressive_removal_of_non_causal_nodes()

        if self.verbose:
            print("######################################")
            print("Final Results (KITMI)")
            print("######################################")
            print("Summary Graph:")
            print(self.graph.edges)
        return self.graph.edges


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
        p_list, q_list = np.where((np.triu(self.graph.edges)-np.diag(np.diag(self.graph.edges))) == 2)
        print(self.graph.edges)
        print(np.triu(self.graph.edges)-np.diag(np.diag(self.graph.edges)))
        print(p_list, q_list)
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
                (self.graph.edges[r, p] == 2) and (self.gamma_matrix[self.names[p]].loc[self.names[r]] >= 0)) or (
                (self.graph.edges[r, q] == 2) and (self.gamma_matrix[self.names[q]].loc[self.names[r]] >= 0)))]

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

    def fit2(self):
        """
        run PCTMI
        :return: graph (CPDAG)
        """
        if self.verbose:
            now = datetime.now()
            print("#######################################")
            print("########### Starting KITMI ###########")
            print("########### " + now.strftime("%H:%M:%S" + " ###########"))
            print("#######################################")

        # initialize skeleton
        self.skeleton_initialize()

        # get separation sets
        self.find_sep_set()

        if self.verbose:
            print("######################################")
            print("Final Results (KITMI)")
            print("######################################")
            print("Summary Graph:")
            print(self.graph.edges)
        return self.graph.edges


if __name__ == "__main__":
    from data.sim_data import generate_v_structure, generate_fork, diamond_generator, generate_mediator, mooij_7ts

    # data = generate_v_structure(2000)
    data = generate_fork(2000)
    # data, _, _ = diamond_generator(2000)
    # data.drop([data.columns[1]], axis=1, inplace=True)

    lag = 5
    # d = len(data.columns)

    n_iters = 1000
    hidden_size = 25  # d*(order-1)*2
    learning_rate = 0.01
    # input_size = 3

    res = tskiko_mv(data, lag, learning_rate, n_iters, noise=True, alpha=0.05)
    print(res)
    # print(res['discovery'])

    # x = mts_order(data, order=order) #[:-order]
    # print(x.loc[2:5])
    # # y = mts_order(data[order+1:], order=order)
    # names_x = x.columns[:-d]
    # names_y = x.columns[-d:]
    # y = x[names_y]
    # x = x[names_x]
    # print(x.shape)
    # print(y.shape)
