import numpy as np
import pandas as pd
from baselines.scripts_python.python_packages.TsKIKO.tools.tigramite.tigramite.independence_tests import CMIknn


class TestMI:
    def __init__(self):
        self.cd = CMIknn(mask_type=None, significance='shuffle_test', fixed_thres=None, sig_samples=10000,
                sig_blocklength=3, knn=10, confidence='bootstrap', conf_lev=0.9, conf_samples=10000,
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
        pvalue = self.cd.get_shuffle_significance(X.T, xyz, value)
        return pvalue, value


def causation_entropy(p, q, r=None):
    qt = q.iloc[1:].values
    pt_1 = p.iloc[:-1].values
    if r is not None:
        rt_1 = r.iloc[:-1].values
    else:
        rt_1 = None

    cmi = TestMI()
    pval, val = cmi.fit(qt, pt_1, rt_1)
    return pval


def ocse(data, sig_level=0.05, verbose=True):
    parents = dict()
    for q in range(data.shape[1]):
        name_q = data.columns[q]
        if verbose:
            print("Select " + name_q)
            print("Aggregative Discovery of Causal Nodes")
        pval = 0
        # parents[name_q] = [name_q]
        parents[name_q] = []
        series_q = data[name_q]
        # not_parents_q = list(set(data.columns) - {name_q})
        not_parents_q = list(data.columns)
        while (pval <= sig_level) and (len(not_parents_q) > 0):
            pval_list = []
            for name_p in not_parents_q:
                print(not_parents_q)
                print(name_p)
                print(parents[name_q])
                series_p = data[name_p]
                series_cond = data[parents[name_q]]
                pval = causation_entropy(series_p, series_q, series_cond)
                pval_list.append(pval)
            print(pval_list)
            pval = np.min(pval_list)
            # p = np.argmin(pval_list)
            p_list = list(np.argwhere(pval_list == pval)[:, 0])
            names_p = []
            for p in p_list:
                names_p.append(not_parents_q[p])
            if verbose:
                print('test indeps :' + str(pval_list))
                print('CE('+str(names_p)+'->'+name_q+'|'+str(parents[name_q])+') = '+str(pval))
            if pval <= sig_level:
                if verbose:
                    print(str(names_p)+'->'+name_q)
                for name_p in names_p:
                    parents[data.columns[q]].append(name_p)
                    not_parents_q.remove(name_p)

        if verbose:
            print("Progressive Removal of Non-Causal Nodes")
        # parents_q = list(set(parents[name_q]) - {name_q})
        parents_q = parents[name_q].copy()
        for name_p in parents_q:
            parents_q_without_p = list(set(parents[data.columns[q]]) - {name_p})
            series_p = data[name_p]
            series_cond = data[parents_q_without_p]
            pval = causation_entropy(series_p, series_q, series_cond)
            if verbose:
                print('CE('+name_p+'->'+name_q+'|'+str(parents_q_without_p)+') = '+str(pval))
            if pval > sig_level:
                if verbose:
                    print('Remove '+name_p+' from parents of '+name_q)
                parents[data.columns[q]].remove(name_p)

    parents_df = pd.DataFrame(np.zeros([data.shape[1], data.shape[1]], dtype=np.int8), columns=data.columns,
                              index=data.columns)
    for name_q in parents.keys():
        parents_df[name_q].loc[name_q] = 1
        for name_p in parents[name_q]:
            if name_q == name_p:
                parents_df[name_q].loc[name_p] = 1
            else:
                parents_df[name_q].loc[name_p] = 2
                if parents_df[name_p].loc[name_q] == 0:
                    parents_df[name_p].loc[name_q] = 1
    print(parents)
    return parents_df


if __name__ == "__main__":
    import os
    structure = "diamond"
    print(os.getcwd())
    path = "../../data/simulated_ts_data/"+str(structure)+"/data_"+str(0)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)
    data = data.loc[:250]

    df = ocse(data, sig_level=0.05, verbose=True)
    print(df)
