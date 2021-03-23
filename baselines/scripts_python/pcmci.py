from baselines.scripts_python.python_packages.tigramite.tigramite.pcmci import PCMCI
from baselines.scripts_python.python_packages.tigramite.tigramite.independence_tests import ParCorr, CMIknn
from baselines.scripts_python.python_packages.tigramite.tigramite import data_processing as pp
import numpy as np
import pandas as pd


def pcmci(data, tau_max=5, cond_ind_test="CMIknn", alpha=0.05):
    if cond_ind_test == "CMIknn":
        cond_ind_test = CMIknn()
    elif cond_ind_test == "ParCorr":
        cond_ind_test = ParCorr()

    data_tigramite = pp.DataFrame(data.values, var_names=data.columns)

    pcmci = PCMCI(
        dataframe=data_tigramite,
        cond_ind_test=cond_ind_test,
        verbosity=1)
    pcmci.run_pcmci(tau_min=0, tau_max=tau_max, pc_alpha=alpha)

    res_dict = dict()
    for effect in pcmci.all_parents.keys():
        res_dict[pcmci.var_names[effect]] = []
        for cause, t in pcmci.all_parents[effect]:
            res_dict[pcmci.var_names[effect]].append((pcmci.var_names[cause], t))

    # res_summary_array = np.zeros([data.shape[1], data.shape[1]])
    #
    # for k in pcmci.all_parents.keys():
    #     temp = pcmci.all_parents[k]
    #     temp = np.unique([x[0] for x in temp])
    #     for i in temp:
    #         if k == i:
    #             res_summary_array[k, i] = 1
    #         else:
    #             if res_summary_array[k, i] == 0:
    #                 res_summary_array[k, i] = 1
    #             res_summary_array[i, k] = 2
    # res_summary_array = pd.DataFrame(res_summary_array, columns=data.columns, index=data.columns)
    return res_dict


if __name__ == "__main__":
    import pandas as pd
    structure = "fork"
    path = "../data/simulated_ts_data/"+str(structure)+"/data_"+str(0)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)
    model = pcmci(data)
