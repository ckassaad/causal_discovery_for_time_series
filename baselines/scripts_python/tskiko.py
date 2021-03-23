import pandas as pd
from baselines.scripts_python.python_packages.TsKIKO.tskiko_mv import tskiko as tsk


def tskiko(data, sig_level=0.05, tau_max=5, verbose=True):
    g_df = tsk(data, max_lag=tau_max, learning_rate=0.01, training_epoch=5000, noise=True, alpha=sig_level,
        cond_ind_test="ParCorr", verbose=verbose)

    return g_df


if __name__ == "__main__":
    import os
    structure = "diamond"
    print(os.getcwd())
    path = "../../data/simulated_ts_data/"+str(structure)+"/data_"+str(0)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)

    df = tskiko(data, sig_level=0.05, maxlag=5, verbose=True)
    print(df)

