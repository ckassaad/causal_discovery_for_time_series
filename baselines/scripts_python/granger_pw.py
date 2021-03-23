import numpy as np
import pandas as pd
from statsmodels.tsa.stattools import grangercausalitytests


def granger_pw(X_train, sig_level=0.05, maxlag=5, verbose=False):
    test = 'ssr_ftest'
    names = X_train.columns
    dataset = pd.DataFrame(np.zeros((len(names), len(names)), dtype=int), columns=names, index=names)
    for c in dataset.columns:
        for r in dataset.index:
            test_result = grangercausalitytests(X_train[[r,c]], maxlag=maxlag, verbose=verbose)
            p_values = [round(test_result[i+1][0][test][1], 4) for i in range(maxlag)]
            min_p_value = np.min(p_values)
            # dataset.loc[c, r] = min_p_value
            if min_p_value < sig_level:
                dataset.loc[c, r] = 2

    for c in dataset.columns:
        for r in dataset.index:
            if dataset.loc[c, r] == 2:
                if dataset.loc[r, c] == 0:
                    dataset.loc[r, c] = 1
            if r == c:
                dataset.loc[r, c] = 1
    # dataset.columns = [var + '_y' for var in variables]
    # dataset.index = [var + '_x' for var in variables]
    return dataset


if __name__ == "__main__":
    import os
    structure = "diamond"
    print(os.getcwd())
    path = "../../data/simulated_ts_data/"+str(structure)+"/data_"+str(0)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)

    df = granger_pw(data, sig_level=0.05, maxlag=5, verbose=False)
    print(df)

