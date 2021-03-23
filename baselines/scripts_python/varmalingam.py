from baselines.scripts_python.python_packages.lingam_master.lingam.varma_lingam import VARMALiNGAM

import numpy as np

def varmalingam(data, tau_max=5, alpha=0.05):
    min_causal_effect = alpha
    split_by_causal_effect_sign = True

    model = VARMALiNGAM(order=(tau_max, tau_max), criterion='bic', prune=True)
    model.fit(data)

    print(model.causal_order_)

    m = model._adjacency_matrices
    am = np.concatenate([*m], axis=1)

    dag = np.abs(am) > min_causal_effect

    if split_by_causal_effect_sign:
        direction = np.array(np.where(dag))
        signs = np.zeros_like(dag).astype('int64')
        for i, j in direction.T:
            signs[i][j] = np.sign(am[i][j]).astype('int64')
        dag = signs

    dag = np.abs(dag)
    names = data.columns
    print(dag)
    res_dict = dict()
    for e in range(dag.shape[0]):
        res_dict[names[e]] = []
    for c in range(dag.shape[0]):
        for te in range(dag.shape[1]):
            if dag[c][te] == 1:
                e = te%dag.shape[0]
                t = te//dag.shape[0]
                res_dict[names[e]].append((names[c], -t))
    return res_dict


if __name__ == "__main__":
    import pandas as pd
    structure = "fork"
    path = "../../data/simulated_ts_data/"+str(structure)+"/data_"+str(0)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)
    res = varmalingam(data.loc[:100], tau_max=1, alpha=0.00)

    # print(res)
