import pandas as pd

from causalnex.structure.dynotears import from_pandas_dynamic


def dynotears(data, tau_max=5, alpha=0.0):
    graph_dict = dict()
    for name in data.columns:
        graph_dict[name] = []

    sm = from_pandas_dynamic(data, p=tau_max, w_threshold=0.01, lambda_w=0.05, lambda_a=0.05)

    # print(sm.edges)
    # print(sm.pred)

    tname_to_name_dict = dict()
    count_lag = 0
    idx_name = 0
    for tname in sm.nodes:
        tname_to_name_dict[tname] = data.columns[idx_name]
        if count_lag == tau_max:
            idx_name = idx_name +1
            count_lag = -1
        count_lag = count_lag +1

    for ce in sm.edges:
        c = ce[0]
        e = ce[1]
        tc = int(c.partition("lag")[2])
        te = int(e.partition("lag")[2])
        t = tc - te
        if (tname_to_name_dict[c], -t) not in graph_dict[tname_to_name_dict[e]]:
            graph_dict[tname_to_name_dict[e]].append((tname_to_name_dict[c], -t))

    # g = sm.to_directed()
    return graph_dict

if __name__ == "__main__":
    structure = "diamond"
    path = "../../data/simulated_ts_data/"+str(structure)+"/data_"+str(0)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)
    data = data.loc[:1000]
    g = dynotears(data, 5)

    print(g)

