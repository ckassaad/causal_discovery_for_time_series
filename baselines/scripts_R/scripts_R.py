from subprocess import Popen, PIPE
import os
import glob
import pandas as pd


def clear_args(dir_path):
    files = glob.glob(dir_path + '/args/*')
    for f in files:
        os.remove(f)


def clear_results(dir_path):
    files = glob.glob(dir_path + '/results/*')
    for f in files:
        os.remove(f)


def run_R(script, arg_list):
    if script not in ["timino", "varlingam", "tsfci"]:
        print("Error: method not available")
        exit(0)
    # Remove all arguments from directory
    dir_path = os.path.dirname(os.path.realpath(__file__))
    script_name = script
    script = dir_path + "/" + script + ".R"
    clear_args(dir_path)
    clear_results(dir_path)
    r_arg_list = []
    # COMMAND WITH ARGUMENTS
    for a in arg_list:
        if isinstance(a[0], pd.DataFrame):
            a[0].to_csv(dir_path + "/args/"+a[1]+".csv", index=False)
            r_arg_list.append(dir_path + "/args/" + a[1] + ".csv")
        if isinstance(a[0], int):
            f = open(dir_path + "/args/"+a[1]+".txt", "w")
            f.write(str(a[0]))
            f.close()
            r_arg_list.append(dir_path + "/args/" + a[1] + ".txt")
        if isinstance(a[0], float):
            f = open(dir_path + "/args/"+a[1]+".txt", "w")
            f.write(str(a[0]))
            f.close()
            r_arg_list.append(dir_path + "/args/" + a[1] + ".txt")

    r_arg_list.append(dir_path)
    cmd = ["Rscript", script] + r_arg_list

    p = Popen(cmd, cwd="./", stdin=PIPE, stdout=PIPE, stderr=PIPE)
    # Return R output or error
    output, error = p.communicate()
    print(output)
    if p.returncode == 0:
        print('R Done')
        g_df = pd.read_csv(dir_path + "/results/result.csv", header=0, index_col=0)
        print(g_df)
        if script_name == "tsfci":
            g_dict = ts_fci_dataframe_to_dict(g_df, arg_list[0][0].columns.tolist(), arg_list[2][0])
            return g_dict, g_df
        else:
            return g_df, g_df
    else:
        print('R Error:\n {0}'.format(error))
        exit(0)


def ts_fci_dataframe_to_dict(df, names, nlags):
    # todo: check if its correct
    for i in range(df.shape[1]):
        for j in range(i+1, df.shape[1]):
            if df[df.columns[i]].loc[df.columns[j]] == 2:
                if df[df.columns[j]].loc[df.columns[i]] == 2:
                    print(df.columns[i] + " <-> " + df.columns[j])

    g_dict = dict()
    for name_y in names:
        g_dict[name_y] = []
    for ty in range(nlags):
        for name_y in names:
            t_name_y = df.columns[ty*len(names)+names.index(name_y)]
            for tx in range(nlags):
                for name_x in names:
                    t_name_x = df.columns[tx * len(names) + names.index(name_x)]
                    if df[t_name_y].loc[t_name_x] == 2:
                        if (name_x, tx-ty) not in g_dict[name_y]:
                            g_dict[name_y].append((name_x, tx - ty))
                    # if (name_x, ty - tx, ">") not in g_dict[name_y]:
                    #     g_dict[name_y].append((name_x, ty-tx, ">"))
                    # elif df[t_name_y].loc[t_name_x] == 3:
                    #     if (name_x, ty - tx, "o") not in g_dict[name_y]:
                    #         g_dict[name_y].append((name_x, ty - tx, "o"))
    print(g_dict)
    return g_dict


if __name__ == "__main__":
    structure = "mediator"
    print(os.getcwd())
    path = "../../data/simulated_ts_data/"+str(structure)+"/data_"+str(0)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)
    print(data)

    # graph = run_R("timino.R", [[data, "data"], [0.05, "alpha"], [5, "nlags"]])
    # run_R("varlingam.R", [[data, "data"], [5, "nlags"]])
    df = run_R("varlingam", [[data, "data"], [0.05, "alpha"], [5, "nlags"]])
    # print(df)
