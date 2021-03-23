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


def run_timino_pw_R(arg_list):
    script = "timino"
    # Remove all arguments from directory
    dir_path = os.path.dirname(os.path.realpath(__file__))
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
        # g_df.index = g_df.columns
        print(g_df)
        return g_df.values
    else:
        print('R Error:\n {0}'.format(error))
        exit(0)


if __name__ == "__main__":
    structure = "mediator"
    print(os.getcwd())
    path = "./data/simulated_ts_data/unscaled/"+str(structure)+"/data_"+str(0)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)

    df = run_timino_pw_R([[data, "data"], [0.05, "alpha"], [5, "nlags"]])
    print(df)
