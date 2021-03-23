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


def run_matlab(script, arg_list):
    if script not in ["granger_mv", "globalmit", "fblg"]:
        print("Error: method not available")
        exit(0)
    # Remove all arguments from directory
    dir_path = os.path.dirname(os.path.realpath(__file__))
    script_name = script
    script = dir_path + "/" + script #+ ".m"
    clear_args(dir_path)
    clear_results(dir_path)
    r_arg_list = []
    # COMMAND WITH ARGUMENTS
    for a in arg_list:
        if isinstance(a[0], pd.DataFrame):
            nodes = a[0].columns
            a[0].to_csv(dir_path + "/args/"+a[1]+".csv", index=False)
            r_arg_list.append(dir_path + "/args/" + a[1] + ".csv")
            # r_arg_list.append(a[1])
        if isinstance(a[0], int):
            f = open(dir_path + "/args/"+a[1]+".txt", "w")
            f.write(str(a[0]))
            f.close()
            r_arg_list.append(dir_path + "/args/" + a[1] + ".txt")
            # r_arg_list.append(a[1])
        if isinstance(a[0], float):
            f = open(dir_path + "/args/"+a[1]+".txt", "w")
            f.write(str(a[0]))
            f.close()
            r_arg_list.append(dir_path + "/args/" + a[1] + ".txt")
            # r_arg_list.append(a[1])

    r_arg_list.append(dir_path)
    # matlab - nodesktop - nosplash - r 'try granger_mv("data", "sig_level", "nlags") ; catch; end; quit'
    cmd = ['matlab', '-nodesktop', '-nosplash', '-r',  'try '+script_name+'("'+str(r_arg_list[0])+'","'+str(r_arg_list[1])+'","' +str(r_arg_list[2])+'","' +str(r_arg_list[3])+'"); catch; end; quit']
    print(cmd)

    working_path = os.getcwd()
    os.chdir(dir_path)
    p = Popen(cmd, cwd="./", stdin=PIPE, stdout=PIPE, stderr=PIPE)
    # Return matlab output or error
    output, error = p.communicate()
    os.chdir(working_path)
    print(output)
    if p.returncode == 0:
        print('matlab Done')
        g_df = pd.read_csv(dir_path + "/results/result.txt")
        g_df.columns = nodes
        g_df.index = nodes
        return g_df
    else:
        print('matlab Error:\n {0}'.format(error))
        exit(0)


if __name__ == "__main__":
    structure = "diamond"
    print(os.getcwd())
    path = "../../data/simulated_ts_data/"+str(structure)+"/data_"+str(0)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)

    df = run_matlab("granger_mv", [[data, "data"], [0.05, "sig_level"], [5, "nlags"]])
    print(df)
