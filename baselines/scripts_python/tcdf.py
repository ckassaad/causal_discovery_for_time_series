from subprocess import Popen, PIPE
import os
import glob
import pandas as pd
import json


def clear_args(dir_path):
    files = glob.glob(dir_path + '/args/*')
    for f in files:
        os.remove(f)


def clear_results(dir_path):
    files = glob.glob(dir_path + '/results/*')
    for f in files:
        os.remove(f)


def tcdf(arg_list):
    # Remove all arguments from directory
    dir_path = os.path.dirname(os.path.realpath(__file__))
    script = dir_path + "/python_packages/TCDF-master/runTCDF" + ".py"
    clear_args(dir_path)
    clear_results(dir_path)
    r_arg_list = []
    # COMMAND WITH ARGUMENTS
    for a in arg_list:
        if isinstance(a[0], pd.DataFrame):
            a[0].to_csv(dir_path + "/args/" + a[1] + ".csv", index=False)
            r_arg_list.append("--"+a[1])
            r_arg_list.append(dir_path + "/args/" + a[1] + ".csv")
        if isinstance(a[0], int):
            # f = open(dir_path + "/args/" + a[1] + ".txt", "w")
            # f.write(str(a[0]))
            # f.close()
            r_arg_list.append("--"+a[1])
            r_arg_list.append(str(a[0]))
        if isinstance(a[0], float):
            # f = open(dir_path + "/args/" + a[1] + ".txt", "w")
            # f.write(str(a[0]))
            # f.close()
            r_arg_list.append("--"+a[1])
            r_arg_list.append(str(a[0]))

    r_arg_list.append("--path")
    r_arg_list.append(dir_path)
    cmd = ["python", script] + r_arg_list

    p = Popen(cmd, cwd="./", stdin=PIPE, stdout=PIPE, stderr=PIPE)
    # Return R output or error
    output, error = p.communicate()
    print(output)
    if p.returncode == 0:
        print('Python Done')
        g_dict = json.load(open(dir_path + "/results/tcdf_result.txt"))
        for key in g_dict.keys():
            key_list = []
            for elem in g_dict[key]:
                key_list.append(tuple(elem))
            g_dict[key] = key_list
        return g_dict
    else:
        print('Python Error:\n {0}'.format(error))
        exit(0)


if __name__ == "__main__":
    structure = "fork"
    print(os.getcwd())
    path = "../../data/simulated_ts_data/"+str(structure)+"/data_"+str(3)+".csv"
    data = pd.read_csv(path, delimiter=',', index_col=0)

    model = tcdf([[data, "data"], [1000, "epochs"], [4, "kernel_size"], [1, "hidden_layers"], [0.01, "learning_rate"],
                  [4, "kernel_size"], [4, "dilation_coefficient"], [0.05, "significance"]])
    print(model)
