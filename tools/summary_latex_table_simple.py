from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd

if __name__ == "__main__":
    measures_list = ["Precision Adjacent:", "Recall Adjacent:",
                     "F-Score Adjacent:", "Precision Oriented:", "Recall Oriented:", "F-Score Oriented:"]
    nb_measures = len(measures_list)
    # path_input = './results'
    path_input = '../experiments/performance_average/summary_other_performance_average/'
    files_input_name = [f for f in listdir(path_input) if isfile(join(path_input, f)) and not f.startswith('.')]
    print(files_input_name)
    # methods_list = []
    # for s in files_input_name:
    #     methods_list.append(s.split("_", 1)[0])
    # methods_list = np.unique(methods_list)
    # nb_methods = len(methods_list)
    # print(methods_list)
    methods_list = ['PCTMIwindow=auto', 'PCMCICMIknn',  'PCMCIParCorr', 'TiMINO', 'TCDF', 'GrangerMV']
    nb_methods = len(methods_list)
    print(methods_list)

    # datasets_list = []
    # for s in files_input_name:
    #     datasets_list.append(s.split("_", 1)[1])
    # datasets_list = np.unique(datasets_list)
    # nb_datasets = len(datasets_list)
    datasets_list = ["v_structure", "fork", "mediator", "diamond"]
    nb_datasets = len(datasets_list)
    print(datasets_list)
    # datasets_list = ["Finance"]

    results_table = pd.DataFrame(np.zeros([nb_measures, nb_methods]), columns=methods_list, index=measures_list)
    for dataset in datasets_list:
        for s in files_input_name:
            # if dataset == s.split("_", 1)[1]:
            method = s.split("_", 1)[0]
            if method in methods_list:
                with open(str(path_input)+"/"+str(method)+"_"+str(dataset)+"_1000", "r") as f:
                    for line in f:
                        for measure in measures_list:
                            if measure in line:
                                nextLine = next(f)
                                nextLine = nextLine.replace("+-", "\pm")
                                nextLine = nextLine.replace("\n", "")
                                nextline_processed = nextLine.split(" ")
                                nextline_processed[0] = round(float(nextline_processed[0]), 2)
                                nextline_processed[2] = round(float(nextline_processed[2]), 2)
                                nextline_processed = '$'+' '.join([str(elem) for elem in nextline_processed])+'$'
                                results_table[method].loc[measure] = nextline_processed
        print(dataset)
        print(results_table)
        index_init = results_table.index
        results_table.index = ['P', 'R', 'F', '$\\overrightarrow{\\text{P}}$', '$\\overrightarrow{\\text{R}}$', '$\\overrightarrow{\\text{F1}}$']
        results_table.to_csv(r'./sssssssssssssssss'+dataset+'.txt', sep='&')
        results_table.index = index_init

