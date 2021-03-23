from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd

if __name__ == "__main__":
    type_of_causes = "other"
    for type_of_causes in ["other", "self", "ts"]:
        if type_of_causes == "other":
            measures_list = ["Precision Adjacent:", "Recall Adjacent:", "F-Score Adjacent:", "Precision Oriented:",
                             "Recall Oriented:", "F-Score Oriented:"]
            # methods_list = ['GrangerPW', 'GrangerMV', 'TCDF', 'PCMCI_CMIknn', 'PCMCI_ParCorr', 'PCTMI', 'VarLiNGAM', 'TiMINO']
            methods_list = ['GrangerPW', 'GrangerMV', 'TCDF', 'PCMCIParCorr', 'TiMINO']
            tab_index = ['P', 'R', 'F', '$\\overrightarrow{\\text{P}}$', '$\\overrightarrow{\\text{R}}$',
                                   '$\\overrightarrow{\\text{F1}}$']
        elif type_of_causes == "self":
            measures_list = ["Precision Self:", "Recall Self:", "F-Score Self:"]
            methods_list = ['GrangerPW', 'GrangerMV', 'TCDF', 'PCMCIParCorr', 'TiMINO']
            tab_index = ['SP', 'SR', 'SF']
        elif type_of_causes == "ts":
            measures_list = ["Precision Temporal:", "Recall Temporal:", "F-Score Temporal:"]
            methods_list = ['TCDF', 'PCMCIParCorr']
            tab_index = ['TP', 'TR', 'TF']

        nb_measures = len(measures_list)
        path_input = '../experiments/performance/results_'+type_of_causes
        files_input_name = [f for f in listdir(path_input) if isfile(join(path_input, f)) and not f.startswith('.')]
        print(files_input_name)
        # methods_list = []
        # for s in files_input_name:
        #     methods_list.append(s.split("_", 1)[0])
        # methods_list = np.unique(methods_list)
        nb_methods = len(methods_list)
        print(methods_list)

        datasets_list = []
        for s in files_input_name:
            datasets_list.append(s.split("_", 1)[1])
        datasets_list = np.unique(datasets_list)
        nb_datasets = len(datasets_list)
        print(datasets_list)

        results_table = pd.DataFrame(np.zeros([nb_measures, nb_methods]), columns=methods_list, index=measures_list)
        for dataset in datasets_list:
            for s in files_input_name:
                if dataset == s.split("_", 1)[1]:
                    method = s.split("_", 1)[0]
                    with open(str(path_input)+"/"+str(method)+"_"+str(dataset), "r") as f:
                        for line in f:
                            for measure in measures_list:
                                if measure in line:
                                    nextLine = next(f)
                                    nextLine = nextLine.replace("+-", "\pm")
                                    nextLine = nextLine.replace("\n", "")
                                    nextline_processed = nextLine.split(" ")
                                    nextline_processed[0] = round(float(nextline_processed[0]), 3)
                                    nextline_processed[2] = round(float(nextline_processed[2]), 3)
                                    nextline_processed = '$'+' '.join([str(elem) for elem in nextline_processed])+'$'
                                    results_table[method].loc[measure] = nextline_processed
            print(dataset)
            print(results_table)
            index_init = results_table.index
            results_table.index = tab_index
            if type_of_causes == "other":
                results_table.to_csv(r'../experiments/performance/summary/'+str(dataset)+'.txt', sep='&')
            elif type_of_causes == "self":
                results_table.to_csv(r'../experiments/performance/summary_self/' + str(dataset) + '.txt', sep='&')
            elif type_of_causes == "ts":
                results_table.to_csv(r'../experiments/performance/summary_ts/' + str(dataset) + '.txt', sep='&')
            results_table.index = index_init

