# THIS SOURCE CODE IS SUPPLIED AS IS WITHOUT WAR
# RANTY OF ANY KIND	 AND ITS AUTHOR AND THE JOURNAL OF
# ARTIFICIAL INTELLIGENCE RESEARCH JAIR AND JAIRS PUB
#LISHERS AND DISTRIBUTORS	 DISCLAIM ANY AND ALL WARRANTIES	
# INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE	 
# AND ANY WARRANTIES OR NON INFRINGEMENT THE USER
# ASSUMES ALL LIABILITY AND RESPONSIBILITY FOR USE OF THIS
# SOURCE CODE	 AND NEITHER THE AUTHOR NOR JAIR	 NOR JAIRS
# PUBLISHERS AND DISTRIBUTORS	 WILL BE LIABLE FOR DAM
# AGES OF ANY KIND RESULTING FROM ITS USE Without limiting
# the generality of the foregoing	 neither the author	 nor JAIR	 nor JAIR's 
# publishers and distributors	 warrant that the Source Code will be errorfree	
# will operate without interruption	 or will meet the needs of the user

import time
from os import listdir
from os.path import isfile, join

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import causal_discovery_class as cd
import networkx as nx

from tools.graph_functions import print_graph, print_temporal_graph, tgraph_to_graph, string_nodes


def two_col_format_to_graphs(nodes, two_col_format):
    tgtrue = nx.DiGraph()
    tgtrue.add_nodes_from(nodes)
    for i in range(two_col_format.shape[0]):
        c = nodes[int(two_col_format[i, 0])]
        e = nodes[int(two_col_format[i, 1])]
        tgtrue.add_edges_from([(c, e)])
    gtrue, ogtrue, sgtrue = tgraph_to_graph(tgtrue)
    return gtrue, ogtrue, sgtrue, tgtrue


def run_on_data(method, verbose):
    save_model = True

    path_input = './data/dairy_markets/'

    cheese = pd.read_csv('./data/dairy_markets/Butter.csv', delimiter=',', index_col=0, header=0)
    butter = pd.read_csv('./data/dairy_markets/Cheese.csv', delimiter=',', index_col=0, header=0)
    milk = pd.read_csv('./data/dairy_markets/Milk.csv', delimiter=',', index_col=0, header=0)

    date_start = '2008-09-30'
    milk = milk.loc[date_start:]
    butter = butter.loc[date_start:]
    cheese = cheese.loc[date_start:]
    print(cheese)

    data = pd.concat([milk, butter, cheese], axis=1, sort=False)
    data.columns = ["milk", "butter", "cheese"]
    data = data.reset_index(drop=True)

    nodes = string_nodes(data.columns)


    if verbose:
        print("d = " + str(data.shape[1]))
        print("T = " + str(data.shape[0]))

    three_col_format_ground_truth = np.loadtxt('./data/dairy_markets/ground_truth/gt_dairy_markets.csv', delimiter=',')
    gtrue, ogtrue, sgtrue, tgtrue = two_col_format_to_graphs(nodes, three_col_format_ground_truth)

    start = time.time()

    if method == "GrangerPW":
        model = cd.GrangerPW(nodes, sig_level=0.05, nlags=5)
        model.infer_from_data(data)
    elif method == "GrangerMV":
        model = cd.GrangerMV(nodes, sig_level=0.05, nlags=5)
        model.infer_from_data(data)
    elif method == "TCDF":
        model = cd.TCDF(nodes, epochs=1000,  kernel_size=4, dilation_coefficient=4, hidden_layers=1, learning_rate=0.01,
                    sig_level=0.05)
        model.infer_from_data(data)
    elif method == "PCMCICMIknn":
        model = cd.PCMCI(nodes, sig_level=0.05, nlags=5, cond_ind_test="CMIknn")
        model.infer_from_data(data)
    elif method == "PCMCIParCorr":
        model = cd.PCMCI(nodes, sig_level=0.05, nlags=5, cond_ind_test="ParCorr")
        model.infer_from_data(data)
    elif method == "oCSE":
        model = cd.OCSE(nodes, sig_level=0.05)
        model.infer_from_data(data)
    elif method == "PCTMI":
        model = cd.PCTMI(nodes, sig_level=0.05, nlags=5)
        model.infer_from_data(data)
    elif method == "NBCB_pw":
        model = cd.NBCB(nodes, sig_level=0.05, nlags=5)
        model.infer_from_data(data)
    elif method == "NBCB":
        model = cd.NBCB(nodes, sig_level=0.05, nlags=5, pairwise=False)
        model.infer_from_data(data)
    elif method == "PWNBCBk":
        model = cd.PWNBCBk(nodes, sig_level=0.05, nlags=5, pairwise=False)
        model.infer_from_data(data)
    elif method == "tsFCI":
        model = cd.TsFCI(nodes, sig_level=0.05, nlags=5)
        model.infer_from_data(data)
    elif method == "FCITMI":
        model = cd.FCITMI(nodes, sig_level=0.05, nlags=5)
        model.infer_from_data(data)
    elif method == "VarLiNGAM":
        model = cd.VarLiNGAM(nodes, sig_level=0.05, nlags=5)
        model.infer_from_data(data)
    elif method == "TiMINO":
        model = cd.TiMINO(nodes, sig_level=0.05, nlags=5)
        model.infer_from_data(data)
    elif method == "tsKIKO":
        model = cd.TsKIKO(nodes, sig_level=0.05, nlags=5)
        model.infer_from_data(data)
    elif method == "Dynotears":
        model = cd.Dynotears(nodes, sig_level=0.05, nlags=5)
        model.infer_from_data(data)
    else:
        model = None
        print("Error: method not found")
        exit(0)

    end = time.time()

    if save_model:
        nx.write_gpickle(model.ghat, "./experiments/graphs/summary_other_and_self_graphs/dairy/"+method)
        nx.write_gpickle(model.oghat, "./experiments/graphs/summary_other_graphs/dairy/"+method)
        nx.write_gpickle(model.sghat, "./experiments/graphs/summary_self_graphs/dairy/"+method)

    # evaluation self and other
    pres_a = model.evaluation(gtrue, evaluation_measure="precision_adjacent")
    rec_a = model.evaluation(gtrue, evaluation_measure="recall_adjacent")
    fscore_a = model.evaluation(gtrue, evaluation_measure="f1_adjacent")
    pres_o = model.evaluation(gtrue, evaluation_measure="precision_oriented")
    rec_o = model.evaluation(gtrue, evaluation_measure="recall_oriented")
    fscore_o = model.evaluation(gtrue, evaluation_measure="f1_oriented")
    # evaluation other
    o_pres_a = model.evaluation(ogtrue, evaluation_measure="other_precision_adjacent")
    o_rec_a = model.evaluation(ogtrue, evaluation_measure="other_recall_adjacent")
    o_fscore_a = model.evaluation(ogtrue, evaluation_measure="other_f1_adjacent")
    o_pres_o = model.evaluation(ogtrue, evaluation_measure="other_precision_oriented")
    o_rec_o = model.evaluation(ogtrue, evaluation_measure="other_recall_oriented")
    o_fscore_o = model.evaluation(ogtrue, evaluation_measure="other_f1_oriented")
    # evaluation self
    s_pres = model.evaluation(sgtrue, evaluation_measure="self_precision")
    s_rec = model.evaluation(sgtrue, evaluation_measure="self_recall")
    s_fscore = model.evaluation(sgtrue, evaluation_measure="self_f1")

    if verbose:
        # print(nodes)
        print('True Graph')
        print_graph(gtrue)
        # print('True Graph Self')
        # print_graph(sgtrue)
        # print('True Graph Other')
        # print_graph(ogtrue)
        print('Inferred Graph')
        model.print_graph()

        print("precision adjacent: " + str(pres_a))
        print("recall adjacent: " + str(rec_a))
        print("f-score adjacent: " + str(fscore_a))
        print("precision oriented: " + str(pres_o))
        print("recall oriented: " + str(rec_o))
        print("f-score oriented: " + str(fscore_o))

        print("other precision adjacent: " + str(o_pres_a))
        print("other recall adjacent: " + str(o_rec_a))
        print("other f-score adjacent: " + str(o_fscore_a))
        print("other precision oriented: " + str(o_pres_o))
        print("other recall oriented: " + str(o_rec_o))
        print("other f-score oriented: " + str(o_fscore_o))

        print("self precision: " + str(s_pres))
        print("self recall self: " + str(s_rec))
        print("self f-score self: " + str(s_fscore))

    print("Computation time: " + str(end - start))
    return pres_a, rec_a, fscore_a, pres_o, rec_o, fscore_o, o_pres_a, o_rec_a, o_fscore_a, o_pres_o, o_rec_o, \
           o_fscore_o, s_pres, s_rec, s_fscore, (end - start)



if __name__ == "__main__":
    import sys

    if len(sys.argv) > 2:
        print(len(sys.argv))
        method = sys.argv[1]  # GrangerPW, GrangerMV, TCDF, PCMCICMIknn, PCMCIParCorr, PCTMI, tsFCI, FCITMI VarLiNGAM, TiMINO
        num_processor = int(sys.argv[2])  # -1 for all
        verbose = bool(int(sys.argv[3]))
        print('Argument List:', str(sys.argv))
    else:
        print('Missing arguments so will take default arguments')
        method = "Dynotears"  # GrangerPW, GrangerMV, TCDF, PCMCICMIknn, PCMCIParCorr, oCSE, PCTMI, tsFCI, FCITMI VarLiNGAM, TiMINO
        num_processor = 1
        verbose = True
        print('Default Argument List:', str(method), num_processor)

    results = run_on_data(method, verbose)

    results = np.array(results).reshape(1, -1)
    pres_a_list = results[:, 0]
    rec_a_list = results[:, 1]
    fscore_a_list = results[:, 2]
    pres_o_list = results[:, 3]
    rec_o_list = results[:, 4]
    fscore_o_list = results[:, 5]
    o_pres_a_list = results[:, 6]
    o_rec_a_list = results[:, 7]
    o_fscore_a_list = results[:, 8]
    o_pres_o_list = results[:, 9]
    o_rec_o_list = results[:, 10]
    o_fscore_o_list = results[:, 11]
    s_pres_list = results[:, 12]
    s_rec_list = results[:, 13]
    s_fscore_list = results[:, 14]
    comput_time_list = results[:, 15]

    # method = method+"window=1"
    method = method+"window=auto"
    with open("./experiments/performance_average/summary_other_and_self_performance_average/" + str(method) + "_dairy", "w+") as file:
        file.write("Precision Adjacent: \n" + str(np.mean(pres_a_list)) + " +- " +
                   str(np.std(pres_a_list)))
        file.write("\n")
        file.write("Recall Adjacent: \n" + str(np.mean(rec_a_list)) + " +- " + str(np.std(rec_a_list)))
        file.write("\n")
        file.write("F-Score Adjacent: \n" + str(np.mean(fscore_a_list)) + " +- " + str(np.std(fscore_a_list)))
        file.write("\n")
        file.write("Precision Oriented: \n" + str(np.mean(pres_o_list)) + " +- " + str(np.std(pres_o_list)))
        file.write("\n")
        file.write("Recall Oriented: \n" + str(np.mean(rec_o_list)) + " +- " + str(np.std(rec_o_list)))
        file.write("\n")
        file.write("F-Score Oriented: \n" + str(np.mean(fscore_o_list)) + " +- " + str(np.std(fscore_o_list)))
        file.write("\n")

        file.write("\n\nComputational Time: " + str(np.mean(comput_time_list)) + " +- " + str(np.std(comput_time_list)))

    with open(
            "./experiments/performance_average/summary_other_performance_average/" + str(method) + "_dairy",
            "w+") as file:
        file.write("Other Precision Adjacent: \n" + str(np.mean(o_pres_a_list)) + " +- " +
                   str(np.std(o_pres_a_list)))
        file.write("\n")
        file.write("Other Recall Adjacent: \n" + str(np.mean(o_rec_a_list)) + " +- " + str(np.std(o_rec_a_list)))
        file.write("\n")
        file.write("Other F-Score Adjacent: \n" + str(np.mean(o_fscore_a_list)) + " +- " + str(np.std(o_fscore_a_list)))
        file.write("\n")
        file.write("Other Precision Oriented: \n" + str(np.mean(o_pres_o_list)) + " +- " + str(np.std(o_pres_o_list)))
        file.write("\n")
        file.write("Other Recall Oriented: \n" + str(np.mean(o_rec_o_list)) + " +- " + str(np.std(o_rec_o_list)))
        file.write("\n")
        file.write("Other F-Score Oriented: \n" + str(np.mean(o_fscore_o_list)) + " +- " + str(np.std(o_fscore_o_list)))
        file.write("\n")

        file.write("\n\nComputational Time: " + str(np.mean(comput_time_list)) + " +- " + str(np.std(comput_time_list)))

    with open(
            "./experiments/performance_average/summary_self_performance_average/" + str(method) + "_dairy",
            "w+") as file:
        file.write("Self Precision: \n" + str(np.mean(s_pres_list)) + " +- " + str(np.std(s_pres_list)))
        file.write("\n")
        file.write("Self Recall: \n" + str(np.mean(s_rec_list)) + " +- " + str(np.std(s_rec_list)))
        file.write("\n")
        file.write("Self F-Score: \n" + str(np.mean(s_fscore_list)) + " +- " + str(np.std(s_fscore_list)))
        file.write("\n")

        file.write("\n\nComputational Time: " + str(np.mean(comput_time_list)) + " +- " + str(np.std(comput_time_list)))

    if verbose:
        print("####################### Final Result #######################")
        print("Precision Adjacent: " + str(np.mean(pres_a_list)) + " +- " + str(np.std(pres_a_list)))
        print("Recall Adjacent: " + str(np.mean(rec_a_list)) + " +- " + str(np.std(rec_a_list)))
        print("F-Score Adjacent: " + str(np.mean(fscore_a_list)) + " +- " + str(np.std(fscore_a_list)))
        print("Precision Oriented: " + str(np.mean(pres_o_list)) + " +- " + str(np.std(pres_o_list)))
        print("Recall Oriented: " + str(np.mean(rec_o_list)) + " +- " + str(np.std(rec_o_list)))
        print("F-Score Oriented: " + str(np.mean(fscore_o_list)) + " +- " + str(np.std(fscore_o_list)))
        print("Other Precision Adjacent: " + str(np.mean(o_pres_a_list)) + " +- " + str(np.std(o_pres_a_list)))
        print("Other Recall Adjacent: " + str(np.mean(o_rec_a_list)) + " +- " + str(np.std(o_rec_a_list)))
        print("Other F-Score Adjacent: " + str(np.mean(o_fscore_a_list)) + " +- " + str(np.std(o_fscore_a_list)))
        print("Other Precision Oriented: " + str(np.mean(o_pres_o_list)) + " +- " + str(np.std(o_pres_o_list)))
        print("Other Recall Oriented: " + str(np.mean(o_rec_o_list)) + " +- " + str(np.std(o_rec_o_list)))
        print("Other F-Score Oriented: " + str(np.mean(o_fscore_o_list)) + " +- " + str(np.std(o_fscore_o_list)))
        print("Self Precision: " + str(np.mean(s_pres_list)) + " +- " + str(np.std(s_pres_list)))
        print("Self Recall: " + str(np.mean(s_rec_list)) + " +- " + str(np.std(s_rec_list)))
        print("Self F-Score: " + str(np.mean(s_fscore_list)) + " +- " + str(np.std(s_fscore_list)))
        print("Computational Time: " + str(np.mean(comput_time_list)) + " +- " + str(np.std(comput_time_list)))

    np.savetxt("./experiments/performance_detail/" + str(method) + "_dairy.csv", results,
               delimiter=';', header="pres_a, rec_a, fscore_a, pres_o, rec_o, fscore_o, o_pres_a, o_rec_a, o_fscore_a, o_pres_o, o_rec_o, \
               o_fscore_o, s_pres, s_rec, s_fscore, computational_time")

