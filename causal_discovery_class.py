import networkx as nx
from baselines.scripts_python.pcmci import pcmci
from baselines.scripts_python.varlingam import varlingam
from baselines.scripts_python.tcdf import tcdf
from baselines.scripts_python.granger_pw import granger_pw
from baselines.scripts_python.pctmi import pctmi, tpctmi
from baselines.scripts_python.fcitmi import fcitmi
from baselines.scripts_python.tskiko import tskiko
from baselines.scripts_python.ocse import ocse

from baselines.scripts_python.nbcb import nbcb

try:
    from baselines.scripts_R.scripts_R import run_R
except:
    print("Could not import R packages for TiMINO and tsFCI")
try:
    from baselines.scripts_matlab.scripts_matlab import run_matlab
except:
    print("Could not import Matlab packages for MVGranger")

from tools.graph_functions import tgraph_to_graph, tgraph_to_list, print_graph, print_temporal_graph, draw_graph, \
    draw_temporal_graph, string_nodes


class GraphicalModel:
    def __init__(self, nodes):
        super(GraphicalModel, self).__init__()
        nodes = string_nodes(nodes)
        self.ghat = nx.DiGraph()
        self.ghat.add_nodes_from(nodes)
        self.oghat = nx.DiGraph()
        self.oghat.add_nodes_from(nodes)
        self.sghat = nx.DiGraph()
        self.sghat.add_nodes_from(nodes)

    def infer_from_data(self, data):
        raise NotImplementedError

    def _dataframe_to_graph(self, df):
        for name_x in df.columns:
            if df[name_x].loc[name_x] > 0:
                self.sghat.add_edges_from([(name_x, name_x)])
                self.ghat.add_edges_from([(name_x, name_x)])
            for name_y in df.columns:
                if name_x != name_y:
                    if df[name_y].loc[name_x] == 2:
                        self.oghat.add_edges_from([(name_x, name_y)])
                        self.ghat.add_edges_from([(name_x, name_y)])

    def draw(self, node_size=300):
        draw_graph(self.ghat, node_size=node_size)

    def print_graph(self):
        print_graph(self.ghat)

    def print_other_graph(self):
        print_graph(self.oghat)

    def print_self_graph(self):
        print_graph(self.sghat)

    def evaluation(self, gtrue, evaluation_measure="f1_oriented"):
        if evaluation_measure == "precision_adjacent":
            res = self._precision(gtrue, method="all_adjacent")
        elif evaluation_measure == "recall_adjacent":
            res = self._recall(gtrue, method="all_adjacent")
        elif evaluation_measure == "f1_adjacent":
            res = self._f1(gtrue, method="all_adjacent")
        elif evaluation_measure == "precision_oriented":
            res = self._precision(gtrue, method="all_oriented")
        elif evaluation_measure == "recall_oriented":
            res = self._recall(gtrue, method="all_oriented")
        elif evaluation_measure == "f1_oriented":
            res = self._f1(gtrue, method="all_oriented")
        elif evaluation_measure == "other_precision_adjacent":
            res = self._precision(gtrue, method="other_adjacent")
        elif evaluation_measure == "other_recall_adjacent":
            res = self._recall(gtrue, method="other_adjacent")
        elif evaluation_measure == "other_f1_adjacent":
            res = self._f1(gtrue, method="other_adjacent")
        elif evaluation_measure == "other_precision_oriented":
            res = self._precision(gtrue, method="other_oriented")
        elif evaluation_measure == "other_recall_oriented":
            res = self._recall(gtrue, method="other_oriented")
        elif evaluation_measure == "other_f1_oriented":
            res = self._f1(gtrue, method="other_oriented")
        elif evaluation_measure == "self_precision":
            res = self._precision(gtrue, method="self")
        elif evaluation_measure == "self_recall":
            res = self._recall(gtrue, method="self")
        elif evaluation_measure == "self_f1":
            res = self._f1(gtrue, method="self")
        else:
            raise AttributeError(evaluation_measure)
        return res

    def _hamming_distance(self, gtrue):
        # todo: check if it's correct (maybe it's not truely hamming distance)
        res = nx.graph_edit_distance(self.ghat, gtrue)
        return 1 - res/max(self.ghat.number_of_edges(), gtrue.number_of_edges())

    def _tp(self, gtrue, method="all_oriented"):  # oriented or adjacent
        if method == "all_oriented":
            tp = nx.intersection(gtrue, self.ghat)
        elif method == "all_adjacent":
            undirected_true = gtrue.to_undirected()
            undirected_hat = self.ghat.to_undirected()
            tp = nx.intersection(undirected_true, undirected_hat)
        elif method == "other_oriented":
            tp = nx.intersection(gtrue, self.oghat)
        elif method == "other_adjacent":
            undirected_true = gtrue.to_undirected()
            undirected_hat = self.oghat.to_undirected()
            tp = nx.intersection(undirected_true, undirected_hat)
        elif method == "self":
            tp = nx.intersection(gtrue, self.sghat)
        else:
            raise AttributeError(method)
        return len(tp.edges)

    def _fp(self, gtrue, method="all_oriented"):  # oriented or adjacent
        if method == "all_oriented":
            fp = nx.difference(self.ghat, gtrue)
        elif method == "all_adjacent":
            undirected_true = gtrue.to_undirected()
            undirected_hat = self.ghat.to_undirected()
            fp = nx.difference(undirected_hat, undirected_true)
        elif method == "other_oriented":
            fp = nx.difference(self.oghat, gtrue)
        elif method == "other_adjacent":
            undirected_true = gtrue.to_undirected()
            undirected_hat = self.oghat.to_undirected()
            fp = nx.difference(undirected_hat, undirected_true)
        elif method == "self":
            fp = nx.difference(self.sghat, gtrue)
        else:
            raise AttributeError(method)
        return len(fp.edges)

    def _fn(self, gtrue, method="all_oriented"):  # oriented or adjacent
        if method == "all_oriented":
            fn = nx.difference(gtrue, self.ghat)
        elif method == "all_adjacent":
            undirected_true = gtrue.to_undirected()
            undirected_hat = self.ghat.to_undirected()
            fn = nx.difference(undirected_true, undirected_hat)
        elif method == "other_oriented":
            fn = nx.difference(gtrue, self.oghat)
        elif method == "other_adjacent":
            undirected_true = gtrue.to_undirected()
            undirected_hat = self.oghat.to_undirected()
            fn = nx.difference(undirected_true, undirected_hat)
        elif method == "self":
            fn = nx.difference(gtrue, self.sghat)
        else:
            raise AttributeError(method)
        return len(fn.edges)

    def _topology(self, gtrue, method="all_oriented"):
        correct = self._tp(gtrue, method)
        added = self._fp(gtrue, method)
        missing = self._fn(gtrue, method)
        return correct/(correct + missing + added)

    def _false_positive_rate(self, gtrue, method="all_oriented"):
        true_pos = self._tp(gtrue, method)
        false_pos = self._fp(gtrue, method)
        if false_pos == 0:
            return 0
        else:
            return false_pos / (true_pos + false_pos)

    def _precision(self, gtrue, method="all_oriented"):
        true_pos = self._tp(gtrue, method)
        false_pos = self._fp(gtrue, method)
        if true_pos == 0:
            return 0
        else:
            return true_pos / (true_pos + false_pos)

    def _recall(self, gtrue, method="all_oriented"):
        true_pos = self._tp(gtrue, method)
        false_neg = self._fn(gtrue, method)
        if true_pos == 0:
            return 0
        else:
            return true_pos / (true_pos + false_neg)

    def _f1(self, gtrue, method="all_oriented"):
        p = self._precision(gtrue, method)
        r = self._recall(gtrue, method)
        if (p == 0) and (r == 0):
            return 0
        else:
            return 2 * p * r / (p + r)


class TemporalGraphicalModel(GraphicalModel):
    def __init__(self, nodes):
        GraphicalModel.__init__(self, nodes)
        nodes = string_nodes(nodes)
        self.tghat = nx.DiGraph()
        self.tghat.add_nodes_from(nodes)

    def infer_from_data(self, data):
        raise NotImplementedError

    def _dict_to_tgraph(self, temporal_dict):
        for name_y in temporal_dict.keys():
            for name_x, t_xy in temporal_dict[name_y]:
                if (name_x, name_y) in self.tghat.edges:
                    self.tghat.edges[name_x, name_y]['time'].append(-t_xy)
                else:
                    self.tghat.add_edges_from([(name_x, name_y)])
                    self.tghat.edges[name_x, name_y]['time'] = [-t_xy]
                # self.TGhat.add_weighted_edges_from([(name_x, name_y, t_xy)])

    def _tgraph_to_graph(self):
        self.ghat, self.oghat, self.sghat = tgraph_to_graph(self.tghat)

    def draw_temporal_graph(self, node_size=300):
        draw_temporal_graph(self.tghat, node_size=node_size)

    def print_temporal_graph(self):
        print_temporal_graph(self.tghat)

    def temporal_evaluation(self, tgtrue, evaluation_measure="f1"):
        if evaluation_measure == "precision":
            res = self._temporal_precision(tgtrue)
        elif evaluation_measure == "recall":
            res = self._temporal_recall(tgtrue)
        elif evaluation_measure == "f1":
            res = self._temporal_f1(tgtrue)
        else:
            raise AttributeError(evaluation_measure)
        return res

    def _temporal_tp(self, tgtrue):
        list_tg_true = tgraph_to_list(tgtrue)
        list_tg_hat = tgraph_to_list(self.tghat)
        tp = set(list_tg_true).intersection(list_tg_hat)
        return len(tp)

    def _temporal_fp(self, tgtrue):
        list_tg_true = tgraph_to_list(tgtrue)
        list_tg_hat = tgraph_to_list(self.tghat)
        fp = set(list_tg_hat).difference(list_tg_true)
        return len(fp)

    def _temporal_fn(self, tgtrue):
        list_tg_true = tgraph_to_list(tgtrue)
        list_tg_hat = tgraph_to_list(self.tghat)
        fn = set(list_tg_true).difference(list_tg_hat)
        return len(fn)

    def _temporal_false_positive_rate(self, tgtrue):
        true_pos = self._temporal_tp(tgtrue)
        false_pos = self._temporal_fp(tgtrue)
        return false_pos / (true_pos + false_pos)

    def _temporal_precision(self, tgtrue):
        true_pos = self._temporal_tp(tgtrue)
        false_pos = self._temporal_fp(tgtrue)
        if true_pos == 0:
            return 0
        else:
            return true_pos / (true_pos + false_pos)

    def _temporal_recall(self, tgtrue):
        true_pos = self._temporal_tp(tgtrue)
        false_neg = self._temporal_fn(tgtrue)
        if true_pos == 0:
            return 0
        else:
            return true_pos / (true_pos + false_neg)

    def _temporal_f1(self, tgtrue):
        p = self._temporal_precision(tgtrue)
        r = self._temporal_recall(tgtrue)
        if (p == 0) and (r == 0):
            return 0
        else:
            return 2 * p * r / (p + r)


class GrangerPW(GraphicalModel):
    def __init__(self, nodes, sig_level=0.05, nlags=5):
        GraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level
        self.nlags = nlags

    def infer_from_data(self, data):
        data.columns = list(self.ghat.nodes)
        g_df = granger_pw(data, sig_level=self.sig_level, maxlag=self.nlags, verbose=False)
        self._dataframe_to_graph(g_df)


class GrangerMV(GraphicalModel):
    def __init__(self, nodes, sig_level=0.05, nlags=5):
        GraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level
        self.nlags = nlags

    def infer_from_data(self, data):
        data.columns = list(self.ghat.nodes)
        g_df = run_matlab("granger_mv", [[data, "data"], [self.sig_level, "sig_level"], [self.nlags, "nlags"]])
        self._dataframe_to_graph(g_df)


class TCDF(TemporalGraphicalModel):
    def __init__(self, nodes, epochs=5000,  kernel_size=4, dilation_coefficient=4, hidden_layers=1, learning_rate=0.01,
                 sig_level=0.05):
        """
        TCDF algorithm
        :param nodes:
        :param epochs:
        :param kernel_size:
        :param dilation_coefficient:
        :param hidden_layers:
        :param learning_rate:
        :param sig_level:
        nlags = 1+\sum_{l=0}^{hidden_layers}(kernel_size - 1).dilation_coefficient^l
        """
        TemporalGraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level
        self.epochs = epochs
        self.kernel_size = kernel_size
        self.dilation_coefficient = dilation_coefficient  # recommended to be equal to kernel size
        self.hidden_layers = hidden_layers
        self.learning_rate = learning_rate

    def infer_from_data(self, data):
        data.columns = list(self.tghat.nodes)
        g_dict = tcdf([[data, "data"], [self.epochs, "epochs"], [self.kernel_size, "kernel_size"],
                       [self.dilation_coefficient, "dilation_coefficient"], [self.hidden_layers, "hidden_layers"],
                       [self.learning_rate, "learning_rate"], [self.sig_level, "significance"]])
        self._dict_to_tgraph(g_dict)
        self._tgraph_to_graph()


class PCMCI(TemporalGraphicalModel):
    def __init__(self, nodes, sig_level=0.05, nlags=5, cond_ind_test="CMIknn"):
        TemporalGraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level
        self.nlags = nlags
        self.cond_ind_test = cond_ind_test

    def infer_from_data(self, data):
        data.columns = list(self.tghat.nodes)
        tg_dict = pcmci(data, tau_max=self.nlags, cond_ind_test=self.cond_ind_test, alpha=self.sig_level)
        self._dict_to_tgraph(tg_dict)
        self._tgraph_to_graph()


class OCSE(GraphicalModel):
    def __init__(self, nodes, sig_level=0.05):
        GraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level

    def infer_from_data(self, data):
        data.columns = list(self.ghat.nodes)
        g_df = ocse(data, sig_level=self.sig_level)
        self._dataframe_to_graph(g_df)


class TsFCI(TemporalGraphicalModel):
    def __init__(self, nodes, sig_level=0.05, nlags=5):
        TemporalGraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level
        self.nlags = nlags

    def infer_from_data(self, data):
        data.columns = list(self.tghat.nodes)
        g_dict, init_obj = run_R("tsfci", [[data, "data"], [self.sig_level, "sig_level"], [self.nlags, "nlags"]])
        self._dict_to_tgraph(g_dict)
        self._tgraph_to_graph()


class VarFCI(TemporalGraphicalModel):
    def __init__(self, nodes):
        TemporalGraphicalModel.__init__(self, nodes)

    def infer_from_data(self, data):
        1


class PCTMI(GraphicalModel):
    def __init__(self, nodes, sig_level=0.05, nlags=5):
        GraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level
        self.nlags = nlags

    def infer_from_data(self, data):
        data.columns = list(self.ghat.nodes)
        g_df = pctmi(data, sig_level=self.sig_level, nlags=self.nlags, verbose=True)
        self._dataframe_to_graph(g_df)

class TPCTMI(TemporalGraphicalModel):
    def __init__(self, nodes, sig_level=0.05, nlags=5):
        TemporalGraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level
        self.nlags = nlags

    def infer_from_data(self, data):
        data.columns = list(self.ghat.nodes)
        tg_dict = tpctmi(data, sig_level=self.sig_level, nlags=self.nlags, verbose=True)
        self._dict_to_tgraph(tg_dict)
        self._tgraph_to_graph()

class NBCB(GraphicalModel):
    def __init__(self, nodes, sig_level=0.05, nlags=5):
        GraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level
        self.nlags = nlags

    def infer_from_data(self, data):
        data.columns = list(self.ghat.nodes)
        g_df = nbcb(data, sig_level=self.sig_level, nlags=self.nlags, verbose=True)
        self._dataframe_to_graph(g_df)


class FCITMI(GraphicalModel):
    def __init__(self, nodes, sig_level=0.05, nlags=5):
        GraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level
        self.nlags = nlags

    def infer_from_data(self, data):
        data.columns = list(self.ghat.nodes)
        g_df = fcitmi(data, sig_level=self.sig_level, nlags=self.nlags, verbose=True)
        self._dataframe_to_graph(g_df)


class VarLiNGAM(TemporalGraphicalModel):
    def __init__(self, nodes, sig_level=0.05, nlags=5):
        TemporalGraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level
        self.nlags = nlags

    def infer_from_data(self, data):
        data.columns = list(self.ghat.nodes)
        # g_df, _ = run_R("varlingam", [[data, "data"], [self.sig_level, "sig_level"], [self.nlags, "nlags"]])
        # self._dataframe_to_graph(g_df)
        # alpha is not a test, it's a causal stength
        tg_dict = varlingam(data, tau_max=self.nlags, alpha=self.sig_level)
        self._dict_to_tgraph(tg_dict)
        self._tgraph_to_graph()



class TiMINO(GraphicalModel):
    def __init__(self, nodes, sig_level=0.05, nlags=5):
        GraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level
        self.nlags = nlags

    def infer_from_data(self, data):
        data.columns = list(self.ghat.nodes)
        g_df, _ = run_R("timino", [[data, "data"], [self.sig_level, "sig_level"], [self.nlags, "nlags"]])
        self._dataframe_to_graph(g_df)


class TsKIKO(GraphicalModel):
    def __init__(self, nodes, sig_level=0.05, nlags=5):
        GraphicalModel.__init__(self, nodes)
        self.sig_level = sig_level
        self.nlags = nlags

    def infer_from_data(self, data):
        data.columns = list(self.ghat.nodes)
        g_df = tskiko(data, tau_max=self.nlags, sig_level=self.sig_level)
        self._dataframe_to_graph(g_df)

# class GIMME(GraphicalModel):
#     def __init__(self):
#         1
#
#
# class FASK(GraphicalModel):
#     def __init__(self):
#         1
#
#
# class TwoStep(GraphicalModel):
#     def __init__(self):
#         1
#
#
# class AITIA(TemporalGraphicalModel):
#     def __init__(self):
#         1


if __name__ == "__main__":
    import pandas as pd
    structure = "fork"
    path = "./data/simulated_ts_data/"+str(structure)+"/data_"+str(0)+".csv"
    dataset = pd.read_csv(path, delimiter=',', index_col=0)
    # model = TiMINO(data.columns)
    # model = TsFCI(data.columns)
    # model = VarLiNGAM(data.columns)
    model = PCMCI(dataset.columns, cond_ind_test="ParCorr")
    # model = TCDF(dataset.columns)
    # model = GrangerMV(dataset.columns)
    # model = GrangerPW(dataset.columns)
    model.infer_from_data(dataset)
    g = model.ghat

    # model.draw(node_size=1000)
    # model.print_graph()

    # Gtrue = nx.DiGraph()
    # Gtrue.add_nodes_from(dataset.columns)
    # Gtrue.add_edges_from([("V1", "V2"), ("V1", "V3")])
    # print(model._f1(Gtrue, method="oriented"))
    #
    # model.temporal_draw(node_size=1000)

    model.tghat.edges["V1", "V3"]['time'] = [1]
    model.print_temporal_graph()
    TGtrue = nx.DiGraph()
    TGtrue.add_nodes_from(dataset.columns)
    TGtrue.add_edges_from([("V1", "V1"), ("V1", "V2"), ("V1", "V3"), ("V2", "V2"), ("V3", "V3")])
    TGtrue.edges["V1", "V1"]['time'] = [1]
    TGtrue.edges["V2", "V2"]['time'] = [1]
    TGtrue.edges["V3", "V3"]['time'] = [1]
    TGtrue.edges["V1", "V2"]['time'] = [1]
    TGtrue.edges["V1", "V3"]['time'] = [2]

    print(model._temporal_f1(TGtrue))


