import networkx as nx
import matplotlib.pyplot as plt

# def dataframe_to_graph(g_df):
#     gtrue = nx.DiGraph()
#     gtrue.add_nodes_from(g_df.columns)
#     for name_x in g_df.columns:
#         for name_y in g_df.columns:
#             if g_df[name_y].loc[name_x] == 2:
#                 gtrue.add_edges_from([(name_x, name_y)])
#     return 1


# def _dict_to_tgraph(g_dict):
#     for name_y in g_dict.keys():
#         for name_x, t_xy in g_dict[name_y]:
#             if (name_x, name_y) in g.edges:
#                 tg.edges[name_x, name_y]['time'].append(-t_xy)
#             else:
#                 tg.add_edges_from([(name_x, name_y)])
#                 tg.edges[name_x, name_y]['time'] = [-t_xy]


def tgraph_to_graph(tg):
    g = nx.DiGraph()
    og = nx.DiGraph()
    sg = nx.DiGraph()
    g.add_nodes_from(tg.nodes)
    og.add_nodes_from(tg.nodes)
    sg.add_nodes_from(tg.nodes)
    for cause, effects in tg.adj.items():
        for effect, _ in effects.items():
            if cause != effect:
                og.add_edges_from([(cause, effect)])
                g.add_edges_from([(cause, effect)])
            else:
                sg.add_edges_from([(cause, effect)])
                g.add_edges_from([(cause, effect)])
    return g, og, sg


def tgraph_to_list(tg):
    list_tg = []
    for cause, effects in tg.adj.items():
        for effect, eattr in effects.items():
            t_list = eattr['time']
            for t in t_list:
                list_tg.append((cause, effect, t))
    return list_tg


def print_graph(g):
    for cause, effects in g.adj.items():
        for effect, eattr in effects.items():
            print('(%s -> %s)' % (cause, effect))


def print_temporal_graph(tg):
    list_tg_hat = tgraph_to_list(tg)
    for edge in list_tg_hat:
        print('(%s -> %s with t= %d)' % (edge[0], edge[1], edge[2]))


def draw_graph(g, node_size=300):
    pos = nx.spring_layout(g, k=0.25, iterations=20)
    nx.draw(g, pos, with_labels=True, font_weight='bold', node_size=node_size)
    # nx.draw_shell(g, nlist=[range(4)], with_labels=True, font_weight='bold')
    plt.show()


def draw_temporal_graph(tg, node_size=300):
    pos = nx.spring_layout(tg, k=0.25, iterations=20)
    nx.draw(tg, pos, with_labels=True, font_weight='bold', node_size=node_size)
    edge_labels = nx.get_edge_attributes(tg, 'time')
    nx.draw_networkx_edge_labels(tg, pos, labels=edge_labels)
    plt.show()


def string_nodes(nodes):
    new_nodes = []
    for col in nodes:
        try:
            int(col)
            new_nodes.append("V" + str(col))
        except ValueError:
            new_nodes.append(col)
    return new_nodes


if __name__ == "__main__":
    import networkx as nx

    g = nx.read_gpickle("../experiments/graphs/fork/GrangerPW_1")
    print_graph(g)
