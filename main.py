import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from sklearn.metrics import auc

from file_methods import *
from network_properties import plot_total_prc, plot_node_auprc, plot_rtf_found
from utils import reader
from pyrwr.ppr import PPR
import networkx as nx

G = nx.Graph()
forward_G = nx.Graph()

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

datas = ["Alpha6Beta4Integrin", "AndrogenReceptor", "BCR", "BDNF", "CRH", "EGFR1", "FSH", "Hedgehog", "IL1",
         "IL2", "IL3", "IL4", "IL5", "IL6", "IL9", "IL-7", "KitReceptor", "Leptin", "Notch", "Oncostatin_M",
         "Prolactin", "RANKL", "TCR", "TGF_beta_Receptor", "TNFalpha", "TSH", "TSLP", "TWEAK", "Wnt"]

input_graph = "data/interactome-weights.txt"
graph_type = "directed"

''' 
These are some boolean variables that define what do we expect from the code 
'''
HAS_CLEANED_PATHWAY = True
RUN_ALGORITHMS = False
# output
PLOT_EDGES_PRC = False
PLOT_NODES_PRC = False
COMPUTE_RTF = True
# writing methods
WRITE_EDGES = False
WRITE_PRC = False
WRITE_NODES_TO_ID_MAP = False

overall_recalls_ours = []
overall_precisions_ours = []
overall_recalls_pl = []
overall_precisions_pl = []
overall_recalls_rwr = []
overall_precisions_rwr = []

overall_node_recalls_ours = []
overall_node_precisions_ours = []
overall_node_recalls_pl = []
overall_node_precisions_pl = []
overall_node_recalls_rwr = []
overall_node_precisions_rwr = []

overall_r_ours = []
overall_tf_ours = []
overall_r_pl = []
overall_tf_pl = []
overall_r_rwr = []
overall_tf_rwr = []


# This method computes the edge weights based on the given method
def compute_new_graph(method, alpha=0.0):
    if method == "edgeLinker":
        return np.array(
            [[int(edge[0]), int(edge[1]), math.pow(2, -alpha * float(edge[2] + length[edge[1]]))] for edge in
             graph])
    if method == "m2":
        return np.array(
            [[int(edge[0]), int(edge[1]), math.pow(2, -alpha * float(length[edge[1]]))] for edge in graph])
    # This is our main weighting method
    if method == "ours":
        return np.array(
            [[int(edge[0]), int(edge[1]), math.pow(2, -alpha * float(edge[2] + length[edge[1]]))] for edge in graph])
    if method == "rwr":
        return np.array([[int(edge[0]), int(edge[1]), math.pow(2, -float(edge[2]))] for edge in graph])


def compute_recall_precision(sorted_edges, known_pathway, k, direction, recall_bound=1.0):
    recalls = []
    precisions = []
    tps = []
    fps = []
    tp = 0
    fp = 0
    for i in range(min(k, len(sorted_edges))):
        edge = sorted_edges[i]
        if [edge[0][0], edge[0][1]] in known_pathway or (not direction and [edge[0][1], edge[0][0]] in known_pathway):
            tp += 1
            tps.append(1)
            fps.append(0)
        else:
            fp += 1
            tps.append(0)
            fps.append(1)
        recalls.append(tp / len(known_pathway))
        precisions.append(tp / (tp + fp))
        if recalls[-1] > recall_bound:
            break
    return recalls, precisions, tps, fps


def compute_edge_probs(adj_list, r):
    edge_probs = []
    for node in nodes.values():
        total_weight = sum([edge[1] for edge in adj_list[node]])
        for edge in adj_list[node]:
            edge_probs.append([[int(node), int(edge[0])], r[int(node)] * edge[1] / total_weight])
    return sorted(edge_probs, key=lambda x: -x[1])


def run_random_walk(new_graph, c=0.15):
    ppr = PPR()
    ppr.read_graph(new_graph, True)
    return ppr.compute(seeds, c=c, max_iters=1000)


def run_algorithm(dataset, method, color, alpha=0.0, c=0.15, k=1000000, recall_bound=1.0, direction=True):
    if RUN_ALGORITHMS:
        new_graph = compute_new_graph(method, alpha)
        # creating the graph using the new weights
        adj_list = {}
        for node in nodes.values():
            adj_list[int(node)] = []
        for i in range(len(new_graph)):
            adj_list[new_graph[i][0]].append([new_graph[i][1], new_graph[i][2]])

        # running random walk to obtain node probabilities
        print("running random walk for " + method)
        r = run_random_walk(new_graph, c)

        # computing the edge probabilities
        print("computing the edge probabilities for " + method)
        sorted_edges = compute_edge_probs(adj_list, r)
        result = [[sorted_edges[i][0][0], sorted_edges[i][0][1]] for i in range(min(k, len(sorted_edges)))]
        if WRITE_EDGES:
            write_edges(result, "results/" + dataset + "edges-" + method + ".txt")
    else:
        result = read_edges("results/" + dataset + "edges-" + method + ".txt")
        sorted_edges = [[[r[0], r[1]], 0] for r in result]

    recalls = []
    precisions = []

    if PLOT_EDGES_PRC:
        # computing the precision and recall
        print("computing recall-precision curve for " + method)
        recalls, precisions, tps, fps = compute_recall_precision(sorted_edges, subpathway, k, direction, recall_bound)
        if WRITE_PRC:
            write_precision_recall(precisions, recalls, "results/" + dataset + "PR-" + method + ".txt")

        name = method
        plt.plot(recalls, precisions, color=color, label=name + " " + str(round(auc(recalls, precisions), 4)))

        print("AUPRC of " + method + ": " + str(auc(recalls, precisions)))

    return result, recalls, precisions


def read_pathway(path):
    paths = []
    with open(path, 'r') as f:
        for line in f:
            splitted = line.split()
            if splitted[0] != "#tail" and splitted[1] in nodes and splitted[0] in nodes:
                paths.append([nodes[splitted[1]], nodes[splitted[0]]])
    return paths


def read_pathlinker_output(path):
    edges = []
    with open(path, 'r') as f:
        for line in f:
            sp = line.split()
            if sp[0] != "#tail":
                edge = [nodes[sp[1]], nodes[sp[0]]]
                if edge not in edges:
                    edges.append(edge)
    return [[edge] for edge in edges]


def add_pathlinker(path, color, k=1000000, direction=True):
    edges = read_pathlinker_output(path)
    recalls, precisions, tps, fps = compute_recall_precision(edges, subpathway, k, direction)

    name = "PathLinker"
    plt.plot(recalls, precisions, color=color, label=name + " " + str(round(auc(recalls, precisions), 4)))

    print("AUPRC of pathlinker: " + str(auc(recalls, precisions)))

    # computing the highest ranked edges
    return [[edges[i][0][0], edges[i][0][1]] for i in range(min(k, len(edges)))], len(edges), recalls, precisions


def clean_pathway(path, dataset):
    subpath = []
    for p in path:
        for edge in graph:
            if edge[0] == p[0] and edge[1] == p[1]:
                subpath.append(p)
    with open("data/cleaned_pathways/" + dataset + "-edges.txt", 'w') as f:
        for p in subpath:
            f.write(id_to_node[p[0]] + " " + id_to_node[p[1]] + "\n")
    return subpath


def read_cleaned_pathway(dataset):
    paths = []
    with open("data/cleaned_pathways/" + dataset + "-edges.txt", 'r') as f:
        for line in f:
            splitted = line.split()
            if splitted[0] != "#tail" and splitted[1] in nodes and splitted[0] in nodes:
                paths.append([nodes[splitted[0]], nodes[splitted[1]]])
    return paths


# reading the graph
_, nodes, graph, id_to_node = reader.read_graph(input_graph, graph_type)
if WRITE_NODES_TO_ID_MAP:
    write_nodes("data/nodes_map.txt", nodes)

# creating the networkx graph
if RUN_ALGORITHMS:
    G.add_nodes_from(nodes)
    for e in graph:
        G.add_edge(int(e[1]), int(e[0]), weight=e[2])

    for e in graph:
        forward_G.add_edge(int(e[0]), int(e[1]), weight=e[2])

total_pathway_lengths = 0
total_pathway_node_lengths = 0
for data in datas:
    rtf_path = "data/NetPath/" + data + "-nodes.txt"
    pathway_path = "data/NetPath/" + data + "-edges.txt"
    pathlinker = "data/PathLinker_output/" + data + "k-2000-ranked-edges.txt"

    print("Pathway: " + data)

    # reading seeds and targets
    seeds, targets = read_source_and_destinations(rtf_path, nodes)

    # reading pathway
    if HAS_CLEANED_PATHWAY:
        subpathway = read_cleaned_pathway(data)
    else:
        pathway = read_pathway(pathway_path)
        subpathway = clean_pathway(pathway, data)

    total_pathway_lengths += len(subpathway)

    # finding the distances of each node from the targets
    if RUN_ALGORITHMS:
        length = nx.multi_source_dijkstra_path_length(G, targets)
        forward_length = nx.multi_source_dijkstra_path_length(forward_G, seeds)

    # running the algorithms and get the pathways, true positives, and false positives
    pathway_pathlinker, pl_edge_len, pl_recalls, pl_precisions = add_pathlinker(pathlinker, color=colors["black"],
                                                                                direction=False)
    pathway_ours, our_recalls, our_precisions = run_algorithm(dataset=data, method="ours", color=colors["deepskyblue"],
                                                              alpha=5, c=0.25, k=pl_edge_len, direction=False)
    pathway_rwr, rwr_recalls, rwr_precisions = run_algorithm(dataset=data, method="rwr", color=colors["silver"],
                                                             k=pl_edge_len, direction=False)
    if PLOT_EDGES_PRC:
        overall_recalls_ours.append(our_recalls)
        overall_precisions_ours.append(our_precisions)
        overall_recalls_pl.append(pl_recalls)
        overall_precisions_pl.append(pl_precisions)
        overall_recalls_rwr.append(rwr_recalls)
        overall_precisions_rwr.append(rwr_precisions)

        plt.legend()
        plt.title("recall-precision for " + data)
        plt.savefig("output/edge-PRC/" + data + ".png")
        plt.close()

    if COMPUTE_RTF:
        our_seeds, our_targets, rwr_seeds, rwr_targets, pl_seeds, pl_targets = plot_rtf_found(seeds, targets,
                                                                                              pathway_ours, pathway_rwr,
                                                                                              pathway_pathlinker, data)
        overall_r_ours.append(our_seeds)
        overall_tf_ours.append(our_targets)
        overall_r_pl.append(pl_seeds)
        overall_tf_pl.append(pl_targets)
        overall_r_rwr.append(rwr_seeds)
        overall_tf_rwr.append(rwr_targets)

    # Node AUPRC
    if PLOT_NODES_PRC:
        our_recalls, our_precisions, rwr_recalls, rwr_precisions, pl_recalls, pl_precisions, node_len = plot_node_auprc(
            subpathway, pathway_ours, pathway_rwr, pathway_pathlinker, data)
        overall_node_recalls_ours.append(our_recalls)
        overall_node_precisions_ours.append(our_precisions)
        overall_node_recalls_pl.append(pl_recalls)
        overall_node_precisions_pl.append(pl_precisions)
        overall_node_recalls_rwr.append(rwr_recalls)
        overall_node_precisions_rwr.append(rwr_precisions)
        total_pathway_node_lengths += node_len

if PLOT_EDGES_PRC:
    plot_total_prc(overall_recalls_ours, overall_precisions_ours, overall_recalls_rwr, overall_precisions_rwr,
                   overall_recalls_pl, overall_precisions_pl, "overall-edge-PRC")

if PLOT_NODES_PRC:
    plot_total_prc(overall_node_recalls_ours, overall_node_precisions_ours, overall_node_recalls_rwr,
                   overall_node_precisions_rwr, overall_node_recalls_pl, overall_node_precisions_pl, "overall-node-PRC")

if COMPUTE_RTF:
    plot_total_prc([[(i + 1) for i in range(len(r_ours))] for r_ours in overall_r_ours], overall_r_ours,
                   [[(i + 1) for i in range(len(r_rwr))] for r_rwr in overall_r_rwr], overall_r_rwr,
                   [[(i + 1) for i in range(len(r_pl))] for r_pl in overall_r_pl], overall_r_pl, "overall-receptors")

    plot_total_prc([[(i + 1) for i in range(len(r_ours))] for r_ours in overall_tf_ours], overall_tf_ours,
                   [[(i + 1) for i in range(len(r_rwr))] for r_rwr in overall_tf_rwr], overall_tf_rwr,
                   [[(i + 1) for i in range(len(r_pl))] for r_pl in overall_tf_pl], overall_tf_pl, "overall-tfs")
