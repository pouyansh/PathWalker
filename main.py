import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from sklearn.metrics import auc

from file_methods import *
from network_properties import plot_total_prc, plot_node_auprc, plot_rtf_found
from pathway_files import read_cleaned_pathway, clean_pathway, read_pathway, read_source_and_destinations, \
    read_pathway_names, clean_receptors_and_tfs
from utils import reader
from pyrwr.ppr import PPR
import networkx as nx

G = nx.Graph()
forward_G = nx.Graph()

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

pallet = [colors['black'], colors['orangered'], colors['grey'], colors['lightgray']]
# pallet = [colors['navy'], colors['orangered'], colors['dodgerblue'], colors['lightsteelblue']]

input_graph = "data/interactome.txt"
graph_type = "directed"

''' 
These are some boolean variables that define what do we expect from the code 
'''
DATABASE = "KEGG"
# methods to run
RUN_ALGORITHMS = False
RUN_EDGE_LINKER = False
JUST_CLEAN = True
# output
PLOT_EDGES_PRC = False
PLOT_NODES_PRC = False
COMPUTE_RTF = False
# writing methods
WRITE_EDGES = False
WRITE_EDGES_EDGE_LINKER = False
WRITE_PRC = False
WRITE_NODES_TO_ID_MAP = False
# which methods to include
INCLUDE_PATHLINKER = False
INCLUDE_RWR = False
INCLUDE_EDGE_LINKER = False
# pathway
HAS_CLEANED_PATHWAY = False

overall_recalls_ours = []
overall_precisions_ours = []
overall_recalls_pl = []
overall_precisions_pl = []
overall_recalls_rwr = []
overall_precisions_rwr = []
overall_recalls_el = []
overall_precisions_el = []

overall_node_recalls_ours = []
overall_node_precisions_ours = []
overall_node_recalls_pl = []
overall_node_precisions_pl = []
overall_node_recalls_rwr = []
overall_node_precisions_rwr = []
overall_node_recalls_el = []
overall_node_precisions_el = []

overall_r_ours = []
overall_tf_ours = []
overall_r_pl = []
overall_tf_pl = []
overall_r_rwr = []
overall_tf_rwr = []
overall_r_el = []
overall_tf_el = []


# This method computes the edge weights based on the given method
def compute_new_graph(method, alpha=0.0):
    if method == "edge_linker":
        return np.array(
            [[int(edge[0]), int(edge[1]), math.pow(2, -float(edge[2] + length[edge[1]] + forward_length[edge[0]]))] for
             edge in graph])
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
    tp = 0
    fp = 0
    for i in range(min(k, len(sorted_edges))):
        edge = sorted_edges[i]
        if [edge[0][0], edge[0][1]] in known_pathway or (not direction and [edge[0][1], edge[0][0]] in known_pathway):
            tp += 1
        else:
            fp += 1
        recalls.append(tp / len(known_pathway))
        precisions.append(tp / (tp + fp))
        if recalls[-1] > recall_bound:
            break
    return recalls, precisions


def compute_edge_probs(adj_list, r):
    edge_probs = []
    for node in node_to_id.values():
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
        for node in node_to_id.values():
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
        recalls, precisions = compute_recall_precision(sorted_edges, subpathway, k, direction, recall_bound)
        if WRITE_PRC:
            write_precision_recall(precisions, recalls, "results/" + dataset + "PR-" + method + ".txt")

        name = method
        plt.plot(recalls, precisions, color=color, label=name + " " + str(round(auc(recalls, precisions), 4)))

        print("AUPRC of " + method + ": " + str(auc(recalls, precisions)))

    return result, recalls, precisions


def compute_recall_precision_pathlinker(sorted_edges, known_pathway, k, direction, recall_bound=1.0):
    recalls = []
    precisions = []
    tp = 0
    fp = 0
    counter = 0
    while counter < min(k, len(sorted_edges)):
        ksp = sorted_edges[counter][2]
        while counter < len(sorted_edges) and sorted_edges[counter][2] == ksp:
            edge = sorted_edges[counter]
            if [edge[0], edge[1]] in known_pathway or (not direction and [edge[1], edge[0]] in known_pathway):
                tp += 1
            else:
                fp += 1
            counter += 1
        recalls.append(tp / len(known_pathway))
        precisions.append(tp / (tp + fp))
        if recalls[-1] > recall_bound:
            break
    return recalls, precisions


def read_pathlinker_output(path):
    edges = []
    with open(path, 'r') as f:
        for line in f:
            sp = line.split()
            if sp[0] != "#tail":
                edge = [node_to_id[sp[1]], node_to_id[sp[0]], int(sp[2])]
                if edge not in edges:
                    edges.append(edge)
    return edges


def add_pathlinker(path, color, k=1000000, direction=True):
    edges = read_pathlinker_output(path)
    recalls, precisions = compute_recall_precision_pathlinker(edges, subpathway, k, direction)

    name = "PathLinker"
    plt.plot(recalls, precisions, color=color, label=name + " " + str(round(auc(recalls, precisions), 4)))

    print("AUPRC of pathlinker: " + str(auc(recalls, precisions)))

    # computing the highest ranked edges
    return [[edges[i][0], edges[i][1]] for i in range(min(k, len(edges)))], len(edges), recalls, precisions


def run_edge_linker(dataset, color, k=1000000, recall_bound=1.0, direction=True):
    if RUN_EDGE_LINKER:
        new_graph = compute_new_graph("edge_linker")

        # computing the edge probabilities
        print("computing the edges for edge_linker")
        edges_probs = [[[edge[0], edge[1]], edge[2]] for edge in new_graph]
        sorted_edges = sorted(edges_probs, key=lambda x: -x[1])
        result = [[sorted_edges[i][0][0], sorted_edges[i][0][1]] for i in range(min(k, len(sorted_edges)))]
        if WRITE_EDGES_EDGE_LINKER:
            write_edges(result, "results/" + dataset + "edges-el.txt")
    else:
        result = read_edges("results/" + dataset + "edges-el.txt")
        sorted_edges = [[[r[0], r[1]], 0] for r in result]

    recalls = []
    precisions = []

    if PLOT_EDGES_PRC:
        # computing the precision and recall
        print("computing recall-precision curve for edge_linker")
        recalls, precisions = compute_recall_precision(sorted_edges, subpathway, k, direction, recall_bound)
        if WRITE_PRC:
            write_precision_recall(precisions, recalls, "results/" + dataset + "PR-el.txt")

        plt.plot(recalls, precisions, color=color, label="edgeLinker " + str(round(auc(recalls, precisions), 4)))

        print("AUPRC of edge_linker: " + str(auc(recalls, precisions)))

    return result, recalls, precisions


# reading the graph
_, node_to_id, graph, id_to_node, graph_map = reader.read_graph(input_graph, graph_type)
if WRITE_NODES_TO_ID_MAP:
    write_node_to_id("data/node_to_id_map.txt", node_to_id)

# creating the networkx graph
if RUN_ALGORITHMS or RUN_EDGE_LINKER:
    G.add_nodes_from(node_to_id)
    for e in graph:
        G.add_edge(int(e[1]), int(e[0]), weight=e[2])

    for e in graph:
        forward_G.add_edge(int(e[0]), int(e[1]), weight=e[2])

total_pathway_lengths = 0
total_pathway_node_lengths = 0

if JUST_CLEAN:
    clean_receptors_and_tfs(DATABASE, node_to_id, id_to_node, graph_map)
else:
    pathway_names = read_pathway_names(DATABASE)
    for pathway_name in pathway_names:
        pathlinker = "data/PathLinker_output/" + pathway_name + "k-2000-ranked-edges.txt"

        print("Pathway: " + pathway_name)

        # reading seeds and targets
        seeds, targets = read_source_and_destinations(DATABASE, pathway_name, node_to_id)

        # reading pathway
        if HAS_CLEANED_PATHWAY:
            subpathway = read_cleaned_pathway(DATABASE, pathway_name, node_to_id)
        else:
            pathway = read_pathway(DATABASE, pathway_name, node_to_id)
            subpathway = clean_pathway(DATABASE, pathway, pathway_name, graph, id_to_node)

        total_pathway_lengths += len(subpathway)

        # finding the distances of each node from the targets
        if RUN_ALGORITHMS or RUN_EDGE_LINKER:
            length = nx.multi_source_dijkstra_path_length(G, targets)
            forward_length = nx.multi_source_dijkstra_path_length(forward_G, seeds)

        # running the algorithms and get the pathways, true positives, and false positives
        if INCLUDE_PATHLINKER:
            pathway_pl, pl_edge_len, pl_recalls, pl_precisions = add_pathlinker(pathlinker, color=pallet[0],
                                                                                direction=False)
        else:
            pl_edge_len = len(graph)
            pl_recalls = []
            pl_precisions = []
            pathway_pl = []
        pathway_ours, our_recalls, our_precisions = run_algorithm(dataset=pathway_name, method="ours", color=pallet[1],
                                                                  alpha=5, c=0.25, k=pl_edge_len, direction=False)
        if INCLUDE_RWR:
            pathway_rwr, rwr_recalls, rwr_precisions = run_algorithm(dataset=pathway_name, method="rwr",
                                                                     color=pallet[2], k=pl_edge_len, direction=False)
        else:
            rwr_recalls = []
            rwr_precisions = []
            pathway_rwr = []
        if INCLUDE_EDGE_LINKER:
            pathway_el, el_recalls, el_precisions = run_edge_linker(dataset=pathway_name, color=pallet[3],
                                                                    k=pl_edge_len, direction=False)
        else:
            el_recalls = []
            el_precisions = []
            pathway_el = []
        if PLOT_EDGES_PRC:
            overall_recalls_ours.append(our_recalls)
            overall_precisions_ours.append(our_precisions)
            overall_recalls_pl.append(pl_recalls)
            overall_precisions_pl.append(pl_precisions)
            overall_recalls_rwr.append(rwr_recalls)
            overall_precisions_rwr.append(rwr_precisions)
            overall_recalls_el.append(el_recalls)
            overall_precisions_el.append(el_precisions)

            plt.legend()
            plt.title("recall-precision for " + pathway_name)
            plt.savefig("output/edge-PRC/" + pathway_name + ".png")
            plt.close()

        if COMPUTE_RTF:
            our_seeds, our_targets, rwr_seeds, rwr_targets, pl_seeds, pl_targets, el_seeds, el_targets = \
                plot_rtf_found(seeds, targets, pathway_ours, pathway_rwr, pathway_pl, pathway_el, pathway_name)
            overall_r_ours.append(our_seeds)
            overall_tf_ours.append(our_targets)
            overall_r_pl.append(pl_seeds)
            overall_tf_pl.append(pl_targets)
            overall_r_rwr.append(rwr_seeds)
            overall_tf_rwr.append(rwr_targets)
            overall_r_el.append(el_seeds)
            overall_tf_el.append(el_targets)

        # Node AU-PRC
        if PLOT_NODES_PRC:
            our_recalls, our_precisions, rwr_recalls, rwr_precisions, pl_recalls, pl_precisions, el_recalls, \
                el_precisions, node_len = plot_node_auprc(subpathway, pathway_ours, pathway_rwr, pathway_pl, pathway_el,
                                                          pathway_name)
            overall_node_recalls_ours.append(our_recalls)
            overall_node_precisions_ours.append(our_precisions)
            overall_node_recalls_pl.append(pl_recalls)
            overall_node_precisions_pl.append(pl_precisions)
            overall_node_recalls_rwr.append(rwr_recalls)
            overall_node_precisions_rwr.append(rwr_precisions)
            overall_node_recalls_el.append(el_recalls)
            overall_node_precisions_el.append(el_precisions)
            total_pathway_node_lengths += node_len

    if PLOT_EDGES_PRC:
        plot_total_prc(overall_recalls_ours, overall_precisions_ours, overall_recalls_rwr, overall_precisions_rwr,
                       overall_recalls_pl, overall_precisions_pl, overall_recalls_el, overall_precisions_el,
                       "overall-edge-PRC")

    if PLOT_NODES_PRC:
        plot_total_prc(overall_node_recalls_ours, overall_node_precisions_ours, overall_node_recalls_rwr,
                       overall_node_precisions_rwr, overall_node_recalls_pl, overall_node_precisions_pl,
                       overall_node_recalls_el, overall_node_precisions_el, "overall-node-PRC")

    if COMPUTE_RTF:
        plot_total_prc([[(i + 1) for i in range(len(r_ours))] for r_ours in overall_r_ours], overall_r_ours,
                       [[(i + 1) for i in range(len(r_rwr))] for r_rwr in overall_r_rwr], overall_r_rwr,
                       [[(i + 1) for i in range(len(r_pl))] for r_pl in overall_r_pl], overall_r_pl,
                       [[(i + 1) for i in range(len(r_el))] for r_el in overall_r_el], overall_r_el,
                       "overall-receptors")

        plot_total_prc([[(i + 1) for i in range(len(r_ours))] for r_ours in overall_tf_ours], overall_tf_ours,
                       [[(i + 1) for i in range(len(r_rwr))] for r_rwr in overall_tf_rwr], overall_tf_rwr,
                       [[(i + 1) for i in range(len(r_pl))] for r_pl in overall_tf_pl], overall_tf_pl,
                       [[(i + 1) for i in range(len(r_el))] for r_el in overall_tf_el], overall_tf_el, "overall-tfs")
