import math
import random

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from sklearn.metrics import auc

from file_methods import *
from network_properties import plot_total_rtf, plot_node_auprc, plot_rtf_found
from pathway_files import read_cleaned_pathway, clean_pathway, read_pathway, read_source_and_destinations, \
    read_pathway_names, clean_receptors_and_tfs
from utils import reader
from pyrwr.ppr import PPR
import networkx as nx

G = nx.Graph()
forward_G = nx.Graph()

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

pallet = [colors['black'], colors['red'], colors['grey'], colors['lightgray'], colors['brown']]
# pallet = [colors['navy'], colors['orangered'], colors['dodgerblue'], colors['lightsteelblue']]

input_graph = "data/interactome.txt"
graph_type = "undirected"

''' 
These are some boolean variables that define what do we expect from the code 
'''
DATABASE = "NetPath"
# methods to run
RUN_ALGORITHMS = False
RUN_EDGE_LINKER = False
JUST_CLEAN = False
# output
PLOT_EDGES_PRC = True
PLOT_NODES_PRC = False
COMPUTE_RTF = False
PLOT_INDIVIDUAL_PATHWAYS = False
# writing methods
WRITE_EDGES = False
WRITE_EDGES_EDGE_LINKER = False
WRITE_PRC = False
WRITE_NODES_TO_ID_MAP = False
# which methods to include
INCLUDE_PATHLINKER = True
INCLUDE_RWR = True
INCLUDE_EDGE_LINKER = True
INCLUDE_GROWING_DAGS = False
# pathway
HAS_CLEANED_PATHWAY = True
# PRC
READ_PRC = True
READ_PRC_PL = True
# direction
DIRECTION = False
# subsampling
SUB_SAMPLE = False
# algorithm parameters
ALPHA = 5
C = 0.25

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

overall_tps_ours = []
overall_tps_pl = []
overall_tps_rwr = []
overall_tps_el = []


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


def compute_recall_precision(sorted_edges, known_pathway, k, direction, recall_bound=1.0, sub_samples=None):
    if sub_samples is None:
        sub_samples = []
    recalls = []
    precisions = []
    tp = 0
    fp = 0
    tps = []
    for i in range(min(k, len(sorted_edges))):
        edge = sorted_edges[i]
        if sub_samples:
            if edge[0][0] in sub_samples and edge[0][1] in sub_samples[edge[0][0]]:
                check = True
            else:
                check = False
        else:
            check = True
        if check:
            if [edge[0][0], edge[0][1]] in known_pathway or (
                    not direction and [edge[0][1], edge[0][0]] in known_pathway):
                tp += 1
                tps.append(1)
            else:
                fp += 1
                tps.append(0)
            recalls.append(tp / len(known_pathway))
            precisions.append(tp / (tp + fp))
            if recalls[-1] > recall_bound:
                break
    return recalls, precisions, tps


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


def run_algorithm(dataset, method, color, alpha=0.0, c=0.15, k=1000000, recall_bound=1.0, direction=True,
                  sub_samples=None):
    sorted_edges = []
    result = []
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
    elif not READ_PRC:
        result = read_edges("results/" + dataset + "edges-" + method + ".txt")
        sorted_edges = [[[r[0], r[1]], 0] for r in result]

    recalls = []
    precisions = []
    tps = []

    if PLOT_EDGES_PRC:
        # computing the precision and recall
        print("computing recall-precision curve for " + method)
        if READ_PRC:
            if SUB_SAMPLE:
                recalls, precisions = read_precision_recall("results/" + dataset + "PR-sub-" + method + ".txt")
            else:
                recalls, precisions = read_precision_recall("results/" + dataset + "PR-" + method + ".txt")
            recalls = [recalls[i] for i in range(min(k, len(recalls)))]
            precisions = [precisions[i] for i in range(min(k, len(precisions)))]
        else:
            recalls, precisions, tps = compute_recall_precision(sorted_edges, subpathway, k, direction, recall_bound,
                                                                sub_samples=sub_samples)
        if WRITE_PRC:
            if SUB_SAMPLE:
                write_precision_recall(precisions, recalls, "results/" + dataset + "PR-sub-" + method + ".txt")
            else:
                write_precision_recall(precisions, recalls, "results/" + dataset + "PR-" + method + ".txt")

        if PLOT_INDIVIDUAL_PATHWAYS:
            name = method
            plt.plot(recalls, precisions, color=color, label=name + " " + str(round(auc(recalls, precisions), 4)))

            print("AUPRC of " + method + ": " + str(auc(recalls, precisions)))

    return result, recalls, precisions, tps


def compute_recall_precision_pathlinker(sorted_edges, known_pathway, k, direction, recall_bound=1.0, sub_samples=None):
    if sub_samples is None:
        sub_samples = []
    recalls = []
    precisions = []
    tp = 0
    fp = 0
    counter = 0
    tps = []
    while counter < min(k, len(sorted_edges)):
        ksp = sorted_edges[counter][2]
        tps.append([0, 0])
        while counter < len(sorted_edges) and sorted_edges[counter][2] == ksp:
            edge = sorted_edges[counter]
            if sub_samples:
                if edge[0] in sub_samples and edge[1] in sub_samples[edge[0]]:
                    check = True
                else:
                    check = False
            else:
                check = True
            if check:
                if [edge[0], edge[1]] in known_pathway or (not direction and [edge[1], edge[0]] in known_pathway):
                    tp += 1
                    tps[-1][0] += 1
                else:
                    fp += 1
                    tps[-1][1] += 1
            counter += 1
        if tp + fp > 0:
            recalls.append(tp / len(known_pathway))
            precisions.append(tp / (tp + fp))
            if recalls[-1] > recall_bound:
                break
    return recalls, precisions, tps


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


def add_pathlinker(dataset, path, color, k=1000000, direction=True, name="PathLinker", sub_samples=None):
    edges = read_pathlinker_output(path)
    tps = []
    if READ_PRC_PL:
        if SUB_SAMPLE:
            recalls, precisions = read_precision_recall("results/" + dataset + "PR-sub-pl.txt")
        else:
            recalls, precisions = read_precision_recall("results/" + dataset + "PR-pl.txt")
    else:
        recalls, precisions, tps = compute_recall_precision_pathlinker(edges, subpathway, k, direction,
                                                                       sub_samples=sub_samples)

    if WRITE_PRC:
        if SUB_SAMPLE:
            write_precision_recall(precisions, recalls, "results/" + dataset + "PR-sub-pl.txt")
        else:
            write_precision_recall(precisions, recalls, "results/" + dataset + "PR-sub-pl.txt")

    if PLOT_INDIVIDUAL_PATHWAYS:
        plt.plot(recalls, precisions, color=color, label=name + " " + str(round(auc(recalls, precisions), 4)))

        print("AUPRC of pathlinker: " + str(auc(recalls, precisions)))

    # computing the highest ranked edges
    return [[edges[i][0], edges[i][1]] for i in range(min(k, len(edges)))], len(edges), recalls, precisions, tps


def run_edge_linker(dataset, color, k=1000000, recall_bound=1.0, direction=True, sub_samples=None):
    sorted_edges = []
    result = []
    if RUN_EDGE_LINKER:
        new_graph = compute_new_graph("edge_linker")

        # computing the edge probabilities
        print("computing the edges for edge_linker")
        edges_probs = [[[edge[0], edge[1]], edge[2]] for edge in new_graph]
        sorted_edges = sorted(edges_probs, key=lambda x: -x[1])
        result = [[sorted_edges[i][0][0], sorted_edges[i][0][1]] for i in range(min(k, len(sorted_edges)))]
        if WRITE_EDGES_EDGE_LINKER:
            write_edges(result, "results/" + dataset + "edges-el.txt")
    elif not READ_PRC:
        result = read_edges("results/" + dataset + "edges-el.txt")
        sorted_edges = [[[r[0], r[1]], 0] for r in result]

    recalls = []
    precisions = []
    tps = []

    if PLOT_EDGES_PRC:
        # computing the precision and recall
        print("computing recall-precision curve for edge_linker")
        if READ_PRC:
            if SUB_SAMPLE:
                recalls, precisions = read_precision_recall("results/" + dataset + "PR-sub-el.txt")
            else:
                recalls, precisions = read_precision_recall("results/" + dataset + "PR-el.txt")
            recalls = [recalls[i] for i in range(min(k, len(recalls)))]
            precisions = [precisions[i] for i in range(min(k, len(precisions)))]
        else:
            recalls, precisions, tps = compute_recall_precision(sorted_edges, subpathway, k, direction, recall_bound,
                                                                sub_samples=sub_samples)
        if WRITE_PRC:
            if SUB_SAMPLE:
                write_precision_recall(precisions, recalls, "results/" + dataset + "PR-sub-el.txt")
            else:
                write_precision_recall(precisions, recalls, "results/" + dataset + "PR-sub-el.txt")

        if PLOT_INDIVIDUAL_PATHWAYS:
            plt.plot(recalls, precisions, color=color, label="edgeLinker " + str(round(auc(recalls, precisions), 4)))

            print("AUPRC of edge_linker: " + str(auc(recalls, precisions)))

    return result, recalls, precisions, tps


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

indices = [i for i in range(len(graph))]

if JUST_CLEAN:
    clean_receptors_and_tfs(DATABASE, node_to_id, id_to_node, graph_map)
else:
    pathway_names = read_pathway_names(DATABASE, cleaned=True)
    for pathway_name in pathway_names:
        pathlinker = "data/PathLinker_output/" + pathway_name + "k-2000-ranked-edges.txt"
        growingDAGs = "data/GrowingDAGs_output/" + pathway_name + "-edges.txt"

        print("Pathway: " + pathway_name)

        # reading seeds and targets
        seeds, targets = read_source_and_destinations(DATABASE, pathway_name, node_to_id)

        # reading pathway
        if HAS_CLEANED_PATHWAY:
            subpathway = read_cleaned_pathway(DATABASE, pathway_name, node_to_id)
        else:
            pathway = read_pathway(DATABASE, pathway_name, node_to_id)
            subpathway = clean_pathway(DATABASE, pathway, pathway_name, graph_map)

        if SUB_SAMPLE:
            random.shuffle(indices)
            sub_graph = {}
            for j in range(50 * len(subpathway)):
                e = graph[indices[j]]
                if e[0] not in sub_graph:
                    sub_graph[e[0]] = []
                if e[1] not in sub_graph[e[0]]:
                    sub_graph[e[0]].append(e[1])
                    if graph_type == "undirected":
                        if e[1] not in sub_graph:
                            sub_graph[e[1]] = []
                        sub_graph[e[1]].append(e[0])
            for e in subpathway:
                if e[0] not in sub_graph:
                    sub_graph[e[0]] = []
                if e[1] not in sub_graph[e[0]]:
                    sub_graph[e[0]].append(e[1])
                    if graph_type == "undirected":
                        if e[1] not in sub_graph:
                            sub_graph[e[1]] = []
                        sub_graph[e[1]].append(e[0])
        else:
            sub_graph = []
        print("here")

        total_pathway_lengths += len(subpathway)

        # finding the distances of each node from the targets
        if RUN_ALGORITHMS or RUN_EDGE_LINKER:
            length = nx.multi_source_dijkstra_path_length(G, targets)
            forward_length = nx.multi_source_dijkstra_path_length(forward_G, seeds)

        edge_len = [100000]
        # running the algorithms and get the pathways, true positives, and false positives
        if INCLUDE_GROWING_DAGS:
            pathway_gd, gd_edge_len, gd_recalls, gd_precisions, gd_tps = \
                add_pathlinker(DATABASE, growingDAGs, color=pallet[4], direction=DIRECTION, name="GrowingDAGs",
                               sub_samples=sub_graph)
            edge_len[0] = gd_edge_len
        if INCLUDE_PATHLINKER:
            pathway_pl, pl_edge_len, pl_recalls, pl_precisions, pl_tps = \
                add_pathlinker(pathway_name, pathlinker, color=pallet[0], direction=DIRECTION, k=edge_len[0],
                               sub_samples=sub_graph)
            if not INCLUDE_GROWING_DAGS:
                edge_len[0] = pl_edge_len

        else:
            pl_recalls = []
            pl_precisions = []
            pathway_pl = []
            pl_tps = []

        pathway_ours, our_recalls, our_precisions, our_tps = \
            run_algorithm(dataset=pathway_name, method="ours", color=pallet[1], alpha=ALPHA, c=C, k=edge_len[0],
                          direction=DIRECTION, sub_samples=sub_graph)
        if INCLUDE_RWR:
            pathway_rwr, rwr_recalls, rwr_precisions, rwr_tps = \
                run_algorithm(dataset=pathway_name, method="rwr", color=pallet[2], k=edge_len[0], direction=DIRECTION,
                              sub_samples=sub_graph)
        else:
            rwr_recalls = []
            rwr_precisions = []
            pathway_rwr = []
            rwr_tps = []
        if INCLUDE_EDGE_LINKER:
            pathway_el, el_recalls, el_precisions, el_tps = \
                run_edge_linker(dataset=pathway_name, color=pallet[3], k=edge_len[0], direction=DIRECTION,
                                sub_samples=sub_graph)
        else:
            el_recalls = []
            el_precisions = []
            pathway_el = []
            el_tps = []

        if PLOT_EDGES_PRC:
            overall_recalls_ours.append(our_recalls)
            overall_precisions_ours.append(our_precisions)
            overall_recalls_pl.append(pl_recalls)
            overall_precisions_pl.append(pl_precisions)
            overall_recalls_rwr.append(rwr_recalls)
            overall_precisions_rwr.append(rwr_precisions)
            overall_recalls_el.append(el_recalls)
            overall_precisions_el.append(el_precisions)

            if PLOT_INDIVIDUAL_PATHWAYS:
                plt.legend()
                plt.title("recall-precision for " + pathway_name)
                plt.savefig("output/edge-PRC/" + pathway_name + ".png")
                plt.close()

        if COMPUTE_RTF:
            our_seeds, our_targets, rwr_seeds, rwr_targets, pl_seeds, pl_targets, el_seeds, el_targets = \
                plot_rtf_found(seeds, targets, pathway_ours, pathway_rwr, pathway_pl, pathway_el, pathway_name,
                               PLOT_INDIVIDUAL_PATHWAYS)
            overall_r_ours.append(our_seeds)
            overall_tf_ours.append(our_targets)
            overall_r_pl.append(pl_seeds)
            overall_tf_pl.append(pl_targets)
            overall_r_rwr.append(rwr_seeds)
            overall_tf_rwr.append(rwr_targets)
            overall_r_el.append(el_seeds)
            overall_tf_el.append(el_targets)

            overall_tps_ours.append(our_tps)
            overall_tps_pl.append(pl_tps)
            overall_tps_rwr.append(rwr_tps)
            overall_tps_el.append(el_tps)

        # Node AU-PRC
        if PLOT_NODES_PRC:
            our_recalls, our_precisions, rwr_recalls, rwr_precisions, pl_recalls, pl_precisions, el_recalls, \
                el_precisions, node_len = plot_node_auprc(subpathway, pathway_ours, pathway_rwr, pathway_pl, pathway_el,
                                                          pathway_name, PLOT_INDIVIDUAL_PATHWAYS, read_data=READ_PRC)
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
        plot_total_rtf(overall_recalls_ours, overall_precisions_ours, overall_recalls_rwr, overall_precisions_rwr,
                       overall_recalls_pl, overall_precisions_pl, overall_recalls_el, overall_precisions_el,
                       DATABASE + "-overall-edge-PRC", DATABASE, "Edge precision-recall curve",
                       1.05, 0.45, "recall", "precision")

    if PLOT_NODES_PRC:
        plot_total_rtf(overall_node_recalls_ours, overall_node_precisions_ours, overall_node_recalls_rwr,
                       overall_node_precisions_rwr, overall_node_recalls_pl, overall_node_precisions_pl,
                       overall_node_recalls_el, overall_node_precisions_el, DATABASE + "-overall-node-PRC", DATABASE,
                       "Node precision-recall curve", 1.05, 1, "recall",
                       "precision")

    if COMPUTE_RTF:
        plot_total_rtf([[(i + 1) for i in range(len(r_ours))] for r_ours in overall_r_ours], overall_r_ours,
                       [[(i + 1) for i in range(len(r_rwr))] for r_rwr in overall_r_rwr], overall_r_rwr,
                       [[(i + 1) for i in range(len(r_pl))] for r_pl in overall_r_pl], overall_r_pl,
                       [[(i + 1) for i in range(len(r_el))] for r_el in overall_r_el], overall_r_el,
                       DATABASE + "-overall-receptors", DATABASE, "percentage of receptors found", 1.05, 2000,
                       "Number of edges", "Percentage", False)

        plot_total_rtf([[(i + 1) for i in range(len(r_ours))] for r_ours in overall_tf_ours], overall_tf_ours,
                       [[(i + 1) for i in range(len(r_rwr))] for r_rwr in overall_tf_rwr], overall_tf_rwr,
                       [[(i + 1) for i in range(len(r_pl))] for r_pl in overall_tf_pl], overall_tf_pl,
                       [[(i + 1) for i in range(len(r_el))] for r_el in overall_tf_el], overall_tf_el,
                       DATABASE + "-overall-tfs", DATABASE, "percentage of transcription factors found", 1.05, 2000,
                       "Number of edges", "Percentage", False)

    print(total_pathway_lengths)
