import math
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from sklearn.metrics import auc

from file_methods import write_edges, write_precision_recall
from network_properties import plot_rtf_found, plot_node_auprc, plot_total_prc
from utils import reader
from pyrwr.ppr import PPR
import networkx as nx

G = nx.Graph()

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

datas = ["Alpha6Beta4Integrin", "AndrogenReceptor", "BCR", "BDNF", "CRH", "EGFR1", "FSH", "Hedgehog", "IL1",
         "IL2", "IL3", "IL4", "IL5", "IL6", "IL9", "IL-7", "KitReceptor", "Leptin", "Notch", "Oncostatin_M",
         "Prolactin", "RANKL", "TCR", "TGF_beta_Receptor", "TNFalpha", "TSH", "TSLP", "TWEAK", "Wnt"]

input_graph = "data/interactome-weights.txt"
graph_type = "directed"

overall_tps_ours = []
overall_fps_ours = []
overall_tps_pl = []
overall_fps_pl = []
overall_tps_rwr = []
overall_fps_rwr = []

overall_node_tps_ours = []
overall_node_fps_ours = []
overall_node_tps_pl = []
overall_node_fps_pl = []
overall_node_tps_rwr = []
overall_node_fps_rwr = []


# This method computes the edge weights based on the given method
def compute_new_graph(method, alpha=0.0):
    if method == "m5":
        return np.array(
            [[int(edge[0]), int(edge[1]), 1 / ((alpha + length[edge[0]]) * (alpha + length[edge[1]]))] for edge in
             graph])
    if method == "m2":
        return np.array(
            [[int(edge[0]), int(edge[1]), math.pow(2, -alpha * float(length[edge[1]]))] for edge in graph])
    if method == "m3":
        return np.array(
            [[int(edge[0]), int(edge[1]), 1 / (length[edge[0]] * length[edge[1]] + alpha)] for edge in graph])
    if method == "m4":
        return np.array(
            [[int(edge[0]), int(edge[1]), 1 / (length[edge[1]] + alpha)] for edge in graph])
    # This is our main weighting method
    if method == "ours":
        return np.array(
            [[int(edge[0]), int(edge[1]), math.pow(2, -alpha * float(edge[2] + length[edge[1]]))] for edge in graph])
    if method == "rwr":
        return np.array([[int(edge[0]), int(edge[1]), math.pow(2, -float(edge[2]))] for edge in graph])


def compute_recall_precision(sorted_edges, known_pathway, k, recall_bound=1.0):
    recalls = []
    precisions = []
    tps = []
    fps = []
    tp = 0
    fp = 0
    for i in range(min(k, len(sorted_edges))):
        edge = sorted_edges[i]
        if [edge[0][0], edge[0][1]] in known_pathway:
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


def run_algorithm(dataset, method, color, alpha=0.0, c=0.15, k=1000000, recall_bound=1.0):
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
    # write_edges(result, "results/" + dataset + "edges-" + method + ".txt")

    # computing the precision and recall
    print("computing recall-precision curve for " + method)
    recalls, precisions, tps, fps = compute_recall_precision(sorted_edges, subpathway, k, recall_bound)
    # write_precision_recall(precisions, recalls, "results/" + dataset + "PR-" + method + ".txt")

    name = method
    plt.plot(recalls, precisions, color=color, label=name + " " + str(round(auc(recalls, precisions), 4)))

    print("AUPRC of " + method + ": " + str(auc(recalls, precisions)))

    # computing the highest ranked edges
    return result, tps, fps


def read_pathway(path):
    paths = []
    with open(path, 'r') as f:
        for line in f:
            splitted = line.split()
            if splitted[0] != "#tail" and splitted[1] in nodes and splitted[0] in nodes:
                paths.append([nodes[splitted[1]], nodes[splitted[0]]])
    return paths


def read_source_and_destinations(path):
    sources = []
    destinations = []
    with open(path, 'r') as f:
        for line in f:
            splitted = line.split()
            if splitted[0] != "#node":
                if splitted[1] == "tf":
                    destinations.append(nodes[splitted[0]])
                if splitted[1] == "receptor":
                    sources.append(nodes[splitted[0]])
    return sources, destinations


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


def add_pathlinker(path, color, k=1000000):
    edges = read_pathlinker_output(path)
    recalls, precisions, tps, fps = compute_recall_precision(edges, subpathway, k)

    name = "PathLinker"
    plt.plot(recalls, precisions, color=color, label=name + " " + str(round(auc(recalls, precisions), 4)))

    print("AUPRC of pathlinker: " + str(auc(recalls, precisions)))

    # computing the highest ranked edges
    return [[edges[i][0][0], edges[i][0][1]] for i in range(min(k, len(edges)))], len(edges), tps, fps


# reading the graph
_, nodes, graph = reader.read_graph(input_graph, graph_type)

# creating the networkx graph
G.add_nodes_from(nodes)
for e in graph:
    G.add_edge(int(e[1]), int(e[0]), weight=e[2])

total_pathway_lengths = 0
total_pathway_node_lengths = 0
for data in datas:
    rtf_path = "data/NetPath/" + data + "-nodes.txt"
    pathway_path = "data/NetPath/" + data + "-edges.txt"
    pathlinker = "data/PathLinker_output/" + data + "k-2000-ranked-edges.txt"

    # reading seeds and targets
    seeds, targets = read_source_and_destinations(rtf_path)

    # reading pathway
    pathway = read_pathway(pathway_path)

    # removing the set of edges in the pathway that are not in the interactome
    subpathway = []
    for p in pathway:
        for e in graph:
            if e[0] == p[0] and e[1] == p[1]:
                subpathway.append(p)

    total_pathway_lengths += len(subpathway)

    # finding the distances of each node from the targets
    length = nx.multi_source_dijkstra_path_length(G, targets)

    # running the algorithms and get the pathways, true positives, and false positives
    pathway_pathlinker, pl_edge_len, pl_tps, pl_fps = add_pathlinker(pathlinker, color=colors["black"])
    pathway_ours, our_tps, our_fps = run_algorithm(dataset=data, method="ours", color=colors["deepskyblue"], alpha=5,
                                                   c=0.25, k=50000)
    pathway_rwr, rwr_tps, rwr_fps = run_algorithm(dataset=data, method="rwr", color=colors["silver"], k=50000)
    overall_tps_ours.append(our_tps)
    overall_fps_ours.append(our_fps)
    overall_tps_pl.append(pl_tps)
    overall_fps_pl.append(pl_fps)
    overall_tps_rwr.append(rwr_tps)
    overall_fps_rwr.append(rwr_fps)

    plt.legend()
    plt.title("recall-precision for " + data)

    # plt.savefig("output/" + data + "/PRC.png")
    plt.savefig("output/edge-PRC/" + data + ".png")
    plt.close()

    # plot_rtf_found(seeds, targets, pathway_ours, pathway_rwr, pathway_pathlinker, data)
    our_tps, our_fps, rwr_tps, rwr_fps, pl_tps, pl_fps, node_len = plot_node_auprc(pathway, pathway_ours, pathway_rwr,
                                                                                   pathway_pathlinker, data)
    overall_node_tps_ours.append(our_tps)
    overall_node_fps_ours.append(our_fps)
    overall_node_tps_pl.append(pl_tps)
    overall_node_fps_pl.append(pl_fps)
    overall_node_tps_rwr.append(rwr_tps)
    overall_node_fps_rwr.append(rwr_fps)
    total_pathway_node_lengths += node_len

plot_total_prc(overall_tps_ours, overall_fps_ours, overall_tps_rwr, overall_fps_rwr,
               overall_tps_pl, overall_fps_pl, total_pathway_lengths, "overall-edge-PRC")

plot_total_prc(overall_node_tps_ours, overall_node_fps_ours, overall_node_tps_rwr, overall_node_fps_rwr,
               overall_node_tps_pl, overall_node_fps_pl, total_pathway_node_lengths, "overall-node-PRC")
