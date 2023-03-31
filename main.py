import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

from utils import reader
from pyrwr.ppr import PPR
import networkx as nx

G = nx.Graph()

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

data = "TGF_beta_Receptor"

input_graph = "data/interactome-weights.txt"
graph_type = "directed"
rtf_path = "data/NetPath/" + data + "-nodes.txt"
pathway_path = "data/NetPath/" + data + "-edges.txt"

if data == "Wnt":
    pathlinker = "data/Wnt/test-files/Wnt-weighted-paths-pathlinker.txt"
elif data == "TNFalpha":
    pathlinker = "data/TNFalpha/test-files/TNFalpha-weighted-paths-pathlinker.txt"
else:
    pathlinker = "data/TGF_beta/test-files/TGF_beta-weighted-paths-pathlinker.txt"


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
    if method == "ours":
        return np.array(
            [[int(edge[0]), int(edge[1]), math.pow(2, -alpha * float(edge[2] + length[edge[1]]))] for edge in graph])
    if method == "rwr":
        return np.array([[int(edge[0]), int(edge[1]), math.pow(2, -float(edge[2]))] for edge in graph])


def compute_recall_precision(sorted_edges, known_pathway, recall_bound=1.0):
    recalls = []
    precisions = []
    tp = 0
    fp = 0
    for edge in sorted_edges:
        if [edge[0][0], edge[0][1]] in known_pathway:
            tp += 1
        else:
            fp += 1
        recalls.append(tp / len(known_pathway))
        precisions.append(tp / (tp + fp))
        if recalls[-1] >= recall_bound:
            break
    return recalls, precisions


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


def run_algorithm(method, color, alpha=0.0, c=0.15):
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

    # computing the precision and recall
    print("computing recall-precision curve for " + method)
    recalls, precisions = compute_recall_precision(sorted_edges, subpathway)

    name = method
    plt.plot(recalls, precisions, color=color, label=name)


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
            splitted = line.split()
            sp = splitted[2].split('|')
            for i in range(len(sp) - 1):
                if sp[i] in nodes and sp[i+1] in nodes:
                    edge = [nodes[sp[i]], nodes[sp[i + 1]]]
                    if edge not in edges:
                        edges.append(edge)
                else:
                    print(sp)
    return [[edge] for edge in edges]


def add_pathlinker(path, color):
    edges = read_pathlinker_output(path)
    recalls, precisions = compute_recall_precision(edges, subpathway)

    name = "PathLinker"
    plt.plot(recalls, precisions, color=color, label=name)


# reading the graph
_, nodes, graph = reader.read_graph(input_graph, graph_type)

# reading seeds and targets
seeds, targets = read_source_and_destinations(rtf_path)

# reading pathway
pathway = read_pathway(pathway_path)

subpathway = []
for p in pathway:
    for e in graph:
        if e[0] == p[0] and e[1] == p[1]:
            subpathway.append(p)

G.add_nodes_from(nodes)
for e in graph:
    G.add_edge(int(e[1]), int(e[0]), weight=e[2])

# finding the distances of each node from the targets
length = nx.multi_source_dijkstra_path_length(G, targets)
max_length = max(length.values())

# run_algorithm(method="m1", color='b', alpha=0.01)
# run_algorithm(method="m1", color='r', alpha=0.1)
# run_algorithm(method="m2", color='k', alpha=1)
# run_algorithm(method="m2", color='c', alpha=10)
run_algorithm(method="ours", color=colors["deepskyblue"], alpha=5, c=0.25)
# run_algorithm(method="m3", color='g', alpha=0.001)
# run_algorithm(method="m4", color='c', alpha=0.001)
run_algorithm(method="rwr", color=colors["silver"])
add_pathlinker(pathlinker, color=colors["black"])
plt.legend()
plt.title("recall-precision for " + data)

plt.savefig("results-" + data + ".png")
plt.close()
