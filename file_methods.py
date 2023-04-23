def write_edges(edges, filename):
    with open(filename, 'w') as f:
        for edge in edges:
            f.write(str(edge[0]) + " " + str(edge[1]) + "\n")


def read_edges(filename):
    edges = []
    with open(filename, 'r') as f:
        for row in f:
            sp = row.split()
            edges.append([int(float(sp[0])), int(float(sp[1]))])
    return edges


def write_precision_recall(precisions, recalls, filename):
    with open(filename, 'w') as f:
        for i in range(len(precisions)):
            f.write(str(precisions[i]) + " " + str(recalls[i]) + "\n")


def read_precision_recall(filename):
    precisions = []
    recalls = []
    with open(filename, 'r') as f:
        for row in f:
            sp = row.split()
            precisions.append(float(sp[0]))
            recalls.append(float(sp[1]))
    return precisions, recalls


def write_node_to_id(filename, nodes):
    with open(filename, 'w') as f:
        for node in nodes.keys():
            f.write(node + " " + str(nodes[node]) + "\n")


def read_nodes(filename):
    nodes = {}
    with open(filename, 'r') as f:
        for line in f:
            sp = line.split()
            nodes[sp[0]] = int(sp[1])
    return nodes
