path1 = "data/2015pathlinker-weighted.txt"
path2 = "data/interactome-weights.txt"
output_path = "data/interactome.txt"

X = []
graph = {}
graph_weighted = {}
count = [0]


def read_graph(path):
    with open(path, 'r') as file:
        for row in file:
            sp = row.split()
            if sp:
                if sp[0] not in graph:
                    graph[sp[0]] = []
                    graph_weighted[sp[0]] = []
                if sp[1] not in graph[sp[0]]:
                    X.append([sp[0], sp[1], float(sp[2])])
                    graph[sp[0]].append(sp[1])
                    graph_weighted[sp[0]].append([sp[1], float(sp[2]), len(X) - 1])
                else:
                    for e in graph_weighted[sp[0]]:
                        if e[0] == sp[1]:
                            if e[1] < float(sp[2]):
                                X[e[2]][2] = float(sp[2])
                            break


read_graph(path1)
read_graph(path2)
print(len(X))
with open(output_path, 'w') as f:
    for edge in X:
        f.write(edge[0] + "\t" + edge[1] + "\t" + str(edge[2]) + "\n")
