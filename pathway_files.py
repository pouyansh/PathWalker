# this method reads the name of the pathways
def read_pathway_names(database, signaling=True):
    names = []
    if database == "KEGG":
        if signaling:
            file = "data/KEGG/_signaling_pathway_list.txt"
        else:
            file = "data/KEGG/_pathway_list.txt"
        with open(file, 'r') as f:
            for line in f:
                if line:
                    sp = line.split()
                    if sp:
                        if not signaling:
                            sp2 = sp[0].split(":")
                            names.append(sp2[1])
                        else:
                            names.append(sp[0])
    if database == "NetPath":
        with open("data/NetPath_signaling_pathways.txt", 'r') as f:
            for line in f:
                if line:
                    sp = line.split()
                    names.append(sp[0])
    return names


def read_pathway(database, pathway_name, node_to_id):
    edges = []
    file = ""
    if database == "KEGG":
        file = "data/KEGG/" + pathway_name + "-expanded-edges.txt"
    if database == "NetPath":
        file = "data/NetPath/" + pathway_name + "-edges.txt"
    with open(file, 'r') as f:
        for line in f:
            sp = line.split()
            if sp[0] != "#tail" and sp[1] in node_to_id and sp[0] in node_to_id:
                if database == "KEGG":
                    edges.append([node_to_id[sp[0]], node_to_id[sp[1]]])
                else:
                    edges.append([node_to_id[sp[1]], node_to_id[sp[0]]])
            else:
                print(line)
    return edges


def clean_pathway(database, pathway_edges, pathway_name, graph, id_to_node):
    sub_edges = []
    counter = 0
    for p in pathway_edges:
        check = False
        for edge in graph:
            if edge[0] == p[0] and edge[1] == p[1]:
                sub_edges.append(p)
                check = True
                break
        if not check:
            # print(id_to_node[p[0]] + " " + id_to_node[p[1]])
            counter += 1
    with open("data/" + database + "/cleaned_pathways/" + pathway_name + "-edges.txt", 'w') as f:
        for p in sub_edges:
            f.write(id_to_node[p[0]] + " " + id_to_node[p[1]] + "\n")
    print(counter / max(1, len(pathway_edges)))
    return sub_edges


def read_cleaned_pathway(database, pathway_name, node_to_id):
    paths = []
    with open("data/" + database + "/cleaned_pathways/" + pathway_name + "-edges.txt", 'r') as f:
        for line in f:
            sp = line.split()
            if sp[0] != "#tail" and sp[1] in node_to_id and sp[0] in node_to_id:
                paths.append([node_to_id[sp[0]], node_to_id[sp[1]]])
    return paths


def read_source_and_destinations(database, pathway_name, node_to_id):
    sources = []
    destinations = []
    file = "data/" + database + "/" + pathway_name + "-node_to_id.txt"
    with open(file, 'r') as f:
        for line in f:
            sp = line.split()
            if sp[0] != "#node":
                if sp[1] == "tf":
                    destinations.append(node_to_id[sp[0]])
                if sp[1] == "receptor":
                    sources.append(node_to_id[sp[0]])
    return sources, destinations


def read_raw_source_and_destinations(database, node_to_id):
    if database == "KEGG":
        receptors = []
        tfs = []
        with open("data/KEGG/receptors/uniprot-target-list.txt", 'r') as f:
            for line in f:
                if line and line in node_to_id:
                    receptors.append(node_to_id[line])
        with open("data/KEGG/transcription-factors/TFcheckpoint/all-tfs.txt", 'r') as f:
            for line in f:
                if line and line in node_to_id:
                    tfs.append(node_to_id[line])
        return receptors, tfs
    return [], []


def compute_intersection(database, pathway_name, receptors, tfs, pathway_edges, id_to_node):
    sub_receptors = []
    sub_tfs = []
    pathway_nodes = []
    for edge in pathway_edges:
        if edge[0] not in pathway_nodes:
            pathway_nodes.append(edge[0])
        if edge[1] not in pathway_nodes:
            pathway_nodes.append(edge[1])
    for receptor in receptors:
        if receptor in pathway_nodes:
            sub_receptors.append(receptor)
    for tf in tfs:
        if tf in pathway_nodes:
            sub_tfs.append(tf)

    with open("data/" + database + "/" + pathway_name + "-nodes.txt", 'w') as f:
        f.write("#node\tnode_symbol\n")
        for receptor in sub_receptors:
            f.write(id_to_node[receptor] + "\treceptor\n")
        for tf in sub_tfs:
            f.write(id_to_node[tf] + "\ttf\n")
        for node in pathway_nodes:
            if node not in sub_receptors and node not in sub_tfs:
                f.write(id_to_node[node] + "\tnone\n")


def clean_receptors_and_tfs(database, graph, node_to_id, id_to_node):
    pathway_names = read_pathway_names(database)

    receptors, tfs = read_raw_source_and_destinations(database, node_to_id)

    for pathway_name in pathway_names:
        print(pathway_name)
        pathway_edges = read_pathway(database, pathway_name, node_to_id)
        sub_edges = clean_pathway(database, pathway_edges, pathway_name, graph, id_to_node)
        compute_intersection(database, pathway_name, receptors, tfs, sub_edges, id_to_node)
