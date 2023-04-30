# this method reads the name of the pathways
def read_pathway_names(database, cleaned=False):
    names = []
    if database == "KEGG":
        if cleaned:
            file = "data/KEGG/cleaned_pathways/_signaling_pathway_list.txt"
        else:
            file = "data/KEGG/_signaling_pathway_list.txt"
        with open(file, 'r') as f:
            for line in f:
                if line:
                    sp = line.split()
                    if sp:
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
    edges_original = []
    file = ""
    count = 0
    total = 0
    if database == "KEGG":
        file = "data/KEGG/" + pathway_name + "-expanded-edges.txt"
    if database == "NetPath":
        file = "data/NetPath/" + pathway_name + "-edges.txt"
    with open(file, 'r') as f:
        for line in f:
            total += 1
            sp = line.split()
            if sp[0] != "#tail" and sp[1] in node_to_id and sp[0] in node_to_id:
                if database == "KEGG":
                    edges.append([node_to_id[sp[0]], node_to_id[sp[1]]])
                    edges_original.append([sp[0], sp[1]])
                else:
                    edges.append([node_to_id[sp[1]], node_to_id[sp[0]]])
                    edges_original.append([sp[1], sp[0]])
            else:
                count += 1
                # print(line)
    return edges, edges_original, count/total


def clean_pathway(database, pathway_edges, pathway_name, graph_map):
    sub_edges = []
    counter = 0
    for p in pathway_edges:
        check = False
        if p[0] in graph_map and p[1] in graph_map[p[0]]:
            sub_edges.append(p)
            check = True
        if not check:
            # print(id_to_node[p[0]] + " " + id_to_node[p[1]])
            counter += 1
    with open("data/" + database + "/cleaned_pathways/" + pathway_name + "-edges.txt", 'w') as f:
        for p in sub_edges:
            f.write(p[0] + " " + p[1] + "\n")
    return sub_edges, counter / max(1, len(pathway_edges))


def read_cleaned_pathway(database, pathway_name, node_to_id):
    paths = []
    file = "data/" + database + "/cleaned_pathways/" + pathway_name + "-edges.txt"
    with open(file, 'r') as f:
        for line in f:
            sp = line.split()
            if sp[0] != "#tail" and sp[1] in node_to_id and sp[0] in node_to_id:
                paths.append([node_to_id[sp[0]], node_to_id[sp[1]]])
    return paths


def read_source_and_destinations(database, pathway_name, node_to_id):
    sources = []
    destinations = []
    file = "data/" + database + "/cleaned_pathways/" + pathway_name + "-nodes.txt"
    with open(file, 'r') as f:
        for line in f:
            sp = line.split()
            if sp[0] != "#node":
                if sp[1] == "tf":
                    destinations.append(node_to_id[sp[0]])
                if sp[1] == "receptor":
                    sources.append(node_to_id[sp[0]])
    return sources, destinations


def read_raw_source_and_destinations(node_to_id):
    receptors = []
    tfs = []
    with open("data/KEGG/receptors/uniprot-target-list.txt", 'r') as f:
        for line in f:
            if line.split() and line.split()[0] in node_to_id:
                receptors.append(node_to_id[line.split()[0]])
    with open("data/KEGG/transcription-factors/TFcheckpoint/db-tfs.txt", 'r') as f:
        for line in f:
            if line.split() and line.split()[0] in node_to_id:
                tfs.append(node_to_id[line.split()[0]])
    return receptors, tfs


def compute_intersection(database, pathway_name, receptors, tfs, pathway_edges, id_to_node, node_to_id):
    sub_receptors = []
    sub_tfs = []
    pathway_nodes = []
    for edge in pathway_edges:
        if node_to_id[edge[0]] not in pathway_nodes:
            pathway_nodes.append(node_to_id[edge[0]])
        if node_to_id[edge[1]] not in pathway_nodes:
            pathway_nodes.append(node_to_id[edge[1]])
    for receptor in receptors:
        if receptor in pathway_nodes:
            sub_receptors.append(receptor)
    for tf in tfs:
        if tf in pathway_nodes:
            sub_tfs.append(tf)
        if tf in receptors:
            print(tf)

    with open("data/" + database + "/cleaned_pathways/" + pathway_name + "-nodes.txt", 'w') as f:
        f.write("#node\tnode_symbol\n")
        for receptor in sub_receptors:
            f.write(id_to_node[receptor] + "\treceptor\n")
        for tf in sub_tfs:
            f.write(id_to_node[tf] + "\ttf\n")
        for node in pathway_nodes:
            if node not in sub_receptors and node not in sub_tfs:
                f.write(id_to_node[node] + "\tnone\n")
    if len(sub_receptors) > 0 and len(sub_tfs) > 0:
        return True
    return False


def clean_receptors_and_tfs(database, node_to_id, id_to_node, graph_map):
    pathway_names = read_pathway_names(database)

    receptors, tfs = read_raw_source_and_destinations(node_to_id)

    with open("data/" + database + "/cleaned_pathways/_signaling_pathway_list.txt", 'w') as f:
        for pathway_name in pathway_names:
            print(pathway_name)
            pathway_edges, pathway_edges_original, percentage1 = read_pathway(database, pathway_name, node_to_id)
            sub_edges, percentage2 = clean_pathway(database, pathway_edges_original, pathway_name, graph_map)
            if compute_intersection(database, pathway_name, receptors, tfs, sub_edges, id_to_node, node_to_id):
                print("overall", percentage1 + (1-percentage1) * percentage2)
                if percentage1 + (1-percentage1) * percentage2 < 0.4:
                    f.write(pathway_name + "\n")
