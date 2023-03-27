#!/usr/bin/env python
# -*- coding: utf-8 -*-

from utils import reader, normalizer
import numpy as np

from utils.reader import read_directed_graph, read_undirected_graph


class PyRWR:
    normalized = False

    def __init__(self):
        pass

    def read_graph(self, graph, weighted):
        '''
        Read a graph from the edge list at input_path

        inputs
            input_path : str
                path for the graph data
            graph_type : str
                type of graph {'directed', 'undirected', 'bipartite'}
        '''

        self.A = read_directed_graph(graph, weighted)
        self.m, self.n = self.A.shape
        self.node_ids = np.arange(0, self.n)
        self.normalize()

    def normalize(self):
        '''
        Perform row-normalization of the adjacency matrix
        '''
        if self.normalized is False:
            nA = normalizer.row_normalize(self.A)
            self.nAT = nA.T
            self.normalized = True
