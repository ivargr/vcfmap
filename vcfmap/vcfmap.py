import logging
import numpy as np


class VcfMap:
    # More simple VcfMap that assumes no overlapping variants, and never more than two edges out for any node
    def __init__(self, from_nodes_to_haplotypes, from_nodes_to_to_nodes, from_nodes_to_n_haplotypes, haplotypes, n_haplotypes,
                 from_nodes_to_missing_haplotypes, from_nodes_to_n_missing_haplotypes, missing_haplotypes, graph_min_node=0):
        self.graph_min_node = 0
        self.from_nodes_to_haplotypes = from_nodes_to_haplotypes
        self.from_nodes_to_to_nodes = from_nodes_to_to_nodes
        self.from_nodes_to_n_haplotypes = from_nodes_to_n_haplotypes
        self.haplotypes = haplotypes
        self.n_haplotypes = n_haplotypes
        self.from_nodes_to_missing_haplotypes = from_nodes_to_missing_haplotypes
        self.from_nodes_to_n_missing_haplotypes = from_nodes_to_n_missing_haplotypes
        self.missing_haplotypes = missing_haplotypes
        self.missing_haplotypes_set = set(missing_haplotypes)

        self.possible_haplotypes = set(range(0, n_haplotypes))

    def get_haplotypes_on_edge(self, from_node, to_node):
        from_node = from_node - self.graph_min_node
        to_node = to_node - self.graph_min_node

        index = self.from_nodes_to_haplotypes[from_node]
        n_haplotypes = self.from_nodes_to_n_haplotypes[from_node]
        if n_haplotypes == 0:
            return None
        haplotypes = self.haplotypes[index:index+n_haplotypes]

        missing_index = self.from_nodes_to_n_missing_haplotypes[from_node]
        n_missing = self.from_nodes_to_n_missing_haplotypes[from_node]
        missing = set(self.missing_haplotypes[missing_index:missing_index+n_missing])

        if self.from_nodes_to_to_nodes[from_node] == to_node + self.graph_min_node:
            # We have a match in the index
            return set(haplotypes) - missing
        else:
            # No match
            return self.possible_haplotypes - set(haplotypes) - missing

    def to_file(self, file_name):
        np.savez(file_name,
                 from_nodes_to_haplotypes=self.from_nodes_to_haplotypes,
                 from_nodes_to_to_nodes=self.from_nodes_to_to_nodes,
                 from_nodes_to_n_haplotypes=self.from_nodes_to_n_haplotypes,
                 haplotypes=self.haplotypes,
                 n_haplotypes=len(self.possible_haplotypes),
                 from_nodes_to_missing_haplotypes=self.from_nodes_to_missing_haplotypes,
                 from_nodes_to_n_missing_haplotypes=self.from_nodes_to_n_missing_haplotypes,
                 missing_haplotypes=self.missing_haplotypes)

    @classmethod
    def from_file(cls, file_name):
        data = np.load(file_name)
        return cls(data["from_nodes_to_haplotypes"],
                   data["from_nodes_to_to_nodes"],
                   data["from_nodes_to_n_haplotypes"],
                   data["haplotypes"],
                   data["n_haplotypes"],
                   data["from_nodes_to_missing_haplotypes"],
                   data["from_nodes_to_n_missing_haplotypes"],
                   data["missing_haplotypes"])

    def get_n_haplotypes_on_node(self, node):
        # Returns all haplotypes that passes that node
        node = node - self.graph_min_node
        n_missing = self.from_nodes_to_n_missing_haplotypes[node]
        return self.n_haplotypes - n_missing

    def allele_frequency(self, from_node, to_node):
        haplotypes = self.get_haplotypes_on_edge(from_node, to_node)
        if haplotypes is None:
            # This means an edge is not a variant
            return 1.0
        return len(haplotypes) / self.get_n_haplotypes_on_node(from_node)

    def interval_allele_frequencies(self, interval):
        # Returns list of allele frequencies for all edges in interval
        rps = interval.region_paths
        return [self.allele_frequency(from_node, to_node) \
                for from_node, to_node in zip(rps[0:-1], rps[1:])]

    def haplotypes_consistent_with_whole_interval(self, interval):
        haplotypes_on_edges = []
        rps = interval.region_paths
        for from_node, to_node in zip(rps[0:-1], rps[1:]):
            haplotypes = self.get_haplotypes_on_edge(from_node, to_node)
            if haplotypes is None:
                continue
            #assert len(haplotypes) > 0, "Edge %d-%d has no haplotypes" % (from_node, to_node)
            haplotypes_on_edges.append(haplotypes)

        if len(haplotypes_on_edges) == 0:
            # No edges, means interval has all haplotypes
            return self.possible_haplotypes

        common = haplotypes_on_edges[0]
        for h in haplotypes_on_edges[1:]:
            common = common.intersection(h)

        return common



class VcfMapComplex:
    def __init__(self, edge_ids, edge_index_to_haplotypes, edge_index_to_n_haplotypes, haplotypes):
        self._edge_ids = edge_ids
        self._edge_index_to_haplotypes = edge_index_to_haplotypes
        self._edge_index_to_n_haplotypes = edge_index_to_n_haplotypes
        self._haplotypes = haplotypes

    def get_haplotypes_on_edge_id(self, edge_id):
        edge_index = np.searchsorted(self._edge_ids, edge_id)
        if self._edge_ids[edge_index] != edge_id:
            # Edge is not in index
            return None

        return self._haplotypes[self._edge_index_to_haplotypes[edge_index]:
                                self._edge_index_to_haplotypes[edge_index] + self._edge_index_to_haplotypes[edge_index]]

    def get_haplotypes_on_edge(self, from_node, to_node):
        edge_id = from_node * 100000000 + to_node
        return self.get_haplotypes_on_edge_id(edge_id)

    def to_file(self, file_name):
        np.savez(file_name,
                 edge_ids=self._edge_ids,
                 edge_index_to_haplotypes=self._edge_index_to_haplotypes,
                 edge_index_to_n_haplotypes=self._edge_index_to_n_haplotypes,
                 haplotypes=self._haplotypes
                 )

    @classmethod
    def from_file(cls, file_name):
        data = np.load(file_name)
        return cls(data["edge_ids"],
                   data["edge_index_to_haplotypes"],
                   data["edge_index_to_n_haplotypes"],
                   data["haplotypes"])

