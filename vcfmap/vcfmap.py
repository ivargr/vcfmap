import logging
import numpy as np


class VcfMap:
    # More simple VcfMap that assumes no overlapping variants, and never more than two edges out for any node
    def __init__(self, from_nodes_to_haplotypes, from_nodes_to_to_nodes, from_nodes_to_n_haplotypes, haplotypes, n_haplotypes):
        self.from_nodes_to_haplotypes = from_nodes_to_haplotypes
        self.from_nodes_to_to_nodes = from_nodes_to_to_nodes
        self.from_nodes_to_n_haplotypes = from_nodes_to_n_haplotypes
        self.haplotypes = haplotypes
        self.possible_haplotypes = set(range(0, n_haplotypes))

    def get_haplotypes_on_edge(self, from_node, to_node):
        index = self.from_nodes_to_haplotypes[from_node]
        n_haplotypes = self.from_nodes_to_n_haplotypes[from_node]
        haplotypes = self.haplotypes[index:index+n_haplotypes]

        if self.from_nodes_to_to_nodes[from_node] == to_node:
            # We have a match in the index
            return set(haplotypes)
        else:
            # No match
            return self.possible_haplotypes - set(haplotypes)

    def to_file(self, file_name):
        np.savez(file_name,
                 from_nodes_to_haplotypes=self.from_nodes_to_haplotypes,
                 from_nodes_to_to_nodes=self.from_nodes_to_to_nodes,
                 from_nodes_to_n_haplotypes=self.from_nodes_to_n_haplotypes,
                 haplotypes=self.haplotypes,
                 n_haplotypes=len(self.possible_haplotypes))

    @classmethod
    def from_file(cls, file_name):
        data = np.load(file_name)
        return cls(data["from_nodes_to_haplotypes"],
                   data["from_nodes_to_to_nodes"],
                   data["from_nodes_to_n_haplotypes"],
                   data["haplotypes"],
                   data["n_haplotypes"])


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

