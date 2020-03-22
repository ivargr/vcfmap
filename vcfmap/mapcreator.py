import logging
from offsetbasedgraph import Graph, NumpyIndexedInterval, SequenceGraph
import pickle
import gzip
from collections import defaultdict
from .vcfmap import VcfMap
import numpy as np
from array import array


def get_variant_type(vcf_line):

    l = vcf_line.split()
    if len(l[3]) == len(l[4]):
        return "SNP"
    elif "VT=SNP" in vcf_line:
        return "SNP"
    elif "VT=INDEL" in vcf_line:
        if len(l[3]) > len(l[4]):
            return "DELETION"
        else:
            return "INSERTION"
    else:
        raise Exception("Unsupported variant type on line %s" % vcf_line)


def get_vcf_sample_names(file_name):
    with gzip.open(file_name) as f:
        for line in f:
            line = line.decode("utf-8")
            if line.startswith("#CHROM"):
                samples = list(line.split())[9:]
                return samples


class MapCreator:
    def __init__(self,  graph, sequence_graph, reference_path, linear_reference_nodes, vcf_file_name):
        logging.info("Initing mapcreator")
        self.graph = graph
        self.sequence_graph = sequence_graph
        self.reference_path = reference_path
        self.linear_reference_nodes = linear_reference_nodes
        self.vcf_file_name = vcf_file_name
        self.n_deletions = 0
        self.n_insertions = 0
        self.n_substitutions = 0
        self.graph_min_node = graph.min_node

        self.sample_names = get_vcf_sample_names(self.vcf_file_name)
        logging.info("There are %d samples in vcf file" % len(self.sample_names))
        self.haplotype_to_edges = defaultdict(list)
        self.edges_to_haplotypes = defaultdict(list)

        self._from_nodes_to_haplotypes = np.zeros(len(graph.blocks), dtype=np.uint32)
        self._from_nodes_to_to_nodes = np.zeros(len(graph.blocks), dtype=np.uint32)
        self._from_nodes_to_n_haplotypes = np.zeros(len(graph.blocks), dtype=np.uint16)
        self._haplotypes = array("I") # np.array([])
        self._n_haplotypes = 0

        self._from_nodes_to_missing_haplotypes = np.zeros(len(graph.blocks), dtype=np.uint32)
        self._from_nodes_to_n_missing_haplotypes = np.zeros(len(graph.blocks), dtype=np.uint16)
        self._missing_haplotypes = array("I")

    def _make_vcfmap(self):
        self.vcfmap = VcfMap(self._from_nodes_to_haplotypes,
                               self._from_nodes_to_to_nodes,
                               self._from_nodes_to_n_haplotypes,
                               np.array(self._haplotypes, dtype=np.uint16),
                               self._n_haplotypes,
                             self._from_nodes_to_missing_haplotypes,
                             self._from_nodes_to_n_missing_haplotypes,
                             self._missing_haplotypes)

    def _store_processed_variant(self, line, edge):
        #print("adding edge %s" % (str(edge)))
        from_node = edge[0]
        to_node = edge[1]
        self._from_nodes_to_to_nodes[from_node-self.graph_min_node] = to_node
        self._from_nodes_to_haplotypes[from_node-self.graph_min_node] = len(self._haplotypes)
        self._from_nodes_to_missing_haplotypes[from_node-self.graph_min_node] = len(self._missing_haplotypes)

        n_haplotypes = 0
        n_missing_haplotypes = 0

        for i, genotype in enumerate(line[9:]):
            genotype = genotype[0:3]
            haplotype0_id = i*2
            haplotype1_id = i*2 + 1

            if genotype == "./.":
                self._missing_haplotypes.append(haplotype0_id)
                self._missing_haplotypes.append(haplotype1_id)
                n_missing_haplotypes += 2
            elif genotype == "0|1":
                self._haplotypes.append(haplotype1_id)
                n_haplotypes += 1
            elif genotype == "1|0":
                self._haplotypes.append(haplotype0_id)
                n_haplotypes += 1
            elif genotype == "1|1":
                self._haplotypes.append(haplotype0_id)
                self._haplotypes.append(haplotype1_id)
                n_haplotypes += 2
            elif genotype == "0|0":
                continue
            else:
                # Unphased. This is okay if homozygous
                if genotype == "1/1":
                    self._haplotypes.append(haplotype0_id)
                    self._haplotypes.append(haplotype1_id)
                elif genotype == "0/0":
                    continue
                else:
                    logging.error("Line %s" % line)
                    raise Exception("Variant not phased: %s" % genotype)

        #self._haplotypes = np.concatenate((self._haplotypes, haplotypes_tmp))

        self._from_nodes_to_n_haplotypes[from_node-self.graph_min_node] = n_haplotypes
        self._from_nodes_to_n_missing_haplotypes[from_node-self.graph_min_node] = n_missing_haplotypes

    def create(self):

        for i, line in enumerate(gzip.open(self.vcf_file_name)):
            if i % 1000 == 0:
                logging.info("%d lines processed. N haplotypes stored: %d" % (i, len(self._haplotypes)))

            line = line.decode("utf-8")
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    self._n_haplotypes = (len(line.split()) - 9) * 2
                    logging.info("There are %d haplotypes in this file" % self._n_haplotypes)
                continue

            variant_type = get_variant_type(line)
            l = line.split()
            ref_allele = l[3]
            variant_allele = l[4].lower()
            ref_offset = int(l[1]) - 1
            assert "," not in variant_allele, "Only biallelic variants are allowed. Line is not bialleleic"

            if variant_type == "SNP":
                edge = self._process_substitution(ref_offset, variant_allele)
            elif variant_type == "DELETION":
                edge = self._process_deletion(ref_offset+1, len(variant_allele)-1)
            else:
                continue

            self._store_processed_variant(l, edge)

        self._make_vcfmap()

    def _process_substitution(self, ref_offset, variant_bases):
        #logging.info("Processing SNP at pos %d" % ref_offset)
        node = self.reference_path.get_node_at_offset(ref_offset)
        node_offset = self.reference_path.get_node_offset_at_offset(ref_offset)

        assert node_offset == 0

        prev_node = self.reference_path.get_node_at_offset(ref_offset - 1)

        #print("Node: %d, offset: %d, node size: %d. Prev node: %d. Read base: %s" %
        #(node, node_offset, self.graph.blocks[node].length(), prev_node, variant_bases))

        # Try to find next node that matches read base
        #print("Next nodes: %s" % self.graph.adj_list[prev_node])
        for potential_next in self.graph.adj_list[prev_node]:
            if potential_next == node:
                continue
            node_seq = self.sequence_graph.get_sequence(potential_next, 0, 1)
            #print("  Next node %d has seq %s" % (potential_next, node_seq))
            if node_seq.lower() == variant_bases.lower():
                # Found a match!
                return (prev_node, potential_next)

        logging.error("Could not parse substitution at offset %d with bases %s" % (ref_offset, variant_bases))
        raise Exception("Parseerrror")

    def _process_deletion(self, ref_offset, deletion_length):
        logging.info("Processing deletion at ref pos %d with size %d" % (ref_offset, deletion_length))
        node = self.reference_path.get_node_at_offset(ref_offset)
        node_offset = self.reference_path.get_node_offset_at_offset(ref_offset)

        print("Processing deltion %s, %d, node offset %d" % (ref_offset, deletion_length, node_offset))
        assert node_offset == 0

        prev_node = self.reference_path.get_node_at_offset(ref_offset - 1)

        # Find next reference node with offset corresponding to the number of deleted base pairs
        next_ref_pos = ref_offset + deletion_length
        next_ref_node = self.reference_path.get_node_at_offset(ref_offset + deletion_length)
        if self.reference_path.get_node_offset_at_offset(next_ref_pos) != 0:
            logging.error("Offset %d is not at beginning of node" % next_ref_pos)
            logging.error("Node at %d: %d" % (next_ref_pos, next_ref_node))
            logging.error("Ref length in deletion: %s" % deletion_length)
            logging.info("Ref pos beginning of deletion: %d" % ref_offset)
            raise Exception("Deletion not in graph")

        self.n_deletions += 1
        return (prev_node, next_ref_node)


    def _process_insertion(self, ref_offset, read_offset):
        base = self.bam_entry.query_sequence[read_offset]
        node = self.reference_path.get_node_at_offset(ref_offset)
        node_offset = self.reference_path.get_node_offset_at_offset(ref_offset)
        node_size = self.graph.blocks[node].length()
        if node_offset != node_size - 1:
            # We are not at end of node, insertion is not represented in the graph, ignore
            return False

        # Find out which next node matches the insertion
        for potential_next in self.graph.adj_list[node]:
            if potential_next in self.linear_reference_nodes:
                continue  # Next should not be in linear ref

            #print("Processing insertion at ref offset %d with base %s" % (ref_offset, base))
            #print("  Node %d with offset %d. Node size: %d" % (node, node_offset, self.graph.blocks[node].length()))

            next_base = self.sequence_graph.get_sequence(potential_next, 0, 1).upper()
            if next_base == base.upper():
                self._variant_edges_detected.add((node, potential_next))
                self.n_insertions += 1
                #print("  Found next node %d with seq %s" % (potential_next, next_base))
                return


