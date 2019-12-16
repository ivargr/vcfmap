import logging
logging.basicConfig(level=logging.INFO)
import argparse
import sys
from offsetbasedgraph import Graph, SequenceGraph, NumpyIndexedInterval
from .mapcreator import MapCreator
from .vcfmap import VcfMap

def main():
    run_argument_parser(sys.argv[1:])


def make_vcfmap(args):
    graph = Graph.from_file(args.graph)
    sequence_graph = SequenceGraph.from_file(args.graph + ".sequences")
    linear_path = NumpyIndexedInterval.from_file(args.linear_path)
    linear_ref_nodes = linear_path.nodes_in_interval()
    logging.info("Done reading graphs")


    creator = MapCreator(graph, sequence_graph, linear_path, linear_ref_nodes, args.vcf_file)
    creator.create()
    creator.vcfmap.to_file(args.out_file_name)


def run_argument_parser(args):
    parser = argparse.ArgumentParser(
        description='VcfMap',
        prog='vcfmap',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=50, width=100))

    subparsers = parser.add_subparsers()

    # Create map
    sub = subparsers.add_parser("make")
    sub.add_argument("-g", "--graph", required=True)
    sub.add_argument("-v", "--vcf_file", required=True)
    sub.add_argument("-l", "--linear_path", required=True)
    sub.add_argument("-o", "--out_file_name", required=True)
    sub.set_defaults(func=make_vcfmap)

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)
    args.func(args)
