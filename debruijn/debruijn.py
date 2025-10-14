#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import statistics
import textwrap
from pathlib import Path
from random import randint
from typing import Iterator, Dict, List
import random
import matplotlib
import matplotlib.pyplot as plt
from networkx import (
    DiGraph,
    all_simple_paths,
    lowest_common_ancestor,
    has_path,
    random_layout,
    draw_networkx_nodes,
    draw_networkx_edges,
)

random.seed(9001)
matplotlib.use("Agg")

__author__ = "Assa DIABIRA"
__copyright__ = "Universite Paris Cité"
__credits__ = ["Assa DIABIRA"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Assa DIABIRA"
__email__ = "assa.diabira@etu.u-paris.fr"
__status__ = "Developpement"


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h")
    parser.add_argument(
        "-i", dest="fastq_file", type=isfile, required=True, help="Fastq file"
    )
    parser.add_argument(
        "-k", dest="kmer_size", type=int, default=22, help="k-mer size (default 22)"
    )
    parser.add_argument(
        "-o",
        dest="output_file",
        type=Path,
        default=Path(os.curdir + os.sep + "contigs.fasta"),
        help="Output contigs in fasta file (default contigs.fasta)",
    )
    parser.add_argument(
        "-f", dest="graphimg_file", type=Path, help="Save graph as an image (png)"
    )
    return parser.parse_args()


def read_fastq(fastq_file: Path) -> Iterator[str]:
    """Extract reads from fastq files.

    :param fastq_file: (Path) Path to the fastq file.
    :return: A generator object that iterate the read sequences.
    """
    with open(fastq_file, "r", encoding="utf-8") as file:
        for _, seq, _, _ in zip(*[file] * 4):
            yield seq.strip()


def cut_kmer(read: str, kmer_size: int) -> Iterator[str]:
    """Cut read into kmers of size kmer_size.

    :param read: (str) Sequence of a read.
    :return: A generator object that provides the kmers (str) of size kmer_size.
    """
    for i in range(len(read) - kmer_size + 1):
        yield read[i : i + kmer_size]


def build_kmer_dict(fastq_file: Path, kmer_size: int) -> Dict[str, int]:
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_dict = {}

    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1

    return kmer_dict


def build_graph(kmer_dict: Dict[str, int]) -> DiGraph:
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    graph = DiGraph()

    for kmer, weight in kmer_dict.items():
        # prefix = first k-1 characters
        # suffix = last k-1 characters
        prefix = kmer[:-1]
        suffix = kmer[1:]

        # Add edge with weight (occurrence count)
        graph.add_edge(prefix, suffix, weight=weight)

    return graph


def remove_paths(
    graph: DiGraph,
    path_list: List[List[str]],
    delete_entry_node: bool,
    delete_sink_node: bool,
) -> DiGraph:
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for path in path_list:
        if not path:
            continue

        # Calculate slice indices
        # If delete_entry_node is False, start from index 1 (skip first)
        # If delete_sink_node is False, end at len-1 (skip last)
        start_idx = 0 if delete_entry_node else 1
        end_idx = len(path) if delete_sink_node else len(path) - 1

        # Remove nodes in the calculated range
        for node in path[start_idx:end_idx]:
            if graph.has_node(node):
                graph.remove_node(node)

    return graph


def select_best_path(
    graph: DiGraph,
    path_list: List[List[str]],
    path_length: List[int],
    weight_avg_list: List[float],
    delete_entry_node: bool = False,
    delete_sink_node: bool = False,
) -> DiGraph:
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    if not path_list or len(path_list) == 1:
        return graph

    # Criterion 1: Compare by weight (standard deviation > 0 means differences exist)
    weight_std = statistics.stdev(weight_avg_list) if len(weight_avg_list) > 1 else 0

    if weight_std > 0:
        # Select path with highest weight
        max_weight = max(weight_avg_list)
        best_indices = [i for i, w in enumerate(weight_avg_list) if w == max_weight]
    else:
        best_indices = list(range(len(path_list)))

    # Criterion 2: If tie, compare by length
    if len(best_indices) > 1:
        length_std = (
            statistics.stdev([path_length[i] for i in best_indices])
            if len(best_indices) > 1
            else 0
        )

        if length_std > 0:
            max_length = max(path_length[i] for i in best_indices)
            best_indices = [i for i in best_indices if path_length[i] == max_length]

    # Criterion 3: If still tie, random selection
    if len(best_indices) > 1:
        best_idx = randint(0, len(best_indices) - 1)
        best_idx = best_indices[best_idx]
    else:
        best_idx = best_indices[0]

    # Remove all paths except the best one
    paths_to_remove = [path_list[i] for i in range(len(path_list)) if i != best_idx]

    return remove_paths(graph, paths_to_remove, delete_entry_node, delete_sink_node)


def path_average_weight(graph: DiGraph, path: List[str]) -> float:
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean(
        [d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)]
    )


def solve_bubble(graph: DiGraph, ancestor_node: str, descendant_node: str) -> DiGraph:
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    # Find all simple paths between ancestor and descendant
    all_paths = list(all_simple_paths(graph, ancestor_node, descendant_node))

    if len(all_paths) <= 1:
        return graph

    # Calculate path lengths (number of nodes)
    path_lengths = [len(path) for path in all_paths]

    # Calculate average weights for each path
    path_weights = [path_average_weight(graph, path) for path in all_paths]

    # Select best path
    # For bubbles: we don't delete the entry (ancestor) or sink (descendant) nodes
    return select_best_path(
        graph,
        all_paths,
        path_lengths,
        path_weights,
        delete_entry_node=False,
        delete_sink_node=False,
    )


def simplify_bubbles(graph: DiGraph) -> DiGraph:
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble = False
    ancestor = None
    node = None

    # Iterate through all nodes to find bubbles
    for node in list(graph.nodes()):
        if not graph.has_node(node):
            continue

        # Get predecessors of current node
        predecessors = list(graph.predecessors(node))

        # If node has multiple predecessors, check for common ancestor
        if len(predecessors) > 1:
            # Check all unique combinations of predecessors
            for i, pred_i in enumerate(predecessors):
                for pred_j in predecessors[i + 1 :]:
                    # Find lowest common ancestor
                    ancestor = lowest_common_ancestor(graph, pred_i, pred_j)

                    if ancestor is not None:
                        # Bubble detected between ancestor and current node
                        bubble = True
                        break

                if bubble:
                    break

        if bubble:
            break

    # Recursive approach: simplify one bubble then check again
    if bubble and ancestor is not None and node is not None:
        graph = simplify_bubbles(solve_bubble(graph, ancestor, node))

    return graph


def solve_entry_tips(graph: DiGraph, starting_nodes: List[str]) -> DiGraph:
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of starting nodes
    :return: (nx.DiGraph) A directed graph object
    """
    tip = False

    # Find nodes with multiple predecessors
    for node in list(graph.nodes()):
        if not graph.has_node(node):
            continue

        predecessors = list(graph.predecessors(node))

        # If node has multiple predecessors
        if len(predecessors) > 1:
            # Check if any predecessors are starting nodes or connected to starting nodes
            paths = []

            for start_node in starting_nodes:
                if not graph.has_node(start_node):
                    continue

                # Check if there's a path from this starting node to current node
                if has_path(graph, start_node, node):
                    # Get all simple paths
                    for path in all_simple_paths(graph, start_node, node):
                        paths.append(path)

            # If we found multiple paths from starting nodes to this node
            if len(paths) > 1:
                # Calculate metrics
                path_lengths = [len(path) for path in paths]
                path_weights = [path_average_weight(graph, path) for path in paths]

                # Remove entry tips: delete entry nodes, keep convergence node
                graph = select_best_path(
                    graph,
                    paths,
                    path_lengths,
                    path_weights,
                    delete_entry_node=True,
                    delete_sink_node=False,
                )
                tip = True
                break

    # Recursive approach
    if tip:
        graph = solve_entry_tips(graph, get_starting_nodes(graph))

    return graph


def solve_out_tips(graph: DiGraph, ending_nodes: List[str]) -> DiGraph:
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :param ending_nodes: (list) A list of ending nodes
    :return: (nx.DiGraph) A directed graph object
    """
    tip = False

    # Find nodes with multiple successors
    for node in list(graph.nodes()):
        if not graph.has_node(node):
            continue

        successors = list(graph.successors(node))

        # If node has multiple successors
        if len(successors) > 1:
            # Check if any successors are ending nodes or connected to ending nodes
            paths = []

            for end_node in ending_nodes:
                if not graph.has_node(end_node):
                    continue

                # Check if there's a path from current node to this ending node
                if has_path(graph, node, end_node):
                    # Get all simple paths
                    for path in all_simple_paths(graph, node, end_node):
                        paths.append(path)

            # If we found multiple paths from this node to ending nodes
            if len(paths) > 1:
                # Calculate metrics
                path_lengths = [len(path) for path in paths]
                path_weights = [path_average_weight(graph, path) for path in paths]

                # Remove out tips: keep divergence node, delete sink nodes
                graph = select_best_path(
                    graph,
                    paths,
                    path_lengths,
                    path_weights,
                    delete_entry_node=False,
                    delete_sink_node=True,
                )
                tip = True
                break

    # Recursive approach
    if tip:
        graph = solve_out_tips(graph, get_sink_nodes(graph))

    return graph


def get_starting_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    return [node for node in graph.nodes() if graph.in_degree(node) == 0]


def get_sink_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    return [node for node in graph.nodes() if graph.out_degree(node) == 0]


def get_contigs(
    graph: DiGraph, starting_nodes: List[str], ending_nodes: List[str]
) -> List:
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contigs = []

    for start_node in starting_nodes:
        for end_node in ending_nodes:
            # Check if path exists between start and end
            if has_path(graph, start_node, end_node):
                # Get all simple paths
                for path in all_simple_paths(graph, start_node, end_node):
                    # Build contig sequence
                    # First node gives us the first k-mer (which is k-1 length node)
                    contig = path[0]

                    # Each subsequent node adds one nucleotide (last character)
                    for node in path[1:]:
                        contig += node[-1]

                    # Add tuple (contig, length)
                    contigs.append((contig, len(contig)))

    return contigs


def save_contigs(contigs_list: List[str], output_file: Path) -> None:
    """Write all contigs in fasta format

    :param contigs_list: (list) List of [contiguous sequence and their length]
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w", encoding="utf-8") as file:
        for i, (contig, length) in enumerate(contigs_list):
            # Write header with contig number and length
            file.write(f">contig_{i} len={length}\n")
            # Write sequence wrapped at 80 characters per line
            wrapped_sequence = textwrap.fill(contig, width=80)
            file.write(wrapped_sequence + "\n")


def draw_graph(graph: DiGraph, graphimg_file: Path) -> None:  # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (Path) Path to the output file
    """
    plt.figure()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d["weight"] > 3]
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d["weight"] <= 3]

    # Draw the graph with networkx
    pos = random_layout(graph)
    draw_networkx_nodes(graph, pos, node_size=6)
    draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    draw_networkx_edges(
        graph, pos, edgelist=esmall, width=6, alpha=0.5, edge_color="b", style="dashed"
    )
    # save image
    plt.savefig(graphimg_file.resolve())
    plt.close()


# ==============================================================
# Main program
# ==============================================================
def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Step 1: Read file and build graph
    print("Step 1: Building k-mer dictionary and De Bruijn graph...")
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    print(f"  K-mers: {len(kmer_dict)}")
    print(f"  Graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")

    # Step 2: Resolve bubbles
    print("\nStep 2: Simplifying bubbles...")
    graph = simplify_bubbles(graph)
    print(f"  Graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")

    # Step 3: Resolve entry and out tips
    print("\nStep 3: Resolving tips...")
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)
    print(f"  Graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")

    # Update starting and ending nodes after simplification
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    print(f"  Entry nodes: {len(starting_nodes)}, Exit nodes: {len(ending_nodes)}")

    # Step 4: Extract and save contigs
    print("\nStep 4: Extracting contigs...")
    contigs = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs, args.output_file)
    print(f"  Found {len(contigs)} contig(s)")

    # Display statistics
    if contigs:
        total_length = sum(length for _, length in contigs)
        max_length = max(length for _, length in contigs)
        min_length = min(length for _, length in contigs)
        avg_length = total_length / len(contigs)

        print("\n  Contig statistics:")
        print(f"    Total length: {total_length} bp")
        print(f"    Average length: {avg_length:.1f} bp")
        print(f"    Min length: {min_length} bp")
        print(f"    Max length: {max_length} bp")

    print(f"\nContigs saved to: {args.output_file}")

    # Draw graph if requested
    if args.graphimg_file:
        draw_graph(graph, args.graphimg_file)
        print(f"Graph image saved to: {args.graphimg_file}")


if __name__ == "__main__":  # pragma: no cover
    main()
