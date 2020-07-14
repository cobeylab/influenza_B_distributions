#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 13:56:28 2020

@author: marcosvieira

Analyzes trees of flu B HA and NA in the 2000s in New Zealand and Australia. 
Determines proximity to B/Victoria and B/Yamagata to detect reassortants.
"""

import sys
from dendropy import Tree
from Bio import AlignIO
import os
import csv
import re
import numpy


def determine_lineage_status(tree, Vic_ref, Yam_ref)



def main(argv):
    HA_tree_path = '../results/NZ-Aus_2000s_seqs_phylogenetics/HA_best_tree.tree'
    NA_tree_path = '../results/NZ-Aus_2000s_seqs_phylogenetics/NA_best_tree.tree'
    
    #parent_dir = os.path.dirname(HA_tree_path)

    HA_tree = Tree.get_from_path(HA_tree_path, schema='newick')
    NA_tree = Tree.get_from_path(NA_tree_path, schema='newick')
    
    alignment = AlignIO.read(node_seqs_path, 'fasta')
    # alignment = remove_mostly_gaps(alignment)
    alignment = list(alignment)

    annotation = annotate_aa_distances(tree, alignment, partition_points_file)
    annotated_tree = annotation['tree']

    mutation_counts = annotation['mutation_counts']

    annotated_tree.write(path=parent_dir + '/annotated_best_tree.nex', schema='nexus')

    with open(parent_dir + '/mutation_counts.csv', 'w') as mutation_counts_file:
        mutation_counts_file.write('mutation,region,n_internal,n_terminal,mean_n_desc,max_n_desc,min_n_desc\n')
        for mutation in mutation_counts.keys():
            region = mutation_counts[mutation]['region']
            n_internal = str(mutation_counts[mutation]['internal'])
            n_terminal = str(mutation_counts[mutation]['terminal'])
            mean_n_desc = str(numpy.mean(mutation_counts[mutation]['n_descendants']))
            max_n_desc = str(max(mutation_counts[mutation]['n_descendants']))
            min_n_desc = str(min(mutation_counts[mutation]['n_descendants']))
            mutation_counts_file.write(','.join([mutation, region, n_internal, n_terminal,
                                                 mean_n_desc, max_n_desc, min_n_desc]))
            mutation_counts_file.write('\n')


if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)
