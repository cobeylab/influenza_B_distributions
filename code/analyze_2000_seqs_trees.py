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
import os


def determine_lineage_status(tree, Vic_ref, Yam_ref):
    Vic_ref_tip = [node for node in tree.leaf_nodes() if Vic_ref in node.taxon.label]
    Yam_ref_tip = [node for node in tree.leaf_nodes() if Yam_ref in node.taxon.label]
    
    assert len(Vic_ref_tip) == 1
    assert len(Yam_ref_tip) == 1
    Vic_ref_tip, Yam_ref_tip = Vic_ref_tip[0], Yam_ref_tip[0]
    
    Vic_parent_node = Vic_ref_tip.parent_node
    Yam_parent_node = Yam_ref_tip.parent_node
    
    Vic_nodes = Vic_parent_node.preorder_iter()
    Yam_nodes = Yam_parent_node.preorder_iter()


    lineage_status = {}
    
    for node in list(Vic_nodes):
        if node.is_leaf():
            
            lineage_status[node.taxon.label.split('|')[0]] = 'B/Victoria'
    
    for node in list(Yam_nodes):
        if node.is_leaf():
            lineage_status[node.taxon.label.split('|')[0]] = 'B/Yamagata'
    
    return lineage_status
            
 
def annotate_tree(tree, lineage_status, segment_name):
    for node in tree.leaf_nodes():
        node.taxon.label
        
        if node.taxon.label.split('|')[0] in lineage_status.keys():
            
            node_lineage = lineage_status[node.taxon.label.split('|')[0]]            
            node.annotations.add_new(segment_name + '_lineage', node_lineage)
    return tree

        
def main(argv):
    HA_tree_path = '../results/NZ-Aus_2000s_seqs_phylogenetics/HA_best_tree.tree'
    NA_tree_path = '../results/NZ-Aus_2000s_seqs_phylogenetics/NA_best_tree.tree'
    
    parent_dir = os.path.dirname(HA_tree_path)

    HA_tree = Tree.get_from_path(HA_tree_path, schema='newick')
    NA_tree = Tree.get_from_path(NA_tree_path, schema='newick')
    
    Vic_ref = 'B/Victoria/2/87'
    Yam_ref = 'B/Yamagata/16/88'
    
    # HA lineage status (which lineage each strain is closest to in HA tree)
    HA_lineage_status = determine_lineage_status(HA_tree, Vic_ref, Yam_ref)
    
    # NA lineage status (which lineage each strain is closest to in NA tree)
    NA_lineage_status = determine_lineage_status(NA_tree, Vic_ref, Yam_ref)
    
    # Annotate both trees with both HA and NA lineage status
    HA_annotated_tree = annotate_tree(HA_tree, HA_lineage_status, 'HA')
    HA_annotated_tree = annotate_tree(HA_annotated_tree, NA_lineage_status, 'NA')
    
    
    NA_annotated_tree = annotate_tree(NA_tree, HA_lineage_status, 'HA')
    NA_annotated_tree = annotate_tree(NA_annotated_tree, NA_lineage_status, 'NA')
    
    HA_annotated_tree.write(path = parent_dir + '/annotated_HA_tree.nex', schema = 'nexus')
    NA_annotated_tree.write(path = parent_dir + '/annotated_NA_tree.nex', schema = 'nexus')
 

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)
