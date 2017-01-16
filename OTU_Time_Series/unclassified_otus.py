#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 13:56:03 2017

To confirm that the effect is still present, we used FastTree 2.1.9 (compiled 
with double precision) to create a new phylogeny of most of the sequences
contained in the original tree.  The PyNAST aligner was used on the fasta 
file listed below to create a new alignment file. 


./fastree/FastTree -nt -gtr -gamma < pynast_aligned/unique.dbOTU.nonchimera_aligned.fasta 
> unique.dbOTU.nonchimera.tree

It is looking like this will take a long time. So we will probably want to 
restrict tree construction to only those sequences that were successfully aligned
by PyNAST. It would be useful at this point to check which sequences were removed
and from which group. Among those that remain, we want to retain all 
unclassifiables, all taxa classified as "bacteria" that are 1 standard deviation 
above the median node depth / distance of the self-same group. Once we have a 
final list of these high branching nodes, we can select an equal number of better
classified taxa. These will be selected at random. 

The alignment file is on MARCC in ~/work/kaw/unclassifiable/pynast_aligned
and is called unique.dbOTU.nonchimera_aligned.fasta. 

@author: login
"""

import pandas as pd
import numpy as np
import sys, os
from otu_ts_support import inferTaxaLevel
import dendropy

# Relevant File Names OTU Table
seq_file='unique.dbOTU.nonchimera.fasta'
otu_matrix_file = 'unique.dbOTU.nonchimera.mat.rdp'

from Bio import SeqIO
records = list(SeqIO.parse(seq_file, "fasta"))

print "\nReading in OTU Table"
otu_table = pd.read_csv(otu_matrix_file, sep="\t", index_col=0)
taxa_series = otu_table.ix[:, -1]

print "\nReading in OTU tree & calculating node distances"

tree_files = ['unique.dbOTU.nonchimera.edit.tree', 'unique.dbOTU.tree']
for tree_file in tree_files:
    print "\nReading in {}".format(tree_file)    
    masterTree = dendropy.Tree.get(path=tree_file, schema='newick', rooting="force-rooted")
    seq_List = masterTree.taxon_namespace
    node_Dists = masterTree.calc_node_root_distances()
    seq_dist_dict = {str(i)[1:-1]:j for i, j in zip(seq_List, node_Dists)}

    print "Calculating the average distance / depth per taxanomic level"
    
    taxa_frame = inferTaxaLevel(taxa_series)
    node_distance = np.full(len(taxa_series), np.nan)
    taxa_frame["Node Distance"] = node_distance
    for i in taxa_frame.index:
        if i in seq_dist_dict.keys():
            taxa_frame.ix[i, "Node Distance"] = seq_dist_dict[i]
    
    seq_levels = {}
    for i in masterTree.leaf_nodes():
        lvl = i.level()
        seq = i.taxon.label
        seq_levels[seq] = lvl
    
    node_lvl = np.full(len(taxa_series), np.nan)
    taxa_frame["Node Level"] = node_lvl
    for i in taxa_frame.index:
        if i in seq_levels.keys():
            taxa_frame.ix[i, "Node Level"] = seq_levels[i]
    
    frac_of_data_assigned_terminal_taxa = []
    for t_d in range(0,7):
        sub_df_bool = taxa_frame["Taxa Depth"] == float(t_d)
        sub_sr = taxa_frame[sub_df_bool].ix[:, "Node Distance"].dropna()
        print "Classification Depth", t_d
        print "Node Distance (median)", np.median(sub_sr)
        print "Node Distance (variance)", np.var(sub_sr)
        assigned_frac = len(sub_sr.index)/float(len(taxa_frame.index))
        print "N", len(sub_sr.index)
        frac_of_data_assigned_terminal_taxa.append(assigned_frac)
    print ""    
    for t_d in range(0,7):
        sub_df_bool = taxa_frame["Taxa Depth"] == float(t_d)
        sub_sr = taxa_frame[sub_df_bool].ix[:, "Node Level"].dropna()
        print "Classification Depth", t_d
        print "Node Level (median)", np.median(sub_sr)
        print "Node Level (variance)", np.var(sub_sr)
        assigned_frac = len(sub_sr.index)/float(len(taxa_frame.index))
        print "N", len(sub_sr.index)

    #TODO: finish this plot
    #plt.bar(range(1,8), frac_of_data_assigned_terminal_taxa)