#! /usr/bin/env python3

######################################################################
# Goal: Automate Dendropy based exploratory operations for personal use
# Author: A3M5
# Date: 10/07/2022
######################################################################

import dendropy
from dendropy.calculate import treecompare as tc
from dendropy.calculate import treemeasure as tm

# custom function for carrying out the RF weighted distances between two trees
def pair_importer(tre1, tre2):
    t1 = dendropy.Tree.get(
        file=open(tre1, 'r'),
        schema="newick"
        )
    t2 = dendropy.Tree.get(
        file=open(tre2, 'r'),
        schema="newick",
        taxon_namespace=t1.taxon_namespace
        )
    t1.encode_bipartitions()
    t2.encode_bipartitions()
    return t1,t2

def pairwise_weightedRF(tre1, tre2):
    t1,t2 = pair_importer(tre1,tre2)
    t1.encode_bipartitions()
    t2.encode_bipartitions()
    return tc.weighted_robinson_foulds_distance(t1, t2)

def pairwise_unweightedRF(tre1, tre2):
    t1,t2 = pair_importer(tre1,tre2)
    t1.encode_bipartitions()
    t2.encode_bipartitions()
    return tc.unweighted_robinson_foulds_distance(t1, t2)

def fp_and_n(tre1, tre2):
    t1,t2 = pair_importer(tre1,tre2)
    t1.encode_bipartitions()
    t2.encode_bipartitions()
    return tc.false_positives_and_negatives(t1, t2)

def colles_tree_imb(tre1):
    tns = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get(
        file=open(tre1, 'r'),
        schema="newick"
        )
    t1.encode_bipartitions()
    return tm.colless_tree_imbalance(t1)

def sackin_ind(tre1):
    tns = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get(
        file=open(tre1, 'r'),
        schema="newick"
        )
    t1.encode_bipartitions()
    return tm.sackin_index(t1)




# Open file for output
outfile = open("dendOut_weightedRF.tsv", "w")
# Import list of query trees generated via ep
path_list = open("/home/jduque2/moon_logic/ep_bestTrees.path", 'r')
# Include path for the base tree used by TTR to generate TTR output
ttr_tree = "/home/jduque2/moon_logic/example_out/simtree.tre"

header = "base_tree\t"+"genTree\t"+"falsePosBip\t"+"falseNegBip\t"+"uRF\t"+"wRF\t"+"sackinBase\t"+"sackinGen"+"\n"
outfile.write(header)
# for loop for assign tree files to query files for wRF operations
for line in path_list:
     outline = ""
     print("UnweightedRF: ", pairwise_unweightedRF(ttr_tree,line.rstrip("\n")))
     outline = "simtree" + "\t" + line.split("/")[4].rstrip("_ep_out")
     outline = outline+"\t"+str(fp_and_n(ttr_tree,line.rstrip("\n"))[0])+"\t"+str(fp_and_n(ttr_tree,line.rstrip("\n"))[1])
     outline = outline+"\t"+str(pairwise_unweightedRF(ttr_tree,line.rstrip("\n")))+"\t"+str(pairwise_weightedRF(ttr_tree,line.rstrip("\n"))) +"\t"+str(sackin_ind(ttr_tree))+"\t"+str(sackin_ind(line.rstrip("\n")))+"\n"
     outfile.write(outline)
path_list.close()
