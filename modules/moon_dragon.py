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
    tns = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get(
        file=open(tre1, 'r'),
        schema="newick",
        taxon_namespace=tns,
	rooting="force-rooted",
	preserve_underscores=True,
        )
    t2 = dendropy.Tree.get(
        file=open(tre2, 'r'),
        schema="newick",
        taxon_namespace=tns,
	rooting="force-rooted",
	preserve_underscores=True
	)
    t1.encode_bipartitions()
    t2.encode_bipartitions()
    return tns, t1, t2

def reroot_pair(tns, t1, t2, root_taxa):
	# rerooting at input mrca edge object need to standardize namespace prior to this
	mrca1 = t1.mrca(taxon_labels=root_taxa)
	mrca2 = t2.mrca(taxon_labels=root_taxa)
	t1.reroot_at_node(mrca1, update_bipartitions=True)
	t2.reroot_at_node(mrca2, update_bipartitions=True)
	return tns,t1,t2

def standardize_namespace(tns, t1, t2):
    """requires two trees with a common taxon namespace"""
    tree_1_taxa = set()
    tree_2_taxa = set()

    for tip in t1.leaf_node_iter():
        tree_1_taxa.add(tip.taxon)

    for tip in t2.leaf_node_iter():
        tree_2_taxa.add(tip.taxon)

    shared_taxa = tree_1_taxa.intersection(tree_2_taxa)

    assert(len(shared_taxa) >= 1)
    #print("These two trees have {s} shared taxa".format(s=len(shared_taxa)))

    t1.retain_taxa(shared_taxa)
    t2.retain_taxa(shared_taxa)
    return(tns, t1, t2)


def pairwise_weightedRF(t1, t2, tns):
    #tns,t1,t2 = pair_importer(tre1,tre2)
    #t1.encode_bipartitions()
    #t2.encode_bipartitions()
    return tc.weighted_robinson_foulds_distance(t1, t2)

def pairwise_unweightedRF(t1, t2, tns):
    #tns,t1,t2 = pair_importer(tre1,tre2)
    #t1.encode_bipartitions()
    #t2.encode_bipartitions()
    return tc.symmetric_difference(t1, t2)

def fp_and_n(tns, t1, t2):
    #tns,t1,t2 = pair_importer(tre1,tre2)
    #t1.encode_bipartitions()
    #t2.encode_bipartitions()
    return tc.false_positives_and_negatives(t1, t2)

def colles_tree_imb(t1):
    #tns = dendropy.TaxonNamespace()
    #t1 = dendropy.Tree.get(
    #    file=open(tre1, 'r'),
    #    schema="newick"
    #    )
    #t1.encode_bipartitions()
    return tm.colless_tree_imbalance(t1)

def sackin_ind(t1,tns):
    #tns = dendropy.TaxonNamespace()
    #t1 = dendropy.Tree.get(
    #    file=open(tre1, 'r'),
    #    schema="newick"
    #    )
    #t1.encode_bipartitions()
    return tm.sackin_index(t1)


def main():
	# Open file for output
	outfile = open("dendOut_weightedRF.tsv", "w")
	# Import list of query trees generated via ep
	path_list = open("/home/josue/Desktop/op_lunar_trees/bestPath.txt", 'r')
	# Include path for the base tree used by TTR to generate TTR output
	ttr_tree = "simtree.tre"

	header = "base_tree\t"+"genTree\t"+"falsePosBip\t"+"falseNegBip\t"+"uRF\t"+"wRF\t"+"sackinBase\t"+"sackinGen"+"\n"
	outfile.write(header)
	# for loop for assign tree files to query files for wRF operations
	for line in path_list:
	     outline = ""
	     t2=line.rstrip("\n")
	     print("UnweightedRF: ", pairwise_unweightedRF(ttr_tree,t2))
	     outline = "simtree" + "\t" + line.rstrip("_bestBiCon.tre\n")
	     outline = outline+"\t"+str(fp_and_n(ttr_tree,t2)[0])+"\t"+str(fp_and_n(ttr_tree,t2)[1])
	     outline = outline+"\t"+str(pairwise_unweightedRF(ttr_tree,t2))+"\t"+str(pairwise_weightedRF(ttr_tree,t2)) +"\t"+str(sackin_ind(ttr_tree))+"\t"+str(sackin_ind(t2))+"\n"
	     outfile.write(outline)
	     tns,t1,t2 = pair_importer(ttr_tree,t2)
	     print("\n___________________________________New Analysis______________________________________________")
	     print(t1.as_ascii_plot())
	     print(t2.as_ascii_plot())
	     print("\n_____________________________________________________________________________________")
	     print("_____________________________________________________________________________________")
	     print("_____________________________________________________________________________________")


	path_list.close()


if __name__ == '__main__':
    main()


### OUT GROUP
### SRR498373
### SRR500494
### SRR498369
### SRR500493


################## TODO LIST #############################
####    Contruct branch trimmer function              ####
####    Contruct namesspace overlap forcer function   ####
##########################################################
