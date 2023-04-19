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

def reroot_by_multiTaxa(tns, t1, root_taxa):
        """Function to make sure two trees have same outgroup.
        Args:
        tns = dendropy taxon_namespace for both trees
        t1 = dendropy tree object
        root_taxa = list of taxon names
        """
        # rerooting at input mrca edge object need to standardize namespace prior to this
        mrca1 = t1.mrca(taxon_labels=root_taxa)
        assert mrca1
        t1.reroot_at_node(mrca1, update_bipartitions=False)
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


