#! usr/bin/python3.11
import dendropy
import unittest

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


class TestRootPair(unittest.TestCase):
	def test_mrcaLeafEquivalence(self):
		tns = dendropy.TaxonNamespace()
		tre = "../bestTrees/SRR498444_bestBiCon.tre"
		tree = dendropy.Tree.get(
			file=open(tre, 'r'),
			schema="newick",
			taxon_namespace=tns,
			preserve_underscores=True
			)
		tree.prune_taxa_with_labels(["ref_SRR498444"])
		root_taxa = ["SRR498373","SRR500494","SRR498369","SRR500493"]## Multiple outgroup taxa
		t_mrca = tree.mrca(taxon_labels=root_taxa)
		print(t_mrca)
		assert t_mrca
		ori_leaves = [leaf.taxon.label for leaf in t_mrca.leaf_iter()]
		print(root_taxa)
		print(ori_leaves)
		self.assertCountEqual(ori_leaves, root_taxa)



"""






example_tree2 = "((D:1.0,F:1.0,((A:1.0,B:1.0):1.0,C:1.0):2.0):0.5,E:0.5);"
tree2 = dendropy.Tree.get(data = example_tree2, schema = "newick")



root_taxa =  ["C","D","E","F"]## Multiple outgroup taxa

for tree in [tree1, tree2]:
    leaves_before = [leaf for leaf in tree.leaf_node_iter()]
    reroot_pair(tns=tree.taxon_namespace, t1=tree1, t2=tree2, root_taxa=root_taxa)

    leaves_after = [leaf for leaf in tree.leaf_node_iter()]

    if not leaves_before == leaves_after:
        print("WHYYYYYYYYY Where did my LEAF GO")
        print("Len leaves_before is {}".format(len(leaves_before)))
        print("Len leaves_after is {}".format(len(leaves_after)))


    root = tree.mrca(leaves_after)
    mrca_og_leaf = tree.mrca(taxon_labels=["A"]+ root_taxa)
    assert(root == mrca_og_leaf)

    internal_node = tree.mrca(taxon_labels=["A","B"])
    assert(root != internal_node)


"""
