#!usr/bin/python

import dendropy as sqd
import unittest


def partistic_dist_matrix(tree):
	"""
	Function to generate pair-wise phylogenetic
	distance matrix from a dendropy tree
	Input:
	tree - path to newick tree
	"""
	tns = sqd.TaxonNamespace()
	tre = open(tree, "r")
	tree = sqd.Tree.get(
		file = tre,
                schema = "newick",
                taxon_namespace=tns,
                preserve_underscore=True
                )
	taxon_list = []
	for idx, taxon1 in enumerate(tree.taxon_namespcae):
		str_taxon1 = str(taxon1.label)
	tre.close()


class TestRefBasisAnalysis(unittest.TestCase):
	def test_patristic_distance_matrix(self):
		tns = sqd.TaxonNamespace()
                tre = open("../test_data/simtree.tre", "r")
                tree = sqd.Tree.get(
                        file = tre,
                        schema = "newick",
                        taxon_namespace=tns,
                        preserve_underscore=True
                        )
		
		self.assertEqual((len(tree_taxon.namespace),len(tree.taxon_namespace))



if __name__ == '__main__':
	unittest.main()
