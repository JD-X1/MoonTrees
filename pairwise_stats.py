#!/usr/bin/env python
import os
import sys
import pathlib
import dendropy
import argparse
import pandas as pd
from dendropy.calculate import treecompare as tc
from dendropy.calculate import treemeasure as tm

#print(sys.path)
path_root = pathlib.Path(__file__).parents[0]
sys.path.append(str(path_root) + '/modules')
from modules.moon_dragon import *

def main():
	parser = argparse.ArgumentParser(description="Generate a dataframe for multiple tree files that share the saem taxa")
	parser.add_argument('-d', '--tree_directory', type=str, required=True, help="Path to tree files generated from EP")
	parser.add_argument('-o', '--out', type=str, required=True, help="Name for output csv file")
	parser.add_argument('-s', '--source_tree', type=str, required=True, help="Path to tree used as input for TTR")
	args = parser.parse_args()
	df = pairwise_tree_stats(args.tree_directory, args.source_tree)
	df.to_csv(args.out, index = False)

def pairwise_tree_stats(tree_dir,  base_tree):
	tree_files = [f for f in os.listdir(tree_dir) if f.endswith(".tre")]
	df = pd.DataFrame()
	btree = base_tree
	tns = dendropy.TaxonNamespace()
	t1 = dendropy.Tree.get(
		file = open(base_tree),
		schema = "newick",
		taxon_namespace = tns
		)


	for file  in tree_files:
		t2_name = str(file)[:-len("_biCons.tre")]
		t2 = dendropy.Tree.get(
                	file = open(base_tree),
                	schema = "newick",
                	taxon_namespace = tns
                	)
		tns,t1,t2 = standardize_namespace(tns, t1, t2)
		t1.encode_bipartitions()
		t2.encode_bipartitions()
		fal_pos_neg = tc.false_positives_and_negatives(t1, t2)
		
		row_data = {
			'BaseTree':btree,
			'genTree':t2_name,
			'falePosBip':fal_pos_neg[0],
			'falseNegBip':fal_pos_neg[1],
			'uRF':tc.symmetric_difference(t1,t2),
			'wRF':tc.weighted_robinson_foulds_distance(t1,t2),
			'sackinBase':tm.sackin_index(t1),
			'sackinGen':tm.sackin_index(t2)
		}
		df = df.append(row_data, ignore_index = True)
	df =df[['BaseTree', 'genTree', 'falePosBip', 'falseNegBip', 'uRF', 'wRF', 'sackinBase', 'sackinGen']]

	return df

if __name__ == "__main__":
    main()
