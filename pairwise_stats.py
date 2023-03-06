#!/usr/bin/env python
import os
import sys
import pathlib
import dendropy

#print(sys.path)
path_root = pathlib.Path(__file__).parents[0]
sys.path.append(str(path_root) + '/moon_quiver')
from moon_quiver.moon_dragon import *

path_list = open("/home/jduque2/moon_logic/bs_path", 'r')
header = "base_tree\t"+"genTree\t"+"falsePosBip\t"+"falseNegBip\t"+"uRF\t"+"wRF\t"+"sackinBase\t"+"sackinGen"+"\n"
print(header)
for line in path_list:
	outline = ""
	t1 = "/home/jduque2/moon_logic/example_out/simtree.tre"
	t2=line.rstrip("\n").split("/")[4]
	tns,t1,t2 = pair_importer(t1,line.rstrip("\n"))
	tns,t1,t2 = standardize_namespace(tns, t1, t2)
	tns,t1,t2 = reroot_pair(tns, t1, t2, "SRR498373,SRR500494,SRR498369,SRR500493")
	print('adjusted:')
	print(t1.as_ascii_plot())
	print(t2.as_ascii_plot())
	outline=outline = "simtree" + "\t" + line.rstrip("_bestBiCon.tre\n").split("/")[5]
	outline = outline+"\t"+str(fp_and_n(tns,t1,t2)[0])+"\t"+str(fp_and_n(tns,t1,t2)[1])
	outline = outline+"\t"+str(pairwise_unweightedRF(t1,t2,tns))+"\t"+str(pairwise_weightedRF(t1,t2,tns)) +"\t"+str(sackin_ind(t1,tns))+"\t"+str(sackin_ind(t2,tns))+"\n"
	print(outline)

path_list.close()
