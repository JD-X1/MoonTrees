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

def fix_names(tree_file, df, sample_name_column):
    """
    Find and replace names if the names in the tree are the SRA numbers
    """
    tree = open(tree_file, 'r').read()

    for idx, row in df.iterrows():
        # print(row['SampleName'])
        # print(row['Run'])
        # print("###")

        sra_regex = row['Run']
        sra_compile = re.compile(sra_regex)
        find_sra = re.findall(sra_compile, tree)
        if find_sra:
            tree = tree.replace(row['Run'], row['SampleName'])

        else:
            print("Problems on regex in: ", tree_file)

    tree_file_name = tree_file.rsplit('/', 1)
    updated_tree_file_name = 'names_updated_' + tree_file_name[1]
    output_tree = open(updated_tree_file_name, 'w')
    output_tree.write(tree)
    output_tree.close()

def pair_importer(tre1, tre2):
    tns = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get(
        file=open(tre1, 'r'),
        schema="newick",
	    taxon_namespace=tns, preserve_underscores=True
	    #root="force-unrooted"
        )
    t2 = dendropy.Tree.get(
        file=open(tre2, 'r'),
        schema="newick",
        taxon_namespace=tns, preserve_underscores=True
        #root="force-unrooted"
	)
    # two lines below are serving very specific functions that is specific
    # to my runs as of Oct 11, 2022 comment these lines out to include all
    # taxa in your tree. It shouldn't really cause problems, but better
    # safe than sorry.
    #prune_target = str(tre2.split("/")[5])
    #prune_target = prune_target.rstrip("_best.tre").lstrip("sim_")
    #print(t2.as_ascii_plot())
    t2.prune_taxa_with_labels([prune_target])
    t1.encode_bipartitions()
    t2.encode_bipartitions()
    return tns,t1,t2

def pairwise_weightedRF(tre1, tre2):
    tns,t1,t2 = pair_importer(tre1,tre2)
    t1.encode_bipartitions()
    t2.encode_bipartitions()
    return tc.weighted_robinson_foulds_distance(t1, t2)

def pairwise_unweightedRF(tre1, tre2):
    tns,t1,t2 = pair_importer(tre1,tre2)
    t1.encode_bipartitions()
    t2.encode_bipartitions()
    return tc.unweighted_robinson_foulds_distance(t1, t2)

def fp_and_n(tre1, tre2):
    tns,t1,t2 = pair_importer(tre1,tre2)
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
     outline = "simtree" + "\t" + line.split("/")[5].rstrip("_best.tre\n")
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
