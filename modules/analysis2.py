#!/usr/bin/env python3

import os
import argparse
import pathlib
import pandas as pd
import subprocess
import datetime
import dateutil
import re
import dendropy
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='Generate a comparison DataFrame for multiple sequence alignments.')

    parser.add_argument('-o', '--ori_aln', type=str, required=True, help='Path to the original alignment in FASTA format.')
    parser.add_argument('-d', '--dir_new_alns', type=str, required=True, help='Path to the directory containing the new alignments in FASTA format.')
    parser.add_argument('-s', '--suffix', type=str, required=True, help='Suffix of the new alignment files.')
    parser.add_argument('-p', '--phylogeny', type=str, required=True, help="Path to the Newick-format tree file.")
    parser.add_argument('-r', '--result', type=str, required=True, help='Path to save the resulting comparison DataFrame as a CSV file.')

    args = parser.parse_args()

    br_df = build_branch_length_matrix(args.phylogeny)
    comp_df = generate_comparison_dataframe(args.ori_aln, args.dir_new_alns, args.suffix, br_df)

    comp_df.to_csv(args.result, index=False)

def build_branch_length_matrix(phylogeny):
    """
    Build a branch length matrix from a newick-format tree file.
    Returns a pandas DataFrame containing the branch lengths.
    """
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(path=phylogeny, schema="newick", taxon_namespace=taxa)
    pdm = tree.phylogenetic_distance_matrix()

    taxon_list = [str(taxon) for taxon in tree.taxon_namespace]

    br_df = pd.DataFrame(index=taxon_list, columns=taxon_list)

    for idx1, taxon1 in enumerate(tree.taxon_namespace):
        str_taxon1 = str(taxon1)
        for idx2 in range(idx1, len(tree.taxon_namespace)):
            taxon2 = tree.taxon_namespace[idx2]
            str_taxon2 = str(taxon2)
            weighted_patristic_distance = pdm.patristic_distance(taxon1, taxon2)
            br_df.at[str_taxon1, str_taxon2] = weighted_patristic_distance
            br_df.at[str_taxon2, str_taxon1] = weighted_patristic_distance

    return br_df


def make_comparison(new_seq, orig_seq, ref_seq):
    """
    Make comparison between the sequences and return counts of the results
    """

    output = [0] * 6

    standard_nucleotides = {'A', 'C', 'G', 'T'}
    gaps_or_degens = {'-', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N'}

    for num, nuc in enumerate(new_seq.strip()):
        nuc = nuc.upper()
        orig_nuc = orig_seq[num].upper()
        ref_nuc = ref_seq[num].upper()

        if nuc in gaps_or_degens:
            output[5] += 1

        elif orig_nuc == ref_nuc and nuc != orig_nuc:
            output[3] += 1

        elif nuc == orig_nuc:
            output[1] += 1

            if ref_nuc == orig_nuc:
                output[2] += 1

        elif nuc == ref_nuc:
            output[0] += 1

    output[4] = len(new_seq) - output[0] - output[1] - output[3] - output[5]

    return output


def generate_comparison_dataframe(ori_aln, dir_new_alns, suffix, br_df):
    print(br_df)
    # Read the original alignment
    ori_records = SeqIO.to_dict(SeqIO.parse(ori_aln, 'fasta'))

    # Get the list of alignment files in the directory
    alignment_files = [f for f in os.listdir(dir_new_alns) if f.endswith(suffix)]

    # Initialize an empty DataFrame for the output
    comp_df = pd.DataFrame()

    # Iterate through the alignment files
    for file in alignment_files:
        # Get the reference taxon from the file name
        reference_taxon = file[:-len(suffix)]

        # Read the new alignment
        new_records = SeqIO.to_dict(SeqIO.parse(os.path.join(dir_new_alns, file), 'fasta'))

        # Iterate through the taxa in the new alignment
        for taxon, new_seq_record in new_records.items():
            new_seq = str(new_seq_record.seq)
            orig_seq = str(ori_records[taxon].seq)
            ref_seq = str(ori_records[reference_taxon].seq)

            # Calculate the comparison using the make_comparison function
            comparison = make_comparison(new_seq, orig_seq, ref_seq)

            # Get the pairwise phylogenetic distance between new_seq and ref_seq
            pairwise_distance = br_df.at["'" + taxon + "'", "'" + reference_taxon + "'"]

            # Append the comparison result and pairwise distance to the DataFrame
            row_data = {
                'Taxon': taxon,
                'Reference': reference_taxon,
                'New_Matches_Ref_Not_Ori': comparison[0],
		'New_Matches_Ori': comparison[1],
		'TriSeqMatch': comparison[2],
		'TriSeqMismatch': comparison[3],
		'New_Mismatch_ORI_and_REF': comparison[4],
		'Degen': comparison[5],
                'Pairwise_Distance': pairwise_distance
            }
            comp_df = comp_df.append(row_data, ignore_index=True)

    # Set the DataFrame column order
    comp_df = comp_df[['Taxon','Reference', 'New_Matches_Ref_Not_Ori', 'New_Matches_Ori', 'TriSeqMatch', 'TriSeqMismatch', 'New_Mismatch_ORI_and_REF', 'Degen', 'Pairwise_Distance']]

    return comp_df


if __name__ == "__main__":
    main()
