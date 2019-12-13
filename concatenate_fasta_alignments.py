#!/usr/bin/env python3
"""
Script for concatenating protein FASTA alignments of different genes in order to make a species phylogeny tree.
Argument 1: path to folder with FASTA alignments that have been processed with Gblocks. Each file should contain a multiple sequence alignment of one 
    protein across multiple species. The FASTA headers should start with a species identifier that is separated from the rest of the header with underscore, e.g. PocGH01_10025200.1-p1.
    All input FASTA files have to have the same number of sequences. The order of the sequences does not have to be the same in all input files.
Output: concatenated protein alignments (FASTA), one per each species.
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License
"""

import hep_project_shared_functions as hpf
import os
from collections import defaultdict
import sys
import argparse


def get_in_files(in_folder):
    """
    Gets a list of files in a folder without including hidden files
    """
    in_files = os.listdir(in_folder)
    in_files = [n for n in in_files if n.startswith(".") == False]
    return in_files


def load_fasta_files_as_dict(in_folder):
    """
    Input: path to folder with Gblocks output files (alignments)
    Output: a defaultdict where the keys are file names and values are lists.
        The first item in each list is a FASTA header and the rest of the items are the corresponding protein sequence lines from the 
        FASTA file
    """
    
    in_files = get_in_files(in_folder)
    files_dict = defaultdict(list)

    for in_file in in_files:
        fasta_data = list(hpf.read_fasta_in_chunks(in_folder + "/" + in_file))
        cropped_fasta_headers = list()
        for pair in fasta_data:
            header = pair[0].split("_")[0]
            cropped_fasta_headers.append(header)
        cropped_fasta_headers.sort()

        sorted_fasta_data = list()
        for cropped_header in cropped_fasta_headers:
            for fasta_item in fasta_data:
                if fasta_item[0].startswith(cropped_header):
                    sorted_fasta_data.append(fasta_item)
        files_dict[in_file] = sorted_fasta_data
    files_dict2 = defaultdict(list)
    for item in files_dict:
        if len(files_dict[item]) > 0:
            files_dict2[item] = files_dict[item]
    return files_dict2


def genes_count_same_in_all_files(files_dict, in_files):
    """
    Function for checking if all alignments in input have an equal number of aligned aligned proteins
    (returns True if the numbers are equal, False if not)
    """
    all_identical_flag = False
    gene_counts = list()
    for in_file in in_files:
        if in_file in files_dict:
            genes_count = len(files_dict[in_file])
            gene_counts.append(genes_count)
    if gene_counts.count(gene_counts[0]) == len(gene_counts):
        all_identical_flag = True
    return all_identical_flag


def main(in_folder):

    in_files = get_in_files(in_folder)
    files_dict = load_fasta_files_as_dict(in_folder)
    
    genes_count = len(files_dict[in_files[0]])
    counts_all_identical_flag = genes_count_same_in_all_files(files_dict, in_files)

    if counts_all_identical_flag == True:
        sys.stderr.write("Processing " + str(len(files_dict)) + " genes from " + str(genes_count) + " samples\n")
        for i in range(0, genes_count):
            out_seq = ""
            for counter, fasta_file in enumerate(files_dict):
                if counter == 0:
                    species_prefix = files_dict[fasta_file][i][0]
                    species_prefix = species_prefix.split("_")[0]
                    print(">" + species_prefix)
                out_seq = out_seq + files_dict[fasta_file][i][1]
            out_seq = out_seq.replace(" ", "")
            hpf.print_with_fixed_row_length(out_seq, 80)
    else:
        sys.stderr.write("Error: the number of proteins is not equal across all input alignments\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_folder", type=str, help="Path to folder with FASTA alignments that have been processed with Gblocks")
    args = parser.parse_args()
    main(args.in_folder)



