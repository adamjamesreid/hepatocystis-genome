#!/usr/bin/env python3
"""
Script for merging htseq-count output files of multiple samples from same study into one table.
Argument 1: path to folder with htseq-count output files (with .txt extension).
Argument 2: path for output file.
Output: tab separated table made from input data.
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License
"""

import hep_project_shared_functions as hpf
import sys
import pandas as pd
import csv
import argparse


def get_all_gene_names(in_files):
    """
    Input: path to a folder with htseq-count output files
    Output: list of unique gene names in input files
    """
    genes_list = list()
    for file_counter, in_file in enumerate(in_files):
        in_data = hpf.l(in_file)
        file_gene_counter = 0
        for line in in_data:
            if line[0:2] != "__":
                split_line = line.split("\t")
                gene_name = split_line[0]
                if gene_name not in genes_list:
                    if file_counter == 0:
                        genes_list.append(gene_name)
                    else:
                        sys.stderr.write("Error: gene name " + gene_name + " in file " + in_file + " mismatches gene names in previously read file(s) in the same folder\n")
                        sys.exit(1)
                file_gene_counter += 1
        if file_gene_counter != len(genes_list):
            sys.stderr.write("Error: mismatch in the number of genes between input files (encountered when reading " + in_file + ")\n")
            sys.exit(1)
    return genes_list


def get_sample_name_from_file_name(file_name):
    """
    Input: input file name
    Output: sample name, extracted from file name
    """
    sample_name = file_name.split("/")[-1]
    sample_name = sample_name.split(".txt")[0]
    return sample_name


def main(in_folder, out_path):
    in_files = hpf.get_file_paths(in_folder, "txt")
    if len(in_files) == 0:
        sys.stderr.write("Error: no files with .txt extension was found in the input folder\n")
        sys.exit(1)

    genes_list = get_all_gene_names(in_files)

    collection_dict = dict()
    for in_file in in_files:
        genes_dict = dict()
        sample_name = get_sample_name_from_file_name(in_file)

        for item in genes_list:
            genes_dict[item] = 0

        in_data = hpf.l(in_file)
        for line in in_data:
            if line[0:2] != "__":
                split_line = line.split("\t")
                count = int(split_line[1])
                gene_name = split_line[0]
                genes_dict[gene_name] = count
        collection_dict[sample_name] = genes_dict
    collection_df = pd.DataFrame(collection_dict)
    collection_df.to_csv(out_path, sep="\t", quoting=csv.QUOTE_NONE)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_folder", type=str, help="path to folder with htseq-count output files (with .txt extension)")
    parser.add_argument("out_path", type=str, help="path for output file (tab separated table)")
    args = parser.parse_args()
    main(args.in_folder, args.out_path)




