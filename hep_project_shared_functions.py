#!/usr/bin/env python3
"""
File for functions that are shared between many scripts in the Hepatocystis project
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License
"""

import csv
import os
from os import listdir
from os.path import isfile, join
import sys


def get_file_paths(in_folder_path, extension):
    """
    Function for getting the paths to all files with a specific extension in a user-specified folder
    in_folder_path: path to the folder with input files
    extension: file extension of input files
    Output: paths to individual files with the specific extension (list)
    """
    onlyfiles = list()
    selected_file_paths = list()
    if os.path.isdir(in_folder_path):
        onlyfiles = [f for f in listdir(in_folder_path) if isfile(join(in_folder_path, f))]
        for file_item in onlyfiles:
            if "." + extension in file_item:
                file_item_split = file_item.split(".")
                if file_item_split[len(file_item_split) - 1] == extension:
                    selected_file_paths.append(in_folder_path + "/" + file_item)
    else:
        sys.stderr.write("Error: folder not found (" + in_folder_path + ")\n")
        sys.exit(1)
    return selected_file_paths


def l(path):
    """
    Function for loading text file as a list and removing line breaks from line ends
    """
    lines = []
    if isfile(path):
        with open(path, "r") as in_file:
            lines = in_file.readlines()
            lines = [x.rstrip() for x in lines]
    else:
        sys.stderr.write("Error: file not found (" + path + ")\n")
        sys.exit(1)  
    return lines


def ll(in_path):
    """
    Function for reading a text file line by line
    """
    if isfile(in_path):
        with open(in_path, "r") as in_file:
            for line in in_file:
                line = line.rstrip()
                yield line
    else:
        sys.stderr.write("Error: file not found (" + in_path + ")\n")
        sys.exit(1)


def export_list_as_comma_separated_file(out_list, out_path):
    """
    Exports a list as a comma separated file, without line breaks between items
    """
    with open(out_path, "w") as out_file:
        for counter, item in enumerate(out_list):
            out_file.write(str(item))
            if counter < len(out_list)-1:
                out_file.write(",")


def export_list_as_line_break_separated_file(out_list, out_path):
    """
    Exports a list to a file, each item on a separate row
    """
    with open(out_path, "w") as out_file:
        for item in out_list:
            out_file.write(str(item))
            out_file.write("\n")


def split_list_to_sublists(in_list, separator):
    """
    Input: a list of strings and a separator string
    Output: list of lists (the input list split into shorter lists). Every time the separator string is found in an input list item, this item and the ones following it are
        moved to a new sublist
    """
    sublists = list()
    current_sublist = list()
    for item in in_list:
        if separator in item:
            if len(current_sublist) > 0:
                sublists.append(current_sublist)
                current_sublist = list()
        current_sublist.append(item)
    if len(current_sublist) > 0:
        sublists.append(current_sublist)
    return sublists


def spl(line, left_splitter, right_splitter, direction=0):
    """
    Function for cropping a string from left and right side
    Direction:  if 0: the string will be cropped first from left and then right
                if 1: the string will be cropped first from right and then left
    Returns None if the splitting cannot be done because a splitter or both splitters are not in the input string
    """
    out_line = None
    if left_splitter in line and right_splitter in line:
        if direction == 0:
            out_line = line.split(left_splitter)[1]
            out_line = out_line.split(right_splitter)[0]
        elif direction == 1:
            out_line = line.split(right_splitter)[0]
            out_line = out_line.split(left_splitter)[1]
    return out_line


def print_with_fixed_row_length(seq, max_length):
    """
    Input: 1) a string 2) maximum line length in output
    Output: the input string printed to STDOUT in chunks with fixed maximum line length
    """
    split_seq = [seq[i:i+max_length] for i in range(0, len(seq), max_length)]
    for line in split_seq:
        print(line)


def split_with_fixed_row_length(seq, max_length):
    """
    Input: 1) a string 2) maximum line length in output
    Output: the input string split in chunks with fixed maximum line length
    """
    split_seq = [seq[i:i + max_length] for i in range(0, len(seq), max_length)]
    return split_seq


def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    reverse_comp = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_comp


def read_fasta_in_chunks(in_path):
    """
    Input: path to FASTA file
    Output (iterator): string tuples where the first element is a FASTA header and the second element is the corresponding FASTA sequence
    """
    in_data = ll(in_path)
    current_seq_header = None
    seq = ""
    for line in in_data:
        if line != "":
            if line[0] == ">":
                if seq != "":
                    yield (current_seq_header, seq)
                seq = ""
                current_seq_header = line[1:len(line)]
            else:
                seq += line
    if seq != "":
        yield (current_seq_header, seq)
