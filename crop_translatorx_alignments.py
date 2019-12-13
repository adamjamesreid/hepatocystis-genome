#!/usr/bin/env python3
"""
Script for processing TranslatorX alignments to remove badly aligned parts (similarly with what Gblocks does).
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License
"""

import hep_project_shared_functions as hpf
import sys
import argparse
import os


def separate_seq_to_subseq(seq):
    """
    Input: DNA sequence (string, with - for gaps)
    Output: input string broken into a list of fragments, so that each list element contains either only gaps or only nucleotides
    """
    seq_list = list()
    current_str = ""
    previous_ch = None
    for i in range(0, len(seq)):
        ch = seq[i]
        if previous_ch != None:
            if ch == "-" and previous_ch != "-" or ch != "-" and previous_ch == "-":
                seq_list.append(current_str)
                current_str = ""
        current_str += ch
        previous_ch = ch
    if len(current_str) != 0:
        seq_list.append(current_str)
    return seq_list


def mask_short_alignment_pieces(seq_list, min_length):
    """
    Input: 1) output of separate_seq_to_subseq function, 2) minimum length of a alignment piece
    Output: input sequence with short alignment fragments replaced with % symbols
    """
    seq_rebuilt = ""
    for item in seq_list:
        if item[0] == "-":
            seq_rebuilt += item
        else:
            if len(item) >= min_length:
                seq_rebuilt += item
            else:
                seq_rebuilt += "%" * len(item)
    return seq_rebuilt


def check_for_in_frame_stop_codons(seq, header):
    """
    Prints a warning if in frame stop codons are detected
    """
    STOP_CODON_SEQ = ("TAG", "TAA", "TGA")
    pos = 0
    seq_len = len(seq)
    stop_codon_count = 0
    if seq_len >= 6:
        while pos < seq_len - 6:
            codon = seq[pos:pos + 3]
            if codon in STOP_CODON_SEQ:
                stop_codon_count += 1
            pos += 3
    if stop_codon_count > 0:
        sys.stderr.write(str(stop_codon_count) + " stop codon(s) found in frame in sequence " + header + "\n")


def load_alignment(fasta_path, min_fragment_length):
    """
    Input: 1) path to an alignment of two transcript sequences (in DNA FASTA format), 2) minimum length of alignment fragments between gaps
    Output: 1) FASTA headers (only the first part before whitespaces is kept), 2) alignment sequences where characters that are not A, T, G, or C have been removed
    """
    fasta_data = list(hpf.read_fasta_in_chunks(fasta_path))
    fasta_collection = list()
    masked_seq_list = list()

    seq_len = None
    for header, seq in fasta_data:
        if seq_len == None:
            seq_len = len(seq)
        else:
            if len(seq) != seq_len:
                sys.stderr.write("FASTA alignment warning: the length of sequence " + header + " in the file '" + fasta_path + "' does not match the length of the other sequence(s) \
                    in the alignment. The sequences will get cropped to the length of the shortest sequence in the alignment\n")
        fasta_dict = dict()
        header = header.split()[0]
        header = header.replace(",", "_")
        fasta_dict["header"] = header
        seq = seq.upper()
        seq_list = separate_seq_to_subseq(seq)
        seq_masked = mask_short_alignment_pieces(seq_list, min_fragment_length)
        fasta_dict["cropped_seq"] = ""
        fasta_collection.append(fasta_dict)
        masked_seq_list.append(seq_masked)

    zipped_masked_seq_list = list(zip(*masked_seq_list))
    for row_tuple in zipped_masked_seq_list:
        row_string = "".join(row_tuple)
        row_ok_flag = False
        if "%" not in row_string and "-" not in row_string:
            row_ok_flag = True
        if row_ok_flag == True:
            for i in range(0, len(row_string)):
                fasta_collection[i]["cropped_seq"] = fasta_collection[i]["cropped_seq"] + row_string[i]

    for fasta_dict in fasta_collection:
        fasta_dict["cropped_seq"] = trim_excess_bases(fasta_dict["cropped_seq"])
        check_for_in_frame_stop_codons(fasta_dict["cropped_seq"], fasta_dict["header"])
    return fasta_collection


def trim_excess_bases(segment):
    """
    If a sequence length is not divisible by 3, the sequence is cropped to make the length divisible by 3
    """
    excess_bases = len(segment) % 3
    if excess_bases % 3 != 0:
        if len(segment) - excess_bases > 1:
            segment = segment[0:len(segment) - excess_bases]
    return segment


def main(alignments_folder, out_folder, min_fragment_length):
    
    os_command = "mkdir -p " + out_folder
    t = os.system(os_command)
    if t != 0:
        sys.exit("Error running command: " + os_command)
    
    in_files = os.listdir(alignments_folder)

    for in_file in in_files:
        out_list = list()
        fasta_path = alignments_folder + "/" + in_file
        
        split_fasta_path = in_file.split(".")
        out_file_path = ".".join(split_fasta_path[0:len(split_fasta_path) - 1]) + "_cropped.fa"
        out_file_path = out_folder + "/" + out_file_path
        
        fasta_collection = load_alignment(fasta_path, min_fragment_length)
        for fasta_dict in fasta_collection:
            out_list.append(">" + fasta_dict["header"])
            out_list.extend(hpf.split_with_fixed_row_length(fasta_dict["cropped_seq"], 80))
        hpf.export_list_as_line_break_separated_file(out_list, out_file_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("alignments_folder", type=str, help="Path to alignment file of transcript sequences (in DNA FASTA format)")
    parser.add_argument("out_folder", type=str, help="Path to folder for cropped alignment files", default="")
    parser.add_argument("--min_fragment_length", type=int, help="Minimum length of alignment fragments between gaps (fragments shorter than this will be left out from alignment). Default: 42", default=42)
    args = parser.parse_args()
    main(args.alignments_folder, args.out_folder, args.min_fragment_length)
        
