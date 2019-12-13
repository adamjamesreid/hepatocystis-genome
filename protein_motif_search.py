#!/usr/bin/env python3
"""
Script for detecting motifs in protein FASTA sequences.
Argument1: path to multiFASTA file. Argument2: motif sequence.
Format of the motif:
Uppercase letters correspond to required amino acids. x corresponds to a position that can contain any amino acid. Positions that contain an amino acid from a limited set of options
should be surrounded with square brackets and within the square brackets the possible amino acids should be separated by slashes.
Example: RxLx[E/Q/D]
This is a motif that has R in the first position, any amino acid in the second position, L in the third position, any amino acid in the fourth position and either E, Q or D in the fifth position.   
Output (STDOUT): tab separated file. Columns: 1) FASTA header of the sequence with motif hit, 2) motif 1-based start coordinate, 3) motif end coordinate.
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License
"""

import sys
import hep_project_shared_functions as hpf
import argparse

def motif_string_to_list(motif_string):
    """
    Converts the input string into a list where each element corresponds to an amino acid position in the motif. 
    The elements of the list can be either single characters or a lists of characters.
    """
    motif_list = list()
    square_brackets_flag = False
    square_bracket_content = list()
    for char in motif_string:
        if char == "[":
            square_brackets_flag = True
            square_bracket_content = list()
        elif char == "]":
            motif_list.append(square_bracket_content)
            square_bracket_content = list()
            square_brackets_flag = False
        elif square_brackets_flag == True and char != "/":
            square_bracket_content.append(char)
        elif square_brackets_flag == False:
            motif_list.append(char)
    return motif_list


def detect_motif_match(seq_fragment, motif_list):
    """
    Compares a protein sequence fragment string with the motif consensus that is encoded as a list. 
    Returns True if the sequence matches the motif and False if it does not match
    """
    match_flag = True
    for i in range(0, len(seq_fragment)):
        seq_char = seq_fragment[i]
        motif_entry = motif_list[i]
        if motif_entry != "x":
            if isinstance(motif_entry, str):
                if seq_char != motif_entry:
                    match_flag = False
                    break
            elif isinstance(motif_entry, list):
                if seq_char not in motif_entry:
                    match_flag = False
                    break
    return match_flag
                   

def check_for_motif(motif_list, fasta_header, seq):
    """
    Input: 1) motif consensus encoded as a list, 2) a header of a protein FASTA entry, 3) sequence of a protein FASTA entry
    Output: For every match to the motif that is found in the FASTA sequence, the FASTA header of the query and the coordinates and the sequence of the match are printed out
    """
    seq = seq.upper()
    seq_len = len(seq)
    motif_len = len(motif_list)
    if seq_len >= motif_len:
        for start_pos in range(0, seq_len - motif_len + 1):
            seq_fragment = seq[start_pos: start_pos + motif_len]
            match_flag = detect_motif_match(seq_fragment, motif_list)
            if match_flag == True:
                out_line = fasta_header + "\t" + str(start_pos + 1) + "\t" + str(start_pos + motif_len) + "\t" + seq_fragment
                print(out_line)


def main(fasta_path, query_motif):
    motif_list = motif_string_to_list(query_motif)
    print("FASTA_header\tmatch_start\tmatch_end\tmatch_sequence")
    fasta_data = hpf.read_fasta_in_chunks(fasta_path)
    for header, seq in fasta_data:
        check_for_motif(motif_list, header, seq)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path protein FASTA file")
    parser.add_argument("query_motif", type=str, help="Query motif sequence")
    args = parser.parse_args()
    main(args.fasta_path, args.query_motif)