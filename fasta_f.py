#!/usr/bin/env python3
"""
Multitool script for processing FASTA files. It collects various functions for performing different operations: breaking scaffolds into contigs, filtering sequences
by length, finding the frequency of stop codons, checking the completeness of transcripts, extracting or removing sequences by name, batch editing of headers,
getting the GC and length distribution of sequences and truncating or deduplicating the sequences.
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License
"""

import hep_project_shared_functions as hpf
import argparse
import re
import sys

def add_quotes_to_header_with_commas(header):
    """
    Surrounds the FASTA header with quotation marks if the header contains commas
    """
    if "," in header:
            header = "\"" + header + "\""
    return header


def print_header_and_seq(header, seq):
    """
    Prints the FASTA header and sequence
    """
    print(">" + header)
    hpf.print_with_fixed_row_length(seq, 80)


def scaffs_to_contigs(args):
    """
    Function for breaking FASTA scaffolds into contigs
    """
    for header, seq in args.fasta_data:
        seq = seq.upper()
        split_seq = seq.split("N")
        header = header.replace(" ", "_")
        subseq_counter = 0
        for subseq in split_seq:
            if subseq != "":
                subseq_header = header + "_" + str(subseq_counter).zfill(5)
                print_header_and_seq(subseq_header, subseq)
                subseq_counter += 1


def filter_fasta_by_length(args):
    """
    Filters sequences in a FASTA file by length. The default mode is to remove sequences shorter than the cutoff value.
    In the low_pass mode, sequences longer than the cutoff value will be removed.
    """
    for header, seq in args.fasta_data:
        seq_len = len(seq)
        seq_len_ok_flag = True
        if args.low_pass == True:
            if seq_len > args.cutoff:
                seq_len_ok_flag = False
        else:
            if seq_len < args.cutoff:
                seq_len_ok_flag = False
        if seq_len_ok_flag == True:
            print_header_and_seq(header, seq)


def get_stop_codon_freq(args):
    """
    Finds the frequency of stop codons (count per sequence base pair) for each sequence in a multiFASTA file
    """
    STOP_CODON_SEQ = ["TAG", "TAA", "TGA"]
    for header, seq in args.fasta_data:
        stop_codons_count = 0
        for stop_codon in STOP_CODON_SEQ:
            query = "(?=" + stop_codon + ")"
            seq = seq.upper()
            rev_comp_seq = hpf.reverse_complement(seq)
            stop_codons_count += len(re.findall(query, seq))
            stop_codons_count += len(re.findall(query, rev_comp_seq))
            stop_codon_freq = stop_codons_count / len(seq)
        header = add_quotes_to_header_with_commas(header)
        print(header + "," + str(stop_codon_freq))


def start_and_stop_codon_presence(args):
    """
    Checks the completeness of transcripts in a DNA FASTA file of transcripts by determining if the sequence starts with a start codon
    and ends with a stop codon. Output: CSV table. Column 1: FASTA headers. Column 2: presence of start codon (Y/N).
    Column 3: presence of stop codon (Y/N).
    """
    print("sequence,start_codon_present,stop_codon_present")
    for header, seq in args.fasta_data:
        seq = seq.upper()
        start_codon_flag = False
        stop_codon_flag = False
        seq = seq.upper()
        if seq.startswith("ATG"):
            start_codon_flag = True
        if (start_codon_flag == True and len(seq) % 3 == 0) or start_codon_flag == False:
            if seq.endswith("TAG") or seq.endswith("TAA") or seq.endswith("TGA"):
                stop_codon_flag = True
        header = add_quotes_to_header_with_commas(header)
        out_line = header + ","
        if start_codon_flag == True:
            out_line += "Y"
        else:
            out_line += "N"
        out_line += ","
        if stop_codon_flag == True:
            out_line += "Y"
        else:
            out_line += "N"
        print(out_line)


def replace_fasta_headers_with_string_plus_sequential_numbers(args):
    """
    Replaces FASTA headers with a string and sequential numbers
    """
    for counter, item in enumerate(args.fasta_data):
        seq = item[1]
        header = args.prefix + "_" + str(counter).zfill(10)
        print_header_and_seq(header, seq)
        

def extract_sequences_from_fasta_by_id(args):
    """
    Function for extracting sequences from a FASTA file by their names. The sequence names are truncated at the first space character
    before they are compared to the query string. The function allows extracting sequences based on 1 query string and also reading a list of 
    query strings from a text file. There is an 'invert' mode to extract all sequences that do not match the query string(s).
    """
    selected_seq_list = None
    if args.string_query == True:
        selected_seq_list = [args.query]
    else:
        selected_seq_list = hpf.l(args.query)
    for header, seq in args.fasta_data:
        fasta_seq_id = header.split()[0]
        seq_to_output = False
        if args.invert == False:
            if fasta_seq_id in selected_seq_list:
                seq_to_output = True
        else:
            if fasta_seq_id not in selected_seq_list:
                seq_to_output = True
        if seq_to_output == True:
            print_header_and_seq(header, seq)


def seq_lengths(args):
    """
    Prints the FASTA headers and lengths of the corresponding sequences in a FASTA file
    """
    for header, seq in args.fasta_data:
        header = add_quotes_to_header_with_commas(header)
        print(header + "," + str(len(seq)))


def gc_percentages(args):
    """
    Prints the FASTA headers and the GC% of the corresponding sequences in a FASTA file
    """
    for header, seq in args.fasta_data:
        counts_dict = {"A": 0, "T": 0, "G": 0, "C": 0}
        header = add_quotes_to_header_with_commas(header)
        seq = seq.upper()
        for char in seq:
            if char in counts_dict:
                counts_dict[char] += 1
        counts_sum = sum(counts_dict.values())
        gc_percentage = None
        if counts_sum > 0:
            gc_percentage = ((counts_dict["G"]+counts_dict["C"])/counts_sum) * 100
        print(header + "," + str(gc_percentage))


def average_gc_percentage(args):
    """
    Finds the average GC% of the FASTA file
    """
    counts_dict = {"A": 0, "T": 0, "G": 0, "C": 0}
    for item in args.fasta_data:
        seq = item[1]
        seq = seq.upper()
        for char in seq:
            if char in counts_dict:
                counts_dict[char] += 1
    counts_sum = sum(counts_dict.values())
    gc_percentage = None
    if counts_sum > 0:
        gc_percentage = ((counts_dict["G"]+counts_dict["C"])/counts_sum) * 100
    print("Average GC percentage of the FASTA file: " + str(gc_percentage))


def truncate_seq(args):
    """
    Truncates each sequence in a FASTA file at a specfic length cutoff
    """
    for header, seq in args.fasta_data:
        if len(seq) > args.max_len:
            seq = seq[0:args.max_len]
        print_header_and_seq(header, seq)


def shorten_headers(args):
    """
    Truncates FASTA headers by splitting them with a delimiter (default: space) and keeping the first element
    """
    for header, seq in args.fasta_data:
        header = header.split(args.splitter)[0]
        print_header_and_seq(header, seq)


def add_string_to_fasta_headers(args):
    """
    Appends a string to all FASTA headers. There is an option to prepend the string instead of appending
    """
    for header, seq in args.fasta_data:
        if args.prepend == False:
            header += "_" + args.new_string
        else:
            header = args.new_string + "_" + header
        print_header_and_seq(header, seq)


def replace_problematic_chars_in_headers(args):
    """
    Sanitises FASTA headers by replacing specific characters with underscores
    """
    bad_text_list = [" ", ":", "#", ",", "-", "|", "'", "[", "]", "(", ")", "*", "__", "___"]
    for header, seq in args.fasta_data:
        for query in bad_text_list:
            header = header.replace(query, "_")
        print_header_and_seq(header, seq)


def remove_identical_sequences(args):
    """
    Removes identical sequences from FASTA. By default, sequences that are the reverse complement of one another are considered to be identical.
    If the fwd_orientation_only flag is used, the function compares the sequences only in FWD orientation
    """
    checked_seq = list()
    out_list = list()
    for header, seq in args.fasta_data:
        seq_upper = seq.upper()
        seq_ok_flag = False
        if args.fwd_orientation_only == True:
            if seq_upper not in checked_seq:
                seq_ok_flag = True
        else:
            if seq_upper not in checked_seq and hpf.reverse_complement(seq_upper) not in checked_seq:
                seq_ok_flag = True
        if seq_ok_flag == True:
            out_list.append((">" + header, seq))
            checked_seq.append(seq_upper)
    for header, seq in out_list:
        print_header_and_seq(header, seq)


def positive_int(value):
    """
    Function for allowing only positive integer values in an Argparse argument
    """
    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(value + " is not a valid positive integer value")
    if value <= 0:
        raise argparse.ArgumentTypeError(value + " is not a valid positive integer value")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers()

    parser_scaffs_to_contigs = subparsers.add_parser("scaffs_to_contigs", description=scaffs_to_contigs.__doc__)
    parser_scaffs_to_contigs.set_defaults(func=scaffs_to_contigs)

    parser_filter_fasta_by_length = subparsers.add_parser("filter_fasta_by_length", description=filter_fasta_by_length.__doc__)
    parser_filter_fasta_by_length.add_argument("cutoff", type=positive_int, help="Sequence length cutoff")
    parser_filter_fasta_by_length.add_argument("-l", "--low_pass", default=False, action="store_true", help="Optional: low_pass mode. In low_pass mode, sequences longer than the cutoff will be removed. Otherwise, by default, sequences shorter than the cutoff will be removed")
    parser_filter_fasta_by_length.set_defaults(func=filter_fasta_by_length)

    parser_get_stop_codon_freq = subparsers.add_parser("get_stop_codon_freq", description=filter_fasta_by_length.__doc__)
    parser_get_stop_codon_freq.set_defaults(func=get_stop_codon_freq)

    parser_start_and_stop_codon_presence = subparsers.add_parser("start_and_stop_codon_presence", description=start_and_stop_codon_presence.__doc__)
    parser_start_and_stop_codon_presence.set_defaults(func=start_and_stop_codon_presence)

    parser_replace_headers = subparsers.add_parser("replace_headers", description=replace_fasta_headers_with_string_plus_sequential_numbers.__doc__)
    parser_replace_headers.add_argument("--prefix", type=str, default="seq", help="Prefix string that will be used in the FASTA headers. Default: 'seq'")
    parser_replace_headers.set_defaults(func=replace_fasta_headers_with_string_plus_sequential_numbers)

    parser_extract_seq = subparsers.add_parser("extract_seq", description=extract_sequences_from_fasta_by_id.__doc__)
    parser_extract_seq.add_argument("query", type=str)
    parser_extract_seq.add_argument("-v", "--invert", default=False, action="store_true", help="Optional argument. If the 'invert' flag is used, the output will be all sequences that do not match the query sequence name(s)")
    parser_extract_seq.add_argument("-s", "--string_query", default=False, action="store_true", help="Optional: a single query can be specified as an argument here instead of loading a list queries from a file")
    parser_extract_seq.set_defaults(func=extract_sequences_from_fasta_by_id)

    parser_seq_lengths = subparsers.add_parser("seq_lengths", description=seq_lengths.__doc__)
    parser_seq_lengths.set_defaults(func=seq_lengths)

    parser_gc_percentages = subparsers.add_parser("gc_percentages", description=gc_percentages.__doc__)
    parser_gc_percentages.set_defaults(func=gc_percentages)

    parser_average_gc = subparsers.add_parser("average_gc", description=average_gc_percentage.__doc__)
    parser_average_gc.set_defaults(func=average_gc_percentage)

    parser_truncate_seq = subparsers.add_parser("truncate_seq", description=truncate_seq.__doc__)
    parser_truncate_seq.add_argument("max_len", type=positive_int, help="Length cutoff in sequence truncation")
    parser_truncate_seq.set_defaults(func=truncate_seq)

    parser_shorten_headers = subparsers.add_parser("shorten_headers", description=shorten_headers.__doc__)
    parser_shorten_headers.add_argument("--splitter", type=str, default=" ", help="Delimiter that will be used when splitting the headers (default: space)")
    parser_shorten_headers.set_defaults(func=shorten_headers)

    parser_add_string_to_headers = subparsers.add_parser("add_string_to_headers", description=add_string_to_fasta_headers.__doc__)
    parser_add_string_to_headers.add_argument("new_string", type=str, help="String that will be added to FASTA headers")
    parser_add_string_to_headers.add_argument("-p", "--prepend", default=False, action="store_true", help="Optional: prepend the string to the FASTA headers instead of appending it")
    parser_add_string_to_headers.set_defaults(func=add_string_to_fasta_headers)

    parser_replace_chars = subparsers.add_parser("replace_chars", description=replace_problematic_chars_in_headers.__doc__)
    parser_replace_chars.set_defaults(func=replace_problematic_chars_in_headers)

    parser_remove_identical_sequences = subparsers.add_parser("dedupe", description=remove_identical_sequences.__doc__)
    parser_remove_identical_sequences.add_argument("-f", "--fwd_orientation_only", default=False, action="store_true", help="Optional: compare the sequences in forward orientation only")
    parser_remove_identical_sequences.set_defaults(func=remove_identical_sequences)

    parser.add_argument("fasta_path", type=str)
    args = parser.parse_args()
    args.fasta_data = hpf.read_fasta_in_chunks(args.fasta_path)
    args.func(args)
    
if __name__ == "__main__":
    main()
