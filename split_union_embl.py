#!/usr/bin/env python3
"""
Script for splitting a union EMBL file (that has been made from scaffolds EMBL file using 'EMBOSS union') back into EMBL files of individual scaffolds.
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation

import hep_project_shared_functions as hpf
import sys
import pandas as pd
import argparse
import re
import os


def simplify_feature_location_coords(feature):
    """
    Input: feature object of SeqIO
    Output: feature start and feature end as integers.
        BeforePosition and AfterPosition are simplified into a single nucleotide position
    """
    feature_start = None
    feature_end = None

    if type(feature.location.start) == SeqFeature.ExactPosition:
        feature_start = int(feature.location.start)
    elif type(feature.location.start) == SeqFeature.BeforePosition:
        feature_start = int(feature.location.start) + 1
    elif type(feature.location.start) == SeqFeature.AfterPosition:
        feature_start = int(feature.location.start) - 1

    if type(feature.location.end) == SeqFeature.ExactPosition:
        feature_end = int(feature.location.end)
    elif type(feature.location.end) == SeqFeature.BeforePosition:
        feature_end = int(feature.location.end) + 1
    elif type(feature.location.end) == SeqFeature.AfterPosition:
        feature_end = int(feature.location.end) - 1
    return feature_start, feature_end


def get_scaffold_coords_by_source_features(in_path):
    """
    Input: path to a union EMBL file where the locations of original scaffolds have been marked with source features
    Output: a pandas dataframe with columns for the start coordinate, end coordinate, header (which will be "scaffold_" and a number) 
        and strand of each scaffold in the union EMBL sequence
    """
    scaffs_collection = list()
    scaffs_counter = 0
    records = list(SeqIO.parse(in_path, "embl"))
    for feature in records[0].features:
        if feature.type == "source":
            scaff_dict = dict()
            if type(feature.location) == FeatureLocation:
                scaff_id = "scaffold_" + str(scaffs_counter).zfill(10)
                scaff_dict["header"] = scaff_id
                scaff_dict["id"] = scaff_id
                scaff_dict["start_coord"] = int(feature.location.start) + 1
                scaff_dict["end_coord"] = int(feature.location.end) + 1
                if int(feature.location.strand) == 1:
                    scaff_dict["strand"] = "+"
                else:
                    scaff_dict["strand"] = "-"
                scaffs_collection.append(scaff_dict)
                scaffs_counter += 1
    scaffs_df = pd.DataFrame(scaffs_collection)
    return scaffs_df


def get_scaffold_coords_by_fasta(union_embl_path, fasta_path):
    """
    Input: 1) path to a union EMBL file, 2) path to a FASTA file with scaffolds that were used to make the union EMBL file.
    This function uses string search to find the sequences from the FASTA file in the union EMBL file.
    Output: a pandas dataframe with columns for the start coordinate, end coordinate, header and strand of each scaffold in the union EMBL sequence 
    """
    records = list(SeqIO.parse(union_embl_path, "embl"))
    union_seq = str(records[0].seq)
    scaffs_data = list(hpf.read_fasta_in_chunks(fasta_path))
    scaffs_collection = list()

    scaffs_counter = 0
    for header, seq in scaffs_data:
        scaff_dict = dict()
        seq = seq.upper()
        strand = "+"
        hit_count = 0
        pos = union_seq.find(seq)
        if pos != -1:
            hit_count = union_seq.count(seq)
        elif pos == -1:
            pos = union_seq.find(hpf.reverse_complement(seq))
            if pos == -1:
                strand = "NA"
            else:
                strand = "-"
                hit_count = union_seq.count(seq)
        if hit_count > 1:
            print("Warning: sequence " + header + " found in more than 1 place in the union file. This sequence is excluded from output", file=sys.stderr)
        else:
            start_coord = -1
            end_coord = -1
            if pos != -1:
                start_coord = pos + 1
                end_coord = start_coord + len(seq)
            scaff_id = "scaffold_" + str(scaffs_counter).zfill(10)
            scaff_dict["header"] = header
            scaff_dict["id"] = scaff_id
            scaff_dict["start_coord"] = start_coord
            scaff_dict["end_coord"] = end_coord
            scaff_dict["strand"] = strand
            scaffs_collection.append(scaff_dict)
        scaffs_counter += 1

    scaffs_df = pd.DataFrame(scaffs_collection)
    return scaffs_df


def process_record_features(records, coords_df, query_start_coord, query_end_coord, my_sequence_record):
    """
    Looks for a feature between query_start_coord and query_end_coord in the union EMBL record and adds it to the new EMBL record (with recalculated coordinates) if found
    """
    for feature in records[0].features:
        
        feature_start, feature_end = simplify_feature_location_coords(feature)
        
        if feature_start >= query_start_coord and feature_end <= query_end_coord: 
            if "organism" not in feature.qualifiers:
                scaffold_gap_entry = False
                if "note" in feature.qualifiers:
                    if feature.qualifiers["note"] == ["*gap_type: between scaffolds"]:
                        scaffold_gap_entry = True
                if scaffold_gap_entry == False:
                    
                    codon_start_value = None
                    
                    if "codon_start" in feature.qualifiers:
                        codon_start_value = int(feature.qualifiers["codon_start"][0])

                    # Conversion of coordinates in union EMBL to local coordinates in a scaffold
                    scaff_row = coords_df[(coords_df["start_coord"] <= feature_start) & (coords_df["end_coord"] >= feature_end)]
                    scaff_coord_start = feature_start - int(scaff_row["start_coord"]) + 1

                    start_shift = scaff_coord_start - feature_start
            
                    my_feature_location = feature.location._shift(start_shift)
                    my_feature_location.strand = feature.location.strand
                    my_feature = feature
                    my_feature.location = my_feature_location
                    my_feature.qualifiers = feature.qualifiers

                    if codon_start_value != None:
                        my_feature.qualifiers["codon_start"] = [str(codon_start_value)]                           
                    my_sequence_record.features.append(my_feature)
    return my_sequence_record



def main(in_path, out_folder, fasta_path, deselected_scaffolds_path):
    coords_df = None
    deselected_scaffolds = []

    if fasta_path == "":
        coords_df = get_scaffold_coords_by_source_features(in_path)
    else:
        coords_df = get_scaffold_coords_by_fasta(in_path, fasta_path)
        if deselected_scaffolds_path != "":
            deselected_scaffolds = hpf.l(deselected_scaffolds_path)

    records = list(SeqIO.parse(in_path, "embl"))
    t = os.system("mkdir -p " + out_folder)
    if t != 0:
        sys.stderr.write("Error occurred when checking for the presence of output folder or creating the output folder ()" + out_folder + ")\n")
        sys.exit(1)

    for selected_scaff in range(0, coords_df.shape[0]):
        coords_df_entry = coords_df.iloc[selected_scaff]
        scaff_name = coords_df_entry["header"]
        scaff_id = coords_df_entry["id"]
        my_sequence_record = None

        if scaff_name not in deselected_scaffolds:
            query_start_coord = int(coords_df_entry.start_coord)
            query_end_coord = int(coords_df_entry.end_coord)
            out_path = out_folder + "/" + scaff_id + ".embl"

            union_seq = str(records[0].seq)
            seq = union_seq[query_start_coord - 1:query_end_coord -1]
            my_sequence = Seq(seq)
            my_sequence_record = SeqRecord(my_sequence, id=scaff_id, name=scaff_name, description="unknown_description", dbxrefs=[])
            my_sequence_record.seq.alphabet = generic_dna
            my_sequence_record.accession = "unknown_accession"

            my_sequence_record = process_record_features(records, coords_df, query_start_coord, query_end_coord, my_sequence_record)

        SeqIO.write(my_sequence_record, out_path, "embl")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_path", help="Path to union EMBL file", type=str)
    parser.add_argument("out_folder", type=str, help="Path to output folder")
    parser.add_argument("--fasta_path", help="Optional: path to a FASTA file that contains the original scaffolds that were used to produce the union EMBL file.\
        If this input is not provided, this script looks for source features in the union EMBL file to determine where the borders between the original scaffolds are.", type=str, default="")
    parser.add_argument("--deselected_scaffolds_path", type=str, help="Optional: path to a text file where each line contains the name of a scaffold that should be excluded from output. This only works when using the --coords_path option.", default="")
    args = parser.parse_args()
    main(args.in_path, args.out_folder, args.fasta_path, args.deselected_scaffolds_path)
