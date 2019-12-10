#!/usr/bin/env python3
"""
Script for counting the number of variants per fixed length segments of the reference genome (default: 100 kb) in a VCF file using a sliding window.
Input: a VCF file derived using GATK by merging GVCF files of multiple samples. The script assumes that the assembly that the reads were mapped to for variant calling was concatenated into one pseudochromosome before mapping.
Output: a CSV file where the rows correspond to samples and the columns correspond to genome sequence bins.
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License
"""

import hep_project_shared_functions as hpf
import numpy as np
import pandas as pd
import math
from collections import defaultdict
import argparse


def main(input_vcf_file, out_path, window_step, exclude_samples):
    exclude_list = list()
    if exclude_samples != None:
        exclude_list = exclude_samples.split(" ")

    in_data = hpf.l(input_vcf_file)

    in_data = [n for n in in_data if n.startswith("#") == False or n.startswith("#CHROM") == True]
    header = in_data[0]
    split_header = header.split("\t")
    sample_names = split_header[9:len(split_header)]
    samples_count = len(sample_names)
    counts_dict = defaultdict(dict)
    in_data = in_data[1:len(in_data)]

    for line in in_data:
        split_line = line.split()
        sample_entries = split_line[9: len(split_line)]
        coord = int(split_line[1])
        bin = math.floor(coord / window_step)
        bin = "bin_" + str(bin).zfill(3)

        for i in range(0, samples_count):
            sample_name = sample_names[i]
            if sample_name not in exclude_list:
                sample_entry = sample_entries[i]
                if sample_name not in counts_dict[bin]:
                    counts_dict[bin][sample_name] = 0
                if sample_entry != ".:0,0:0:.:0,0":
                    allele = sample_entry.split(":")[0]
                    if allele != "0" and allele != ".":
                        counts_dict[bin][sample_name] += 1

    df = pd.DataFrame(counts_dict)
    df.to_csv(out_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_vcf_file", type=str, help="Path to input VCF file")
    parser.add_argument("out_path", type=str, help="Path for output CSV file")
    parser.add_argument("-w", "--window_step", type=int, help="Sliding window step size (bp), default: 100000", default=100000)
    parser.add_argument("-e", "--exclude_samples", help="Samples to exclude (list of sample names surrounded by quotation marks, delimited by space. Example: \"SAMN07757853 SAMN07757863 SAMN07757870\")", type=str)
    args = parser.parse_args()
    main(args.input_vcf_file, args.out_path, args.window_step, args.exclude_samples)
