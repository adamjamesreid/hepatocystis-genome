#!/usr/bin/env python3
"""
Script for finding the average number of variants per 10 kb of reference genome in a VCF file (derived using GATK by merging GVCF files of multiple samples)
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License
"""

import hep_project_shared_functions as hpf
import numpy as np
import argparse

def main(input_vcf_file, genome_size, exclude_samples):
    exclude_list = list()
    if exclude_samples != None:
        exclude_list = exclude_samples.split(" ")

    in_data = hpf.l(input_vcf_file)

    in_data = [n for n in in_data if n.startswith("#") == False or n.startswith("#CHROM") == True]
    header = in_data[0]
    split_header = header.split("\t")
    sample_names = split_header[9:len(split_header)]
    samples_count = len(sample_names)
    counts_dict = dict()
    for sample_name in sample_names:
        counts_dict[sample_name] = 0
    in_data = in_data[1:len(in_data)]

    print("Sample\tMean_nr_of_variants_per_10_kb")

    for line in in_data:
        split_line = line.split()
        sample_entries = split_line[9: len(split_line)]

        for i in range(0, samples_count):
            sample_name = sample_names[i]
            sample_entry = sample_entries[i]
            if sample_entry != ".:0,0:0:.:0,0":
                allele = sample_entry.split(":")[0]
                if allele != "0" and allele != ".":
                    counts_dict[sample_name] += 1

    freq_list = list()
    for sample in counts_dict:
        if sample not in exclude_list:
            sample_count = counts_dict[sample]
            sample_freq = sample_count / (genome_size / 10000)
            print(sample + "\t" + str(sample_freq))
            freq_list.append(sample_freq)

    print("-----")
    print("Averaged numbers across all samples")
    print("Mean number of variants per 10 kb:", np.mean(freq_list))
    print("Median number of variants per 10 kb:", np.median(freq_list))
    print("Standard deviation of the number of variants per 10 kb:", np.std(freq_list))
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_vcf_file", type=str, help="Path to input VCF file")
    parser.add_argument("genome_size", type=int, help="Reference genome size (in base pairs)")
    parser.add_argument("-e", "--exclude_samples", help="Samples to exclude (list of sample names surrounded by quotation marks, delimited by space. Example: \"SAMN07757853 SAMN07757863 SAMN07757870\")", type=str)
    args = parser.parse_args()
    main(args.input_vcf_file, args.genome_size, args.exclude_samples)
