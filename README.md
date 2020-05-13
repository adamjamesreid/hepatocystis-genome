# Hepatocystis genome (Aunin et al.)
Data and code relating to the Hepatocystis ex. Piliocolobus tephrosceles genome and transcriptome paper







## Primary genome assembly and annotation files
PRJEB32891_scaffolds.fasta.gz - assembly sequence in fasta format

PRJEB32891_union.embl.gz - assembly and annotation in single EMBL record (union)

PRJEB32891_scaffolds.embl.gz - assembly and annotation in EMBL format

PRJEB32891_proteins.faa - protein sequences of predicted genes in fasta format

PRJEB32891_transcripts.fa - spliced nucleotide sequences of predicted genes in fasta format

## Alignments for phylogenetic trees
S1_dataset_genes_for_phylogenetic_trees.xlsx - Table of the IDs of sequences that were used to generate the phylogenetic trees

apicoplast_concat_alignments.faa - alignment of apicoplast sequences

cytochrome_b_alignment.fa - cytochrome b alignment

11_genes_with_hepatocystis_epomophori_concat_alignments.fa - 11 nuclear gene alignment

mitoch_proteins_concat_alignments.faa - alignment of mitochondrial sequences

nuclear_genome_proteins_concat_alignments.faa - alignment of nuclear genes

concatenate_fasta_alignments.py - Script for concatenating protein FASTA alignments of different genes in order to make a species phylogeny tree

## Deconvolution of Hepatocystis bulk RNA-seq using Malaria Cell Atlas data and CIBERSORT
merge_htseq-count_output_files.py - Script for merging htseq-count output files of multiple samples from same study into one table

MCA_pseudobulk.R - R code describing how to generate stage-specific pseudobulk samples as a reference for bulk RNA-seq deconvolution

generate_mixtures.py - generate mixtures of pseudobulk life stages at known percentages to test CIBERSORT for accuracy

mca_pseudobulk_meroringschizfilt_cpm.dat - Output from running generate_mixtures.py, e.g. pseudobulk for MCA life stages

## Examination of missing genes in Hepatocystis relative to Plasmodium
heps_per_mca_cluster.py - find orthologue groups shared between P. berghei and either P. ovale or P. vivax, but absent from Hepatocystis and see whether Malaria Cell Atlas gene clusters are enriched for these.

```python heps_per_mca_cluster.py hepatocystis_orthomcl.out mca_gene_clusters.dat Pberghei.prot.desc```

hepatocystis_orthomcl.out - orthoMCL clusters for protein sequences from Hepatocystis DNA assembly (Hepatocystis_DNA), RNA assembly (Hepatocystis_RNA) and various Plasmodium species.

## Evolutionary analysis of genes
codeml_batch.py - Wrapper script for running codeml (http://envgen.nox.ac.uk/bioinformatics/docs/codeml.html) as batch.
Requires codeml to be installed and in path.

gatk_count_variants_per_sample.py - Script for finding the average number of variants per 10 kb of reference genome in a VCF file (derived using GATK by merging GVCF files of multiple samples)

gatk_count_variants_per_sample_sliding_window.py - Script for counting the number of variants per fixed length segments of the reference genome (default: 100 kb) in a VCF file using a sliding window.
Input: a VCF file derived using GATK by merging GVCF files of multiple samples. The script assumes that the assembly that the reads were mapped to for variant calling was concatenated into one pseudochromosome before mapping.
Output: a CSV file where the rows correspond to samples and the columns correspond to genome sequence bins.

mca_gene_clusters.dat - File describing which genes are in which Malaria Cell Atlas clusters

hep_pberghei_povale_3-way_codeml_dn_results.txt - Results of codeml dN analysis

hepatocystis_dn_and_mca_clusters.Rmd - This R Markdown file contains the R code for generating Figure S11 from the manuscript on the draft genome of Hepatocystis sp. ex Piliocolobus tephrosceles

crop_translatorx_alignments.py - Script for processing TranslatorX alignments to remove badly aligned parts (similarly to what Gblocks does)

hep_pfam_domains_in_top_dn.py	- Script for checking for the enrichment of PFAM domains in the Hepatocystis proteins with top dN values

## Analysis of *Haemoproteus tartakovskyi* genome

hepatocystis_orthomcl_with_haemoproteus.out	- Output file of an OrthoMCL run that includes the proteome of Haemoproteus tartakovskyi in addition to the proteomes of Hepatocystis and Plasmodium species

htartakovskyi_proteins_companion.faa - Proteins of Haemoproteus tartakovskyi, annotated using Companion software 

## Processing of 10X Chromium sequencing reads

chromium_reads_remove_barcodes.cpp - C++ program to remove barcodes and linkers from FASTQ files of Chromium reads

chromium_reads_barcode_frequencies.cpp - C++ program for finding frequencies of barcodes in FASTQ files of Chromium reads

chromium_reads_extract_barcodes.cpp - C++ program for extracting Chromium barcode sequences of Chromium reads that have been selected for assembly

## Sequence analysis
fasta_f.py - Multitool script for processing FASTA files. It collects various functions for performing different operations: breaking scaffolds into contigs, filtering sequences by length, finding the frequency of stop codons, checking the completeness of transcripts, extracting or removing sequences by name, batch editing of headers, getting the GC and length distribution of sequences and truncating or deduplicating the sequences

protein_motif_search.py - Script for detecting motifs in protein FASTA sequences

split_union_embl.py - Script for splitting a union EMBL file (that has been made from scaffolds EMBL file using 'EMBOSS union') back into EMBL files of individual scaffolds

detect_low_complexity_contaminant_contigs.R	- Script that was used to detect low complexity contaminant contigs in the Hepatocystis genome assembly

hep_proteins_pfam_domains.txt	- List of PFAM domains in Hepatocystis proteins

## Other
hep_project_shared_functions.py - File for functions that are shared between many scripts in the Hepatocystis project



