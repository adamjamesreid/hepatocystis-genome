# Hepatocystis genome (Aunin et al.)
Data and code relating to the Hepatocystis ex. Piliocolobus tephrosceles genome and transcriptome paper

## Primary genome assembly and annotation files
PRJEB32891_scaffolds.fasta.gz - assembly sequence in fasta format

PRJEB32891_union.embl.gz - assembly and annotation in single EMBL record (union)

PRJEB32891_scaffolds.embl.gz - assembly and annotation in EMBL format

PRJEB32891_proteins.faa - protein sequences of predicted genes in fasta format

PRJEB32891_transcripts.fa - spliced nucleotide sequences of predicted genes in fasta format

## Alignments for phylogenetic trees
apicoplast_concat_alignments.faa - alignment of apicoplast sequences

cytochrome_b_alignment.fa - cytochrome b alignment

11_genes_with_hepatocystis_epomophori_concat_alignments.fa - 11 gene, multi-organelle alignment

mitoch_proteins_concat_alignments.faa - alignment of mitochondrial sequences

nuclear_genome_proteins_concat_alignments.faa - alignment of nuclear genes

## Deconvolution of Hepatocystis bulk RNA-seq using Malaria Cell Atlas data and CIBERSORT
MCA_pseudobulk.R - R code describing how to generate stage-specific pseudobulk samples as a reference for bulk RNA-seq deconvolution

generate_mixtures.py - generate mixtures of pseudobulk life stages at known percentages to test CIBERSORT for accuracy

mca_pseudobulk_meroringschizfilt_cpm.dat - Output from running generate_mixtures.py, e.g. pseudobulk for MCA life stages


## Examination of missing genes in Hepatocystis relative to Plasmodium
heps_per_mca_cluster.py - find orthologue groups shared between P. berghei and either P. ovale or P. vivax, but absent from Hepatocystis and see whether Malaria Cell Atlas gene clusters are enriched for these.

hepatocystis_orthomcl.out - orthoMCL clusters for protein sequences from Hepatocystis DNA assembly (Hepatocystis_DNA), RNA assembly (Hepatocystis_RNA) and various Plasmodium species.

## Evolutionary analysis of genes
codeml_batch.py - Wrapper script for running codeml (http://envgen.nox.ac.uk/bioinformatics/docs/codeml.html) as batch.
Requires codeml to be installed and in path.

gatk_count_variants_per_sample.py - Script for finding the average number of variants per 10 kb of reference genome in a VCF file (derived using GATK by merging GVCF files of multiple samples)

gatk_count_variants_per_sample_sliding_window.py - Script for counting the number of variants per fixed length segments of the reference genome (default: 100 kb) in a VCF file using a sliding window.
Input: a VCF file derived using GATK by merging GVCF files of multiple samples. The script assumes that the assembly that the reads were mapped to for variant calling was concatenated into one pseudochromosome before mapping.
Output: a CSV file where the rows correspond to samples and the columns correspond to genome sequence bins.

