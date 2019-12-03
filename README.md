# hepatocystis-genome
Data and code relating to the Hepatocystis ex. Piliocolobus tephrosceles genome and transcriptome paper

## Primary genome assembly and annotation files
PRJEB32891_scaffolds.fasta.gz - assembly sequence in fasta format

PRJEB32891_union.embl.gz - assembly and annotation in single EMBL record (union)

PRJEB32891_scaffolds.embl.gz - assembly and annotation in EMBL format

PRJEB32891_proteins.faa - protein sequences of predicted genes in fasta format

PRJEB32891_transcripts.fa - spliced nucleotide sequences of predicted genes in fasta format


## Deconvolution of Hepatocystis bulk RNA-seq using Malaria Cell Atlas data and CIBERSORT
MCA_pseudobulk.R - R code describing how to generate stage-specific pseudobulk samples as a reference for bulk RNA-seq deconvolution

generate_mixtures.py - generate mixtures of pseudobulk life stages at known percentages to test CIBERSORT for accuracy

mca_pseudobulk_meroringschizfilt_cpm.dat - Output from running generate_mixtures.py, e.g. pseudobulk for MCA life stages


## Examination of missing genes in Hepatocystis relative to Plasmodium
heps_per_mca_cluster.py - find orthologue groups shared between P. berghei and either P. ovale or P. vivax, but absent from Hepatocystis and see whether Malaria Cell Atlas gene clusters are enriched for these.

hepatocystis_orthomcl.out - orthoMCL clusters for protein sequences from Hepatocystis DNA assembly (Hepatocystis_DNA), RNA assembly (Hepatocystis_RNA) and various Plasmodium species.

## Alignments for phylogenetic trees
apicoplast_concat_alignments.faa - alignment of apicoplast sequences

cytochrome_b_alignment.fa - cytochrome b alignment

11_genes_with_hepatocystis_epomophori_concat_alignments.fa - 11 gene, multi-organelle alignment

mitoch_proteins_concat_alignments.faa - alignment of mitochondrial sequences

nuclear_genome_proteins_concat_alignments.faa - alignment of nuclear genes
