# hepatocystis-genome
Data and code relating to the Hepatocystis ex. Piliocolobus tephrosceles genome and transcriptome paper

## Deconvolution of Hepatocystis bulk RNA-seq using Malaria Cell Atlas data and CIBERSORT
MCA_pseudobulk.R - R code describing how to generate stage-specific pseudobulk samples as a reference for bulk RNA-seq deconvolution

generate_mixtures.py - generate mixtures of pseudobulk life stages at known percentages to test CIBERSORT for accuracy

data


## Examination of missing genes in Hepatocystis relative to Plasmodium

heps_per_mca_cluster.py - find orthologue groups shared between P. berghei and either P. ovale or P. vivax, but absent from Hepatocystis and see whether Malaria Cell Atlas gene clusters are enriched for these.

OrthoMCL data
