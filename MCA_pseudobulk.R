# Generate pseudobulk dataset from MCA for deconvolution of bulk RNA-seq samples

# AUTHOR: Adam Reid
# Copyright (C) 2019 Genome Research Ltd.
# This program is distributed under the terms of the GNU General Public License

library(pheatmap)
library(scater)
library(M3Drop)
library(RColorBrewer)

setwd('/Users/ar11/Work/Writing/Papers/Malaria\ Cell\ Atlas/')

# Start with most recent scater object
mca.qc<-readRDS("data/mca.qc.tmm_ppt_20180907.rds")

# Make pseudobulk of summed single cells
cell_types = rownames(table(colData(mca.qc)$ShortenedLifeStage2))
cell_list <- rownames(colData(mca.qc)[colData(mca.qc)$ShortenedLifeStage2=="bbSpz",])
df <- data.frame("bbSpz" = rowSums(counts(mca.qc)[,cell_list]), row.names=rownames(counts(mca.qc)))
cell_list <- rownames(colData(mca.qc)[colData(mca.qc)$ShortenedLifeStage2=="EEF",])
df$EEF <- rowSums(counts(mca.qc)[,cell_list])
cell_list <- rownames(colData(mca.qc)[colData(mca.qc)$ShortenedLifeStage2=="Female",])
df$Female <- rowSums(counts(mca.qc)[,cell_list])
cell_list <- rownames(colData(mca.qc)[colData(mca.qc)$ShortenedLifeStage2=="Male",])
df$Male <- rowSums(counts(mca.qc)[,cell_list])

# For merozoites, exclude potential rings by taking values of ppt (pseudo-pseudotime) < 300
cell_list <- rownames(colData(mca.qc)[colData(mca.qc)$ShortenedLifeStage2=="Merozoite",])
cell_list_filt = cell_list[(colData(mca.qc)[cell_list,]$ppt < 300)]
df$Merozoite <- rowSums(counts(mca.qc)[,cell_list_filt])

cell_list <- rownames(colData(mca.qc)[colData(mca.qc)$ShortenedLifeStage2=="oocyst",])
df$oocyst <- rowSums(counts(mca.qc)[,cell_list])
cell_list <- rownames(colData(mca.qc)[colData(mca.qc)$ShortenedLifeStage2=="ook",])
df$ook <- rowSums(counts(mca.qc)[,cell_list])

# Subset rings to middle range of pseudotime to avoid confusion between what is a ring and what is a troph
# This might not actually be much of a problem, but certainly makes sense
cell_list <- rownames(colData(mca.qc)[colData(mca.qc)$ShortenedLifeStage2=="Ring",])
cell_list_filt = cell_list[(colData(mca.qc)[cell_list,]$ppt > 400) & (colData(mca.qc)[cell_list,]$ppt < 550)]
df$Ring <- rowSums(counts(mca.qc)[,cell_list_filt])

cell_list <- rownames(colData(mca.qc)[colData(mca.qc)$ShortenedLifeStage2=="ookoo",])
df$ookoo <- rowSums(counts(mca.qc)[,cell_list])

# Subset schizonts to middel range of pseudotime
cell_list <- rownames(colData(mca.qc)[colData(mca.qc)$ShortenedLifeStage2=="Schizont",])
cell_list_filt = cell_list[(colData(mca.qc)[cell_list,]$ppt > 750) & (colData(mca.qc)[cell_list,]$ppt < 900)]
df$Schizont <- rowSums(counts(mca.qc)[,cell_list_filt])

cell_list <- rownames(colData(mca.qc)[colData(mca.qc)$ShortenedLifeStage2=="sgSpz",])
df$sgSpz <- rowSums(counts(mca.qc)[,cell_list])

cell_list <- rownames(colData(mca.qc)[colData(mca.qc)$ShortenedLifeStage2=="Trophozoite",])
df$Trophozoite <- rowSums(counts(mca.qc)[,cell_list])

# Calculate CPMs
cpm = t(t(df*1000000)/colSums(df))

write.table(cpm, file="mca_pseudobulk_meroringschizfilt_cpm.dat", sep='\t', quote=FALSE)
