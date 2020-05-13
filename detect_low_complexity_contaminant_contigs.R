#!/usr/bin/env Rscript
# Script that was used for detecting contaminants in Hepatocystis genome assembly based on PCA and k-means clustering of a dataset that 
#  contains the GC%, tandem repeats percentage, frequency of stop codons and low complexity sequence content of each scaffold in the assembly.
# Input: path to a CSV file that contains the following:
# Column 1: scaffold names. Column 2: GC% of scaffolds. Column 3: tandem repeats percentage in scaffolds (based on the output of Tandem Repeats Finder).
# Column 4: frequency of stop codons. Column 5: low complexity sequence content as a percentage (based on the output of Dustmasker).
# Output: names of the scaffolds in the Hepatocystis assembly that are likely contaminants.
# Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.

library(ggfortify)
library(cluster)
library(fpc)

args = commandArgs(trailingOnly=TRUE)
df_in_path <- args[1]
df <- read.table(df_in_path, header=TRUE, sep=",", stringsAsFactors=FALSE, row.names=1)
df

pca_data <- prcomp(df, center=TRUE, scale.=TRUE)
plot(pca_data)
summary(pca_data)

autoplot(pca_data, loadings=TRUE, loadings.label=TRUE) + ggtitle("PCA of Hepatocystis assembly scaffolds") + theme(plot.title = element_text(hjust = 0.5))

scaled_df <- scale(df, center=TRUE, scale=TRUE)
df_cluster <- kmeans(scaled_df, centers=2)
plotcluster(scaled_df, df_cluster$cluster)

cluster1_names <- names(which(df_cluster$cluster == 1))
cluster2_names <- names(which(df_cluster$cluster == 2))

scaffs_to_remove <- NULL
if(length(cluster1_names) > length(cluster2_names)) {
  scaffs_to_remove <- cluster2_names
} else {
  scaffs_to_remove <- cluster1_names
}

print("scaffs_to_remove:")
for(scaff in scaffs_to_remove) {
  cat(paste0(scaff, "\n"))
}
