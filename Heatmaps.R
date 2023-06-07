# Plots 

library(pheatmap)
library(ggplot2)
library(dplyr)
library(seriation)
library(dendextend)

# H3K4me3 ####

Normalised_H3K4me3_Region <- read.delim(file = "~/Data/GM12878/ENCODE/Joined_Genomes/Normalised_H3K4me3_ENCODE.tsv", header = T)

Large_heatmap <- Normalised_H3K4me3_Region %>% filter(bin_number > 2)
Large_heatmap <- Large_heatmap[,4:62]

jpeg(file='~/Data/GM12878/ENCODE/Plots/ENCODE_H3K4me3_Large_heatmap_RowRegion.jpeg', width = 1500, height = 800)
par(mfrow = c(1, 1))
pheatmap::pheatmap(Large_heatmap, fontsize_row = 9, fontsize_col = 1, cluster_rows = T, cluster_cols = F, cutree_rows = 5, transpose = T)
dev.off()

# H3K27ac ##library(seriation)
library(dendextend)##

Normalised_H3K27ac_Region <- read.delim(file = "~/Data/GM12878/ENCODE/Joined_Genomes/Normalised_H3K27ac_ENCODE.tsv", header = T)

Large_heatmap <- Normalised_H3K27ac_ENCODE %>% filter(bin_number > 2)
Large_heatmap <- Large_heatmap[,4:62]

jpeg(file='~/Data/GM12878/ENCODE/Plots/ENCODE_H3K27ac_Large_heatmap_RowRegion.jpeg', width = 1500, height = 800)
par(mfrow = c(1, 1))
pheatmap::pheatmap(Large_heatmap, fontsize_row = 9, fontsize_col = 1, cluster_rows = T, cluster_cols = F, cutree_rows = 5, transpose = T)
dev.off()

# Hybrids ####

Normalised_Hybrids_Region <- read.delim(file = "~/Data/GM12878/ENCODE/Joined_Genomes/Normalised_Hybrids_Region.tsv", header = T)

Large_heatmap <- Normalised_Hybrids_Region %>% filter(bin_number > 2)
Large_heatmap <- Large_heatmap[,4:62]

jpeg(file='~/Data/GM12878/ENCODE/Plots/ENCODE_Hybrid_Large_heatmap_RowRegion.jpeg', width = 1500, height = 800)
par(mfrow = c(1, 1))
pheatmap::pheatmap(Large_heatmap, fontsize_row = 9, fontsize_col = 1, cluster_rows = T, cluster_cols = F, cutree_rows = 5, transpose = T)
dev.off()

rm(Large_heatmap)

phtmap <- pheatmap(Large_heatmap)
# rotating (or more accurately, not transposing) show that EP300 is clustering with the big boys 

