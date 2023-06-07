library(dplyr)

# H3K4me3 ####
Normalised_H3K4me3_Region <- read.delim("~/Data/GM12878/ENCODE/Joined_Genomes/Joined_H3K4me3_GM12878.tsv", sep = "\t", header = T) %>%  
  filter(H3K4me3 > 0)

Normalised_H3K4me3_Region$width <- Normalised_H3K4me3_Region$end - Normalised_H3K4me3_Region$start

Normalised_H3K4me3_Region$bin_number <- round(Normalised_H3K4me3_Region$width/200)

tmp <- Normalised_H3K4me3_Region[,4:62] # only numeric

# Normalised_H3K4me3_Region$bin_number <- temp$bin_number

tmp <- tmp[,1:59]/Normalised_H3K4me3_Region$bin_number  #normalise all bar bin_number and width 
Normalised_H3K4me3_Region <- cbind(Normalised_H3K4me3_Region[,1:3],tmp, Normalised_H3K4me3_Region[,c(63,64)])
Normalised_H3K4me3_Region %>% 
  filter(H3K4me3 >1 | H3K4me3 <1 ) # no values for H3K4me3 are below 1 or above 1. Meaning that each regions has been normalised to the bins of the origin mark, in this case H3K4me3

readr::write_tsv(x = Normalised_H3K4me3_Region, file = "~/Data/GM12878/ENCODE/Joined_Genomes/Normalised_H3K4me3_ENCODE.tsv", col_names = T)

# H3K27ac ####

Normalised_H3K27ac_Region <- read.delim("~/Data/GM12878/ENCODE/Joined_Genomes/Joined_H3K27ac_GM12878.tsv", sep = "\t", header = T) %>%  
  filter(H3K27ac > 0)

Normalised_H3K27ac_Region$width <- Normalised_H3K27ac_Region$end - Normalised_H3K27ac_Region$start

Normalised_H3K27ac_Region$bin_number <- round(Normalised_H3K27ac_Region$width/200)

tmp <- Normalised_H3K27ac_Region[,4:62] # only numeric

# Normalised_H3K27ac_Region$bin_number <- temp$bin_number

tmp <- tmp[,1:59]/Normalised_H3K27ac_Region$bin_number  #normalise all bar bin_number and width 

Normalised_H3K27ac_Region <- cbind(Normalised_H3K27ac_Region[,1:3],tmp, Normalised_H3K27ac_Region[,c(63,64)])

Normalised_H3K27ac_Region %>% 
  filter(H3K27ac >1 | H3K27ac <1 ) # no values for H3K27ac are below 1 or above 1. Meaning that each regions has been normalised to the bins of the origin mark, in this case H3K27ac 

readr::write_tsv(x = Normalised_H3K27ac_Region, file = "~/Data/GM12878/ENCODE/Joined_Genomes/Normalised_H3K27ac_ENCODE.tsv", col_names = T)

# Hybrids ####

Normalised_Hybrids_Region <- read.delim("~/Data/GM12878/ENCODE/Joined_Genomes/Joined_Hybrids_ENCODE.tsv", sep = "\t", header = T) %>%  
  filter(H3K27ac > 0 | H3K4me3 > 0)

Normalised_Hybrids_Region$width <- Normalised_Hybrids_Region$end - Normalised_Hybrids_Region$start

Normalised_Hybrids_Region$bin_number <- round(Normalised_Hybrids_Region$width/200)

tmp <- Normalised_Hybrids_Region[,6:64] # only numeric

tmp <- tmp[,1:59]/Normalised_Hybrids_Region$bin_number  #normalise all bar bin_number and width 

Normalised_Hybrids_Region <- data.frame(Normalised_Hybrids_Region[,1:3],tmp, Normalised_Hybrids_Region[,65])

Normalised_Hybrids_Region <- Normalised_Hybrids_Region %>% rename(bin_number = Normalised_Hybrids_Region...65.)

readr::write_tsv(x = Normalised_Hybrids_Region, file = "~/Data/GM12878/ENCODE/Joined_Genomes/Normalised_Hybrids_Region.tsv", col_names = T)

