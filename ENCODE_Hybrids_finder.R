# Hybrid maker 

library(plyranges)
library(GenomicRanges)
library(dplyr)

# replicate start end ends for checks and creating new ranges for hybrids :) 

Joined_H3K4me3_GM12878$start_H3K4me3 <- Joined_H3K4me3_GM12878$start
  
Joined_H3K4me3_GM12878$end_H3K4me3 <- Joined_H3K4me3_GM12878$end  
  
Joined_H3K27ac_GM12878$start_H3K27ac <- Joined_H3K27ac_GM12878$start
  
Joined_H3K27ac_GM12878$end_H3K27ac <- Joined_H3K27ac_GM12878$end

Hybrids_final <- data.frame( unique( join_overlap_inner( makeGRangesFromDataFrame( Joined_H3K4me3_GM12878), makeGRangesFromDataFrame( Joined_H3K27ac_GM12878, keep.extra.columns = TRUE)))) # general overlap WITH k4 FIRST 

starts    <- cbind(data.frame(Hybrids_final$start),data.frame(Hybrids_final$start_H3K27ac)) #extract all the  starts
ends      <- cbind(data.frame(Hybrids_final$end),data.frame(Hybrids_final$end_H3K27ac)) #extract all the ends

head(data.frame(start = apply(starts,1,which.min), end = apply(ends,1,which.max)), 7) 
# start end
# 1     1   1
# 2     2   1
# 3     1   1
# 4     2   1
# 5     1   1
# 6     2   1
# 7     2   1

starts    <- apply(starts,1,min)
ends      <- apply(ends,1,max) 
chr       <- data.frame(Hybrids_final$seqnames)
Hybrids_final     <- data.frame(seqnames = chr, start = starts, end = ends)
Hybrids_final     <- makeGRangesFromDataFrame(Hybrids_final, seqnames.field = "Hybrids_final.seqnames")

# overlap with the original to get all the feature info 

GM12878_ENCODE_features <- read.delim( "~/Data/GM12878/ENCODE/GM12878_ENCODE_features.tsv", header = T)

# make granges to overlap with the original 200bp bins 

Hybrids_ENCODE <- makeGRangesFromDataFrame( GM12878_ENCODE_features, keep.extra.columns = T)

Hybrids_ENCODE <- data.frame( unique( join_overlap_inner( Hybrids_ENCODE, Hybrids_final)))

Hybrids_ENCODE <- Hybrids_ENCODE %>% filter( H3K4me3 != 0 | H3K27ac != 0)

l <- list()
i = 1
while (i <= nrow(Hybrids_ENCODE)) {
  if (i %% 5000 == 0) print(paste0("We are at row ", i))
  iFirst = i
  while (((i < nrow(Hybrids_ENCODE) - 1) & Hybrids_ENCODE[i, "H3K4me3"] > 0 | Hybrids_ENCODE[i, "H3K27ac"] > 0) & Hybrids_ENCODE[i, "end"] == Hybrids_ENCODE[i + 1, "start"] & Hybrids_ENCODE[i, "seqnames"] == Hybrids_ENCODE[i + 1, "seqnames"]) i = i + 1
  iLast = i
  # now we copy the last row we have checked
  newRow = Hybrids_ENCODE[i,]
  # we change the start to match the start in the beginning of the chain
  newRow[1, "start"] = Hybrids_ENCODE[iFirst,"start"]
  # we take the maximum for each of the rows that contain 0 or 1
  for (j in 6:64) newRow[1, j] = sum(Hybrids_ENCODE[iFirst:iLast, j]) # now 6:64 because GRanges has added width and strand
  # adds this row to the end of our resulting list
  l[[length(l) + 1]] <- newRow
  # increase the idex i to be ready for the next chain
  i = i + 1
}

Hybrids_ENCODE <- data.frame(data.table::rbindlist(l))
Hybrids_ENCODE$width <- Hybrids_ENCODE$end - Hybrids_ENCODE$start
Hybrids_ENCODE$bin_number <- Hybrids_ENCODE$width/200
rm(l)
rm(Hybrids_final)

readr::write_tsv(x = Hybrids_ENCODE, file = "~/Data/GM12878/ENCODE/Joined_Genomes/Joined_Hybrids_ENCODE.tsv", col_names = T)
