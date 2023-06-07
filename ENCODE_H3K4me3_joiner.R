library(dplyr)

GM12878_ENCODE_features <- read.delim("~/Data/GM12878/ENCODE/GM12878_ENCODE_features.tsv", header = T)
# remove every row with a 0 in every column, these are not informative and will help speed up the joiner

#long winded filter but couldn't get filter(across()) to work so filtered then subset oringinal by rownames
#copy
Original1 <- data.frame(GM12878_ENCODE_features)
Original1 <- data.frame(sapply(Original1[,4:62], as.numeric(Original1[,4:62])))
Original1 <- data.frame(GM12878_ENCODE_features[,1:3], Original1)
#add rownames
rownames(Original1) <- paste(Original1$seqnames, Original1$start, sep = "_")
#filter all features if EVERY column is 0
Original1 <- filter_all(Original1[4:62], any_vars(. != 0))
# check using row sums that nothing = 0
sum(rowSums(Original1) == 0 ) 
# rowSums(Original) useful because it shows that some rows only a 1, meaning that the code is working and only removing rows that have no 1s in. 

# need to merge back to get the seqnmaes, start and end 

H3K4me3_regions <- data.frame(GM12878_ENCODE_features)
#add rownames
rownames(H3K4me3_regions) <- paste(H3K4me3_regions$seqnames, H3K4me3_regions$start, sep = "_")

#use the rownames from the filtered dataframe to get the rest of the data back
H3K4me3_regions <- H3K4me3_regions[rownames(Original1),] 

l <- list()
i = 1
while (i < nrow(H3K4me3_regions)) {
  if (i %% 5000 == 0) print(paste0("We are at row ", i))
  iFirst = i
  while (H3K4me3_regions[i, "H3K4me3"] == 1 & H3K4me3_regions[i + 1, "H3K4me3"] == 1 & H3K4me3_regions[i, "end"] == H3K4me3_regions[i + 1, "start"] & H3K4me3_regions[i, "seqnames"] == H3K4me3_regions[i + 1, "seqnames"]) i = i + 1
  iLast = i
  # now we copy the last row we have checked
  newRow = H3K4me3_regions[i,]
  # we change the start to match the start in the beginning of the chain
  newRow[1, "start"] = H3K4me3_regions[iFirst,"start"]
  
  # we take the maximum for each of the rows that contain 0 or 1
  for (j in 4:62) newRow[1, j] = sum(H3K4me3_regions[iFirst:iLast, j])
  
  # adds this row to the end of our resulting dataframe
  l[[i]] <- newRow
  # increase the idex i to be ready for the next chain
  i = i + 1
}

Joined_H3K4me3_GM12878 <- data.table::rbindlist(l)
rm(l)

readr::write_tsv(x = Joined_H3K4me3_GM12878, file = "~/Data/GM12878/ENCODE/Joined_Genomes/Joined_H3K4me3_GM12878.tsv", col_names = T)
Rprof("test.out")
Rprof(NULL)
summaryRprof("test.out")

save.image()
