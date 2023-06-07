# GM binary repeat of juan et al., data 

library(dplyr)
library(GenomicRanges)
library(tidyverse)
library(stringr)

metadata <- read.delim(file = "~/Data/GM12878/ENCODE/GM12878_ENCODE/metadata.tsv", header = T, sep = "\t")

metadata %>% 
  filter(Output.type == "peaks" & File.assembly == "hg19") %>% 
  select(File.accession, Experiment.target, Biological.replicate.s.)

ChIP_without_replicate <- metadata %>% 
  filter(Output.type == "peaks" & File.assembly == "hg19") %>% 
  filter(Biological.replicate.s. =="") %>% 
  dplyr::select(File.accession, Experiment.target)

Replicated <- metadata %>% 
  filter(Output.type == "peaks" & File.assembly == "hg19") %>% 
  filter(Biological.replicate.s. =="1, 2") %>% 
  dplyr::select(File.accession, Experiment.target)

ENCODE <- rbind(ChIP_without_replicate, Replicated)

# used something like this, not exactly because had to convert back to bed.gz but used a loop to convert all files that had a match in ENCODE.txt to the ChIP-seq name (e.g MYC instead of long ENCODE identifier) 

# for(i in 1:nrow(ENCODE)){
#   # Get the file name 
#   file_name <- ENCODE[i,1]
#   # Construct the full file path
#   file_path <- file.path("~/Data/GM12878/ENCODE/GM12878_ENCODE/CopyOfIn_use/", paste(file_name,"bed", sep = "."))
#   # Check if the file exists
#   if(file.exists(file_path)){
#     # Construct the new file path
#     new_file_path <- file.path("/home/dankent/Data/GM12878/ENCODE/In_use/", paste(ENCODE[i,2],"bed.gz", sep = "."))
#     }
# }

# Then this in command line 
# while IFS=$'\t' read -r oldname newname; do
# mv "$oldname" "$newname"
# done < ~/Data/MouseESC/Juan et al., data/GM12878_ENCODE/ENCODE.txt

# make this for creating binary files with ChromHMM 
ENCODE.txt <- read.delim(file = "~/Data/GM12878/ENCODE/In_use/ENCODE.txt", header = F, sep = "\t")
ENCODE.txt$V2 <- paste0(ENCODE.txt$V2, ".bed")
ENCODE.txt <- ENCODE.txt %>% mutate(GM = "GM12878", .before = V1) %>% 
  select(GM, V2)
ENCODE.txt <- ENCODE.txt %>% mutate(name = sub("-.*", "", ENCODE.txt$V2), .before = V2)
ENCODE.txt <- readr::write_tsv(ENCODE.txt, file = "/home/dankent/Data/GM12878/ENCODE/In_use/myDataFile_ENCODE.txt", col_names = F)

# run in command line with

# java -jar ChromHMM.jar BinarizeBed -peaks CHROMSIZES/hg19.txt ~/Data/GM12878/ENCODE/In_use ~/Data/GM12878/ENCODE/In_use/myDataFile_ENCODE.txt ~/Data/GM12878/ENCODE/In_use/binary/ 

# look
BHLHE40 <- read.delim(file = "~/Data/GM12878/ENCODE/In_use/BHLHE40-human.bed", header = F, sep = "\t")

# read in to R and add appropriate columns 
Chr1_binary <- read.delim(file = "~/Data/GM12878/ENCODE/In_use/binary/GM12878_chr1_binary.txt", header = F, sep = "\t")
Chr1_binary <- Chr1_binary[-1,]
colnames(Chr1_binary) <- Chr1_binary[1,]
Chr1_binary <- Chr1_binary[-1,]

# nrow(Chr1_binary)
# [1] 1246253

# read in genome size to create ranges
hg19_size <- read.delim("~/progs/ChromHMM-1.18/CHROMSIZES/hg19.txt", header = F)
hg19_size <- hg19_size[1:24,]

chr1 <- data.frame(start = seq(0, hg19_size[1,2], by = 200))
chr1 <- data.frame(seqnames = "chr1", start = (head(chr1,-1)), end = chr1[-1,])
chr1 <- makeGRangesFromDataFrame(chr1)

# length(chr1)
# [1] 1246253 

# SAME AS ABOVE SO EACH BIN MUST MATCH UP PERFECTLY! WOOOOOO

# now convert the above into a loop-able function

#split into chr for lapply function to create ranges
hg19_size_list <- split(hg19_size, f = hg19_size$V1)

# function to make ranges 
Create_chr_ranges <- function(list){
  x <- list[] # extract each member
 y <- data.frame(start = seq(0, x[,2], by = 200)) # seq along the entire length of the chrom and split into 200bp bins
 y <- data.frame(seqnames = x[,1], start = head(y, -1), end = y[-1,]) # make a data frame of the chr, start (- last start as this is the last row end), end (take away first end as this is 0 and the first start)
}

Chr_ranges <- lapply(hg19_size_list, Create_chr_ranges)

# read in all files from chromHMM

Chr10_binary <- read.delim(file = "~/Data/GM12878/ENCODE/In_use/binary/GM12878_chr10_binary.txt", header = F, sep = "\t")
Chr1_binary <- Chr1_binary[-1,]
colnames(Chr1_binary) <- Chr1_binary[1,]
Chr1_binary <- Chr1_binary[-1,]

# list files in directory 
list_files <- list.files("~/Data/GM12878/ENCODE/In_use/binary", pattern="*.txt", full.names=TRUE)
#read 'em all! and remove the weird names that ChromHMM adds to the first row
for (i in seq_along(1:length(list_files))) {
  
  # read the data
    data_import <- read.delim( paste0( "", list_files[i]), na = c( "NA", ""), header = F) # create object name based on filename
    data_import <- data_import[-1,]
    colnames(data_import) <- data_import[1,]
    data_import <- data_import[-1,]
    
      dataframe_name <- tolower( paste0("binary_", str_match(  # function to find a matching string
             list_files[i], "8_(.*?)_b"  # extract string between "8_" and "_b" in filename
           )[, 2]  # extract second element of the str_match() output
    )
  )
  
  # give the temporary data_import object the dataframe_name
    assign(dataframe_name, data_import)
    # output the filename and object name in the console
  print(list_files[i])
  print(dataframe_name)
}

head(binary_chr1) # wahooooo! It works 

## now add the chr, start and end from hg19_sizes to the binary files 

Chr1_combined <- cbind(Chr_ranges$chr1, binary_chr1)
Chr2_combined <- cbind(Chr_ranges$chr2, binary_chr2)
Chr3_combined <- cbind(Chr_ranges$chr3, binary_chr3)
Chr4_combined <- cbind(Chr_ranges$chr4, binary_chr4)
Chr5_combined <- cbind(Chr_ranges$chr5, binary_chr5)
Chr6_combined <- cbind(Chr_ranges$chr6, binary_chr6)
Chr7_combined <- cbind(Chr_ranges$chr7, binary_chr7)
Chr8_combined <- cbind(Chr_ranges$chr8, binary_chr8)
Chr9_combined <- cbind(Chr_ranges$chr9, binary_chr9)
Chr10_combined <- cbind(Chr_ranges$chr10, binary_chr10)
Chr11_combined <- cbind(Chr_ranges$chr11, binary_chr11)
Chr12_combined <- cbind(Chr_ranges$chr12, binary_chr12)
Chr13_combined <- cbind(Chr_ranges$chr13, binary_chr13)
Chr14_combined <- cbind(Chr_ranges$chr14, binary_chr14)
Chr15_combined <- cbind(Chr_ranges$chr15, binary_chr15)
Chr16_combined <- cbind(Chr_ranges$chr16, binary_chr16)
Chr17_combined <- cbind(Chr_ranges$chr17, binary_chr17)
Chr18_combined <- cbind(Chr_ranges$chr18, binary_chr18)
Chr19_combined <- cbind(Chr_ranges$chr19, binary_chr19)
Chr20_combined <- cbind(Chr_ranges$chr20, binary_chr20)
Chr21_combined <- cbind(Chr_ranges$chr21, binary_chr21)
Chr22_combined <- cbind(Chr_ranges$chr22, binary_chr22)

GM12878_ENCODE_features <- rbind(Chr1_combined, 
                                 Chr2_combined, 
                                 Chr3_combined, 
                                 Chr4_combined,
                                 Chr5_combined,                             
                                 Chr6_combined, 
                                 Chr7_combined, 
                                 Chr8_combined, 
                                 Chr9_combined,
                                 Chr10_combined, 
                                 Chr11_combined, 
                                 Chr12_combined,
                                 Chr13_combined,
                                 Chr14_combined,
                                 Chr15_combined,
                                 Chr16_combined,
                                 Chr17_combined,
                                 Chr18_combined,
                                 Chr19_combined,
                                 Chr20_combined,
                                 Chr21_combined,
                                 Chr22_combined)

rm(list=ls(pattern="*_combined"))
rm(list=ls(pattern="binary*"))

readr::write_tsv(x = GM12878_ENCODE_features, file = "~/Data/GM12878/ENCODE/GM12878_ENCODE_features.tsv", col_names = T)

