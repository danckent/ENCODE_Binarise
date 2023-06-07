# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggridges)
library(tidyverse)
library(viridis)
library(gridExtra)
library(factoextra)
############### 
##### RAW ###
##############

Hybrids_ENCODE <- read.delim("Joined_Genomes/Joined_Hybrids_ENCODE.tsv") 

# Read in data and select relevant columns
data <- Hybrids_ENCODE %>%
  select(H3K4me3, H3K27ac, width, bin_number)

# Group data by region
data_grouped <- data %>%
  group_by(width)

# Count number of regions in each bin
data_grouped <- group_by(data_grouped, width) 

#plots 
par(mfrow = c(1, 4))

#histograms 
data %>% 
  select(H3K4me3, H3K27ac, everything()) %>%  # reorder the columns using "select"
  gather() %>%  # "gather" the data into a long format
  ggplot(aes(x = value)) +  # specify the x-axis value
  geom_histogram(bins = 20) +  # create a histogram with 20 bins
  facet_wrap(~ key, scales = "free_x")  # create separate plots for each column

# scatter H3K4me3 vs H3K27ac raw
data %>% 
  ggplot(aes(x = H3K4me3, y = H3K27ac)) +
  geom_point() +  # create a scatter plot
  geom_smooth(method = "lm") + # add a linear regression line
  ylim(0,160) +
  xlim(0,160)

ggplot(data, aes(x = H3K4me3, y = H3K27ac)) +
  geom_hex(aes(alpha = ..density..)) +  # compute the density of the points and map it to the alpha aesthetic
  # geom_point(color = "black") +  # add a scatter plot
  geom_smooth(method = "lm") # add a linear regression line

################
## Normalised ##
################

Normalised_Hybrids_Region <- read.delim("Joined_Genomes/Normalised_Hybrids_Region.tsv") 

# Read in data and select relevant columns
data <- Normalised_Hybrids_Region %>%
  select(H3K4me3, H3K27ac, bin_number)

# Group data by region
data_grouped <- data %>%
  group_by(bin_number)

# Count number of regions in each bin
data_grouped <- group_by(data_grouped, bin_number) 

jpeg(file= "~/Data/GM12878/ENCODE/Plots/Normalised_Hybrids_Region_histagrams.jpeg", width = 1500, height = 800)

par(mfrow = c(1, 2))
# histograms ####
data %>% 
  select(H3K4me3, H3K27ac, everything()) %>%  # reorder the columns using "select"
  gather() %>%  # "gather" the data into a long format
  ggplot(aes(x = value)) +  # specify the x-axis value
  geom_histogram(bins = 20) +  # create a histogram with 20 bins
  facet_wrap(~ key, scales = "free_x")  # create separate plots for each column

dev.off()

# scatter H3K4me3 vs H3K27ac raw ####
data %>% 
  ggplot(aes(x = H3K4me3, y = H3K27ac)) +
  geom_point() +  # create a scatter plot
  geom_smooth(method = "lm")  # add a linear regression line

# tile ####
ggplot(data, aes(x = H3K4me3, y = H3K27ac)) +
  stat_bin2d(aes(fill = ..count..)) +  # count the number of points in each tile and map it to the fill aesthetic
  scale_fill_gradient(low = "white", high = "blue")  # set the color gradient for the fill aesthetic

# ridgeline ####

points <- data.frame(data_grouped %>%
                       count(bin_number) %>%
                       filter(n < 2) %>%
                       select(bin_number) %>%
                       inner_join(data_grouped, by = "bin_number"))

data_grouped$bin_number <- factor(data_grouped$bin_number, levels = c(1:160))

plot1 <- ggplot(data_grouped, aes(x = H3K4me3, y = as.factor(bin_number), fill = stat(x))) +
  geom_density_ridges_gradient(aes(y = as.factor(bin_number)), scale = 3, rel_min_height = 0.05) +
  scale_fill_viridis_c(name = "Enrichment", option = "C") +
  geom_point(data = points, aes(x = H3K4me3, y = as.factor(bin_number)), show.legend = F) +
  labs(title = 'H3K4me3') +
  ylab("Region Size (200bp bins)")

plot2 <- ggplot(data_grouped, aes(x = H3K27ac, y = as.factor(bin_number), fill = stat(x))) +
  geom_density_ridges_gradient(aes(y = as.factor(bin_number)), scale = 3, rel_min_height = 0.05) +
  scale_fill_viridis_c(name = "Enrichment", option = "C") +
  geom_point(data = points, aes(x = H3K27ac, y = as.factor(bin_number)), show.legend = F) +
  labs(title = 'H3K27ac') +
  ylab("Region Size (200bp bins)")

jpeg(file="~/Data/GM12878/ENCODE/Plots/Hybrid_Ridgeplot.jpeg", width = 1500, height = 800)
par(mfrow = c(1, 2))
grid.arrange(plot1, plot2, nrow = 1)
dev.off()

### this just makes a plot Normalised values, averaged and plotted ####
delete <- data.frame(Mean_enrichment = colMeans(x = Large_heatmap))
delete$Feature <- rownames(delete)

delete <- delete %>% arrange(desc(Mean_enrichment)) 

levels(delete$Feature) <- factor(delete$Feature)
levels(delete$Mean_enrichment) <- factor(delete$Mean_enrichment)

barplot(delete[order(delete[,1],decreasing=T),
][,1],names.arg=delete[order(delete[,1],decreasing=T),]
[,2], las = 2)

##########


########## copied from Juan data, not edited yet ####
pca <- Normalised_Hybrids

pca <- pca %>% mutate(label = case_when(bin_number > 0 & bin_number <= 5 ~ "0-1kb",
                                        bin_number > 5 & bin_number <= 10 ~ "1-2kb",
                                        bin_number > 10 & bin_number <= 20 ~ "2-4kb",
                                        bin_number > 20 ~ ">4kb"))

pca1 <- pca[,1:82] # bin_number is the grouping
# Perform PCA
pca1 <- prcomp(pca1, scale = FALSE)

#
jpeg(file='Normalised Hybrids Coloured by Size Range.jpeg', width = 1500, height = 800)
par(mfrow = c(1, 1))
fviz_pca_ind(pca1, geom = "point", col.ind = pca$label, )
dev.off()

eig.val <- get_eigenvalue(pca1)
fviz_eig(pca1, addlabels = TRUE, ylim = c(0, 50), linecolor = "orange")


# remove regions below 5 bins 


pca <- Normalised_Hybrids %>% filter(bin_number >5)
pca <- pca %>% mutate(label = case_when(
  bin_number > 5 & bin_number <= 10 ~ "1-2kb",
  bin_number > 10 & bin_number <= 20 ~ "2-4kb",
  bin_number > 20 ~ ">4kb"))

pca1 <- pca[,1:82] # bin_number is the grouping
# Perform PCA
pca1 <- prcomp(pca1, scale = FALSE)

#
jpeg(file='Normalised Hybrids Coloured by Size Range.jpeg', width = 1300, height = 800)
par(mfrow = c(1, 1))
fviz_pca_ind(pca1, geom = "point", col.ind = pca$label, )
dev.off()

eig.val <- get_eigenvalue(pca1)
fviz_eig(pca1, addlabels = TRUE, ylim = c(0, 50), linecolor = "orange")


df <- data.frame(pca1$rotation[, 1])

df %>% 
  ggplot(aes(x = rownames(df), y = pca1.rotation...1.)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ylab("PCA Component 1") + 
  xlab("")


