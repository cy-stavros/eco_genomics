library(DESeq2)
library(ggplot2)
library(WGCNA)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(Rmisc)

options(bitmapType = "cairo") #helps resolve plottting issues when they crop up

#setwd("~/projects/eco_genomics2024/transcriptomics/") #idk why but it said "cannot change wd", 
#so I just clicked on the thingy in the session tab

#Step 1: import out data (counts matrix)

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = T, row.names = 1)
countsTableRound <- round(countsTable) #b/c DESeq2 doesn't like decimals
tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)
head(conds)

traitData <- read.table("/gpfs1/cl/pbio3990/Trait_Data.txt", header = T, row.names = 1)

#filter the matrix to just BASE data (cuz those are the data for which we have the traits measured)
filtered_counts_matrixBASEonly <- countsTable[, conds$FinalTemp == "BASE"]
filtered_sampled_metadata_BASEonly <- conds[conds$FinalTemp == "BASE", ]
rounded_filtered_count_matrix <- round(filtered_counts_matrixBASEonly)

#Step 2: Detecting outliers
#decet outlier genes
gsg <- goodSamplesGenes(t(rounded_filtered_count_matrix))
summary(gsg)

table(gsg$goodGenes)
#8203 good genes, 37235 bad genes???
#goodness and badness determined by stats things like overdispersion ig

# FALSE  TRUE 
# 37235 82203

table(gsg$goodSamples) #all good!

# filter out bad genes
data_WGCNA <- rounded_filtered_count_matrix[gsg$goodGenes==TRUE,]

dim(data_WGCNA)

#use clustering with a tree dendrogram to ID outlier samples
htree <- hclust(dist(t(data_WGCNA)), method="average")
plot(htree)
#looks like sN2C5 is the outgroup here, but according to "goodsamples" or whatever it's still good
#ig this will be a hw option

#PCA - outlier detection method
pca <- prcomp(t(data_WGCNA))
pca_data <- pca$x
#make a df
pca_data <- as.data.frame(pca_data)

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

ggplot(pca_data, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca_data)) +
  labs(x = paste0("PC1: ", pca.var.percent[1], "%"),
       y = paste0("PC2: ", pca.var.percent[2]))
#as you can see the outlier is still a big outlier here!

#Step 3 - Normalization

colData <- row.names(filtered_sampled_metadata_BASEonly)

dds_WGCNA <- DESeqDataSetFromMatrix(countData = data_WGCNA,
                                    colData = filtered_sampled_metadata_BASEonly,
                                    design =~1)
#design=~1 means that there are no specified groups!

dds_WGCNA_75 <- dds_WGCNA[rowSums(counts(dds_WGCNA) >= 15) >=6,]

#the sum of the counts for each of the genes has be at least 15, and it has to be present in at least 6 samples

nrow(dds_WGCNA_75)
#filtered down to 29,559 transcripts

dds_norm <- vst(dds_WGCNA_75)
#vst() is variance transformation
#perform variance stabilization

norm.counts <- assay(dds_norm) %>% 
  t()

#STEP 4: Network construction!!!!

#choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

#call the network topology analysis function (fakes a couple minutes to run)
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed", #meaning we're excluding negative corrs...
                         verbose = 5) 

#...positive correlations are more likely to be biologically relevant, negative corrs are more likely due to chance??? google this

sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x="Power", y="Scale free topology model fit, signed R^2") +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x="Power", y="Mean Connectivity") +
  theme_classic()

grid.arrange(a1, a2, nrow=2)

#basically the top graph shows how increasing power thresholds increase the statistical power of the model,
#but at the same time it decreases the avg number of connections each node has

