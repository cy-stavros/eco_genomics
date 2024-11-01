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
#82203 good genes, 37235 bad genes???
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
# we want to maximize mean connectivity (avg # of connections while also maximizing scale free topology model fit)
#selecting threshold of 26

soft_power <- 26
temp_cor <- cor
cor <- WGCNA::cor #this sets the temp_cor function to use WGCNA's correlation function

norm.counts[] <- sapply(norm.counts, as.numeric)

#the comman below creates the network and identifies modules based on the parameters that we chose
bwnet26 <- blockwiseModules(norm.counts,
                            maxBlockSize = 30000,
                            TOMType = "signed", #only focusing on positive correlations
                            power = soft_power,
                            mergeCutHeight = .25,
                            numericLabels = FALSE,
                            randomSeed=1234,
                            verbose = 3)

#takes a hot minute ^^^^, in hw request more memory

#run again! L in labels wasn't capitalized and an extra s in random seed :/


cor <- temp_cor #this resets cor function to base R's cor function of using WCGNA's

#STEP 5: Explore Module Eigengenes

module_eigengenes <- bwnet26$MEs

head(module_eigengenes)
dim(module_eigengenes)
#lets goooo! 51 modules for .2.6. [incase i acccidentally find and replace]
#get the number of genes for each module
table(bwnet26$colors)


#plot the dendrogram and the module colors:
plotDendroAndColors(bwnet26$dendrograms[[1]], cbind(bwnet26$unmergedColors, bwnet26$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

saveRDS(bwnet26, file = "outputs/bwnet26.rds")

#to load the bwnet file in later

#(make sure you're in the transcriptomics directory when you do this)

#step 6: Correlation of modules with traits!

#define the numbers of genes & samples
nSamples <- nrow(norm.counts) #7
nGenes <- ncol(norm.counts) #29559 ? (sample script has 1750)
#i think there is a serious issue here bc there are far more genes that colors! (slateblue.1 problem!)
#on second thought, it's not like there's 1750 colors either. 
#still weird it's different than the sample script, but there must be some other reason my colors aren't showing up

# test for a correlation between module eigengenes and trait data
module.trait.corr <- cor(module_eigengenes, traitData, use = 'p') #'p' is pearson's correlation

#calculating pvalies for each correlation

module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap 
heatmap.data <- merge(module_eigengenes, traitData, by = 'row.names')
head(heatmap.data)

#address error of row.names not being numeric
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

names(heatmap.data)

#make pretty heatmap of correlations

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[52:54], #these values might need to change based on the # of eigengenes
             y = names(heatmap.data)[1:41],
             col = c("blue2", "skyblue", "white", "pink", "red"))

#cool! make sure to slap some labels on here.
#trying to figure out what the heck the three final columns actually mean

#also note to self: save some other way besides zooming & right click b/c they look corrupted or something on github
