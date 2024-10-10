library(ggplot2)
library(DESeq2)

options(bitmapType = "cairo")

setwd("~/projects/eco_genomics2024/transcriptomics/") #idk why but it said "cannot change wd", 
#so I just clicked on the thingy in the session tab

#import counts matrix

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = T, row.names = 1)
countsTableRound <- round(countsTable) #b/c DESeq2 doesn't like decimals
tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)
head(conds)

########################## Explore counts matrix ##############################

#seeing how many reads we g=have from each sample
colSums(countsTableRound) #millions of counts for each sample!
mean(colSums(countsTableRound)) #the average count number is 18,454,529!
#20 million reads per sample is the gold standard of what you shoot for.

barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound), 
        cex.names = 0.5, las = 2, ylim = c(0,30000000))
abline(h=mean(colSums(countsTableRound)), col = "blue4", lwd=2)

#the average number of counts per gene
rowSums(countsTableRound)
#note that the transcriptome was not made for this study, so many genes are gonna have 0 counts
mean(rowSums(countsTableRound)) #3244.739
median(rowSums(countsTableRound)) #64! a few genes that have a lot of expression, a lot that have a lil 
#(skewed right? tail is on right)

apply(countsTableRound,2,mean) #not super sure what this is the mean of
#251.something is the highest group here
#"gives a sense of variation in sequencing effort across samples"
#appears to be somehow related to the transcript count barplot

##################### starting analysis in DESeq2 ###########################


dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~DevTemp + FinalTemp)
dim(dds)

dds <- dds[rowSums(counts(dds) >= 10) >= 15,] 
#filtering, only looking at transcripts with at least 10 reads, things less than that might not give good info
#transcripts must be present in at least 15 samples

nrow(dds)
#35527 transcripts less (not just genes) (way more transcripts than genes b/c slice variation)

#run the DESeq model to test for global differential gene expression
dds <- DESeq(dds)
#different read numbers per sample (look back at the barplot!) are normalized by ~size factors~

#listing the results you've generated
resultsNames(dds)

#"Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28"  "FinalTemp_BASE_vs_A28"

# visualizing our global gene expression patterns using PCA!
#First we need to transform the data for plotting using variance stabilization

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("DevTemp", "FinalTemp"), returnData = TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))

final_temp_colors <- c("BASE" = "grey", "A28" = "hotpink", "A33" = "red")
shapes_choose <- c("D18" = 16, "D22" = 18)

p <- ggplot(pcaData, aes(PC1, PC2, color = FinalTemp, shape = DevTemp)) +
  geom_point(size = 5) +
  scale_shape_manual(values = shapes_choose) +
  scale_color_manual(values = final_temp_colors) +
  labs(x = paste0('PC1: ', percentVar[1], ' %'),
       y = paste0('PC2: ', percentVar[2], ' %')) +
  theme_bw(base_size = 16)

#cool! looks like most clustering is happening between finaltemp, not dev temp which is kinda interesting
    