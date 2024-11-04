library(ggplot2)
library(DESeq2)

options(bitmapType = "cairo")

#setwd("~/projects/eco_genomics2024/transcriptomics/")

#importing counts matrix

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = T, row.names = 1)
# rounding b/c DESeq2 doesn't like decimals
countsTableRound <- round(countsTable) 
tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)
head(conds)

#seeing how many reads we have from each sample
colSums(countsTableRound) #millions of counts for each sample!
mean(colSums(countsTableRound)) #the average count number is 18,454,529!
#20 million reads per sample is the gold standard of what you shoot for.

#note that the transcriptome was not made for this study, so many genes are gonna have 0 counts
mean(rowSums(countsTableRound)) #3244.739
median(rowSums(countsTableRound)) #64! a few genes that have a lot of expression, a lot that have a lil 
#(skewed right)

apply(countsTableRound,2,mean)
#gives a sense of variation in sequencing effort across samples

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

#skipping pca stuff
#skipping all of script two (i think)


######################## Contrasts ###############################
#set up groups within DESeq object
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group
dds <- DESeq(dds)
dim(dds)
resultsNames(dds)

# [1] "Intercept"               "group_D18A33_vs_D18A28"  *"group_D18BASE_vs_D18A28"* "group_D22A28_vs_D18A28" 
# [5] "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

