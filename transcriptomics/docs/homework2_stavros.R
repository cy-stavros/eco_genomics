library(ggplot2)
library(DESeq2)
library(dplyr)
library(tidyr)
library(eulerr)

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


######################## Contrasts - D18 ###############################
#set up groups within DESeq object
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group #unsure whether or not i need to modify this
dds <- DESeq(dds)
dim(dds)
resultsNames(dds)

# [1] "Intercept"               "group_D18A33_vs_D18A28"  *"group_D18BASE_vs_D18A28"* "group_D22A28_vs_D18A28" 
# [5] "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

#1. compare gene expression between D18BASE and D18A28
#(reversing order so that upregulation means higher expression in 28!!!)
res_D18_BASE_D18_A28 <- results(dds, contrast=c("group", "D18A28", "D18BASE"), alpha = 0.05) #filtering for significance
res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[!is.na(res_D18_BASE_D18_A28$padj),] #filtering out NA pvals
res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[order(res_D18_BASE_D18_A28$padj),] #ordering based on pvals
head(res_D18_BASE_D18_A28)
summary(res_D18_BASE_D18_A28)

# out of 35524 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 30, 0.084%
# LFC < 0 (down)     : 11, 0.031%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

#note to self here: come back and be sure you know which group up and down reg is in reference to.
#okay for future reference the format is contrasts = c("condition", "treated", "untreated")

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_BASE_D18_A28 <- row.names(res_D18_BASE_D18_A28[res_D18_BASE_D18_A28$padj < .05,])
plotMA(res_D18_BASE_D18_A28, ylim=c(-4,4))

#2. compare gene expression between D18BASE and D18A33
#switching the factors to make a proper comparison again
res_D18_BASE_D18_A33 <- results(dds, contrast=c("group", "D18A33", "D18BASE"), alpha = 0.05) #filtering for significance
res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[!is.na(res_D18_BASE_D18_A33$padj),] #filtering out NA pvals
res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[order(res_D18_BASE_D18_A33$padj),] #ordering based on pvals
head(res_D18_BASE_D18_A33)
summary(res_D18_BASE_D18_A33)

# out of 35524 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 240, 0.68%
# LFC < 0 (down)     : 92, 0.26%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_BASE_D18_A33 <- row.names(res_D18_BASE_D18_A33[res_D18_BASE_D18_A33$padj < .05,])
plotMA(res_D18_BASE_D18_A33, ylim=c(-4,4))

length(degs_D18_BASE_D18_A28)
#41 differentially expressed genes
length(degs_D18_BASE_D18_A33)
#332 DEGs

#look at the overlaps in which genes are differentially expressed in both contrasts
length(intersect(degs_D18_BASE_D18_A28, degs_D18_BASE_D18_A33))
#34 genes differentially expressed in both contrasts

#making euler plots

41-34 #7 transcripts unique to BASE vs 28
332-34 #298 transcripts unique to BASE vs 33

myEulerD18 <- euler(c("A28" = 7, "A33" = 298, "A28&A33" = 34))

EuD18 <- plot(myEulerD18, lty=1:3, quantities=TRUE)

#man, euler plots don't look as cool with just two circles (unless i wanna add in d22?)
#not sure if that would be appropriate tho b/c its comparing to a different baseline?

###### make a scatterplot of responses to A28/33 when copepods develop at 18 #####
#going back in and rerunning w/ factors in proper order now

# contrast D18_A28 vs D18_A33
res_D18_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group","D18A28", "D18BASE"), alpha = 0.05))

#contrast D22_A28 vs BASE
res_D18_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group", "D18A33", "D18BASE"), alpha = 0.05))

#merge dfs
res_dfD18 <- merge(res_D18_BASEvsA28, res_D18_BASEvsA33, by = "row.names", suffixes = c(".28",".33"))

#putting row.names back as true rownames and then deleting
rownames(res_dfD18) <- res_dfD18$Row.names
res_dfD18 <- res_dfD18[,-1]

#define colormapping logic with the mutate function

res_dfD18 <- res_dfD18 %>%  
  mutate(fill = case_when(
    padj.28 < 0.05 & stat.28 < 0 ~ "turquoise2",
    padj.28 < 0.05 & stat.28 > 0 ~ "magenta",
    padj.33 < 0.05 & stat.33 < 0 ~ "blue2",
    padj.33 < 0.05 & stat.33 > 0 ~ "red"
  ))

ggplot(res_dfD18, aes(x = log2FoldChange.28, y=log2FoldChange.33, color = fill)) +
  geom_point(alpha=0.8) +
  scale_color_identity() +
  labs(x= "Log2FoldChange 28 vs. BASE",
       y= "Log2FoldChange 33 vs. BASE",
       title= "How does response vary by acute temp?") +
  theme_minimal()

#cool!

color_counts <- res_dfD18 %>% 
  group_by(fill) %>% 
  summarise(count=n())

label_positions <- data.frame(
  fill= c("blue", "magenta", "red", "turquoise2"),
  x_pos = c(1,5,0,-7.5),
  y_pos = c(-5, 0, 9, 3)
)

label_data <- merge(color_counts, label_positions, by = "fill")

plotD18 <- ggplot(res_dfD18, aes(x = log2FoldChange.28, y=log2FoldChange.33, color = fill)) +
  geom_point(alpha=0.8) +
  scale_color_identity() +
  geom_text(data = label_data, aes(x=x_pos, y= y_pos, label = count, color = fill), size=5)+
  geom_abline(intercept = 0, slope=1, linetype = "dashed", color ="black")+
  geom_abline(intercept = 0, slope=-1, linetype= "dashed", color = "grey")+
  xlim(-10,10) + ylim(-10,10)+
  labs(x= "Log2FoldChange 28 vs. BASE",
       y= "Log2FoldChange 33 vs. BASE",
       title= "How does response vary by acute temp? (D18)") +
  theme_minimal()

plotD18

#cool, but like before, the lower label is just not showing up. also no scripts were posted from this day :/
#looks like this is will be a job for ppt text boxes lol

#interp: upper left quadrant: upregged in 33, downregged in 28
# upper right: upregged in both 33 and 28
# lower right: downregged in 33, upregged in 28
# lower left: downregged in 33, downregged in 33

#be sure about the polarity of this! it might be opposite.
#okay after a lotta digging it looks like the opposite is true. 
#I'm hitting the hay for now, but in the morning I'll have to 
#reverse the order on all the contrast statments.
#fixed!!!

#all of my work on 22 didn't save :(


######################## Contrasts - D22 ###############################
#set up groups within DESeq object
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group #unsure whether or not i need to modify this
dds <- DESeq(dds)
dim(dds)
resultsNames(dds)

# [1] "Intercept"               "group_D18A33_vs_D18A28"  *"group_D18BASE_vs_D18A28"* "group_D22A28_vs_D18A28" 
# [5] "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

#1. compare gene expression between D22BASE and D22A28
#(reversing order so that upregulation means higher expression in 28!!!)
res_D22_BASE_D22_A28 <- results(dds, contrast=c("group", "D22A28", "D22BASE"), alpha = 0.05) #filtering for significance
res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[!is.na(res_D22_BASE_D22_A28$padj),] #filtering out NA pvals
res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[order(res_D22_BASE_D22_A28$padj),] #ordering based on pvals
head(res_D22_BASE_D22_A28)
summary(res_D22_BASE_D22_A28)

# out of 30703 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 15, 0.049%
# LFC < 0 (down)     : 274, 0.89%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

#note to self here: come back and be sure you know which group up and down reg is in reference to.
#okay for future reference the format is contrasts = c("condition", "treated", "untreated")

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D22_BASE_D22_A28 <- row.names(res_D22_BASE_D22_A28[res_D22_BASE_D22_A28$padj < .05,])
plotMA(res_D22_BASE_D22_A28, ylim=c(-4,4))

#2. compare gene expression between D22BASE and D22A33
#switching the factors to make a proper comparison again
res_D22_BASE_D22_A33 <- results(dds, contrast=c("group", "D22A33", "D22BASE"), alpha = 0.05) #filtering for significance
res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[!is.na(res_D22_BASE_D22_A33$padj),] #filtering out NA pvals
res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[order(res_D22_BASE_D22_A33$padj),] #ordering based on pvals
head(res_D22_BASE_D22_A33)
summary(res_D22_BASE_D22_A33)

# out of 32768 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 388, 1.2%
# LFC < 0 (down)     : 1176, 3.6%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D22_BASE_D22_A33 <- row.names(res_D22_BASE_D22_A33[res_D22_BASE_D22_A33$padj < .05,])
plotMA(res_D22_BASE_D22_A33, ylim=c(-4,4))

length(degs_D22_BASE_D22_A28)
#289 differentially expressed genes
length(degs_D22_BASE_D22_A33)
#1564 DEGs

#look at the overlaps in which genes are differentially expressed in both contrasts
length(intersect(degs_D22_BASE_D22_A28, degs_D22_BASE_D22_A33))
#144 genes differentially expressed in both contrasts

#making euler plots

289-144 #145 transcripts unique to BASE vs 28
1564-144 #1420 transcripts unique to BASE vs 33

myEulerD22 <- euler(c("A28" = 145, "A33" = 1420, "A28&A33" = 144))

EuD22 <- plot(myEulerD22, lty=1:3, quantities=TRUE)

#man, euler plots don't look as cool with just two circles (unless i wanna add in d22?)
#not sure if that would be appropriate tho b/c its comparing to a different baseline?

###### make a scatterplot of responses to A28/33 when copepods develop at 22 #####
#going back in and rerunning w/ factors in proper order now

# contrast D22_A28 vs D22_A33
res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group","D22A28", "D22BASE"), alpha = 0.05))

#contrast D22_A28 vs BASE
res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group", "D22A33", "D22BASE"), alpha = 0.05))

#merge dfs
res_dfD22 <- merge(res_D22_BASEvsA28, res_D22_BASEvsA33, by = "row.names", suffixes = c(".28",".33"))

#putting row.names back as true rownames and then deleting
rownames(res_dfD22) <- res_dfD22$Row.names
res_dfD22 <- res_dfD22[,-1]

#define colormapping logic with the mutate function

res_dfD22 <- res_dfD22 %>%  
  mutate(fill = case_when(
    padj.28 < 0.05 & stat.28 < 0 ~ "turquoise2",
    padj.28 < 0.05 & stat.28 > 0 ~ "magenta",
    padj.33 < 0.05 & stat.33 < 0 ~ "blue2",
    padj.33 < 0.05 & stat.33 > 0 ~ "red"
  ))

ggplot(res_dfD22, aes(x = log2FoldChange.28, y=log2FoldChange.33, color = fill)) +
  geom_point(alpha=0.8) +
  scale_color_identity() +
  labs(x= "Log2FoldChange 28 vs. BASE",
       y= "Log2FoldChange 33 vs. BASE",
       title= "How does response vary by acute temp?") +
  theme_minimal()

#cool!

color_counts <- res_dfD22 %>% 
  group_by(fill) %>% 
  summarise(count=n())

label_positions <- data.frame(
  fill= c("blue", "magenta", "red", "turquoise2"),
  x_pos = c(1,5,0,-7.5),
  y_pos = c(-5, 0, 9, 3)
)

label_data <- merge(color_counts, label_positions, by = "fill")

plotD22 <- ggplot(res_dfD22, aes(x = log2FoldChange.28, y=log2FoldChange.33, color = fill)) +
  geom_point(alpha=0.8) +
  scale_color_identity() +
  geom_text(data = label_data, aes(x=x_pos, y= y_pos, label = count, color = fill), size=5)+
  geom_abline(intercept = 0, slope=1, linetype = "dashed", color ="black")+
  geom_abline(intercept = 0, slope=-1, linetype= "dashed", color = "grey")+
  xlim(-10,10) + ylim(-10,10)+
  labs(x= "Log2FoldChange 28 vs. BASE",
       y= "Log2FoldChange 33 vs. BASE",
       title= "How does response vary by acute temp? (D22)") +
  theme_minimal()

plotD22

#cool, but like before, the lower label is just not showing up. also no scripts were posted from this day :/
#looks like this is will be a job for ppt text boxes lol

#interp: upper left quadrant: upregged in 33, downregged in 28
# upper right: upregged in both 33 and 28
# lower right: downregged in 33, upregged in 28
# lower left: downregged in 33, downregged in 33

####### making combined graphs:########

library(gridExtra)

combined_plot <- grid.arrange(plotD18, plotD22, ncol = 2)

ggsave("~/projects/eco_genomics/transcriptomics/figures/HWcombined_scatter_plot.png", 
       combined_plot, width = 12, height = 6)



png(filename = "~/projects/eco_genomics/transcriptomics/figures/HWcombined_euler.png")

combined_euler <- grid.arrange(EuD18, EuD22, ncol = 2)

dev.off()
