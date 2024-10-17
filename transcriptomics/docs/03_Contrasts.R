#load library
library(eulerr)

#set up groups within DESeq object
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group
dds <- DESeq(dds)
dim(dds)
resultsNames(dds)

# [1] "Intercept"               "group_D18A33_vs_D18A28"  "group_D18BASE_vs_D18A28" "group_D22A28_vs_D18A28" 
# [5] "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

#^^^^different levels to each sample or something
#don't totally get this

#by grouping them with this analysis you can compare between and within factors?
#asking first: does expression at acute temps differ between developmental temperature treatment groups

#1. compare baseline gene expression between dev treatment groups
res_D18_BASE_D22_BASE <- results(dds, contrast=c("group", "D18BASE", "D22BASE"), alpha = 0.05) #filtering for significance
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[!is.na(res_D18_BASE_D22_BASE$padj),] #filtering out NA pvals
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[order(res_D18_BASE_D22_BASE$padj),] #ordering based on pvals
head(res_D18_BASE_D22_BASE)
summary(res_D18_BASE_D22_BASE)

#make a list of which genes in our comparisons of interest are differnetially expressed (list of DEGs)
degs_D18_BASE_D22_BASE <- row.names(res_D18_BASE_D22_BASE[res_D18_BASE_D22_BASE$padj < .05,])
plotMA(res_D18_BASE_D22_BASE, ylim=c(-4,4))

#2. compare gene expression between developmental temperature treatment at A28
res_D18_A28_D22_A28 <- results(dds, contrast=c("group", "D18A28", "D22A28"), alpha = 0.05) #filtering for significance
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[!is.na(res_D18_A28_D22_A28$padj),] #filtering out NA pvals
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[order(res_D18_A28_D22_A28$padj),] #ordering A28d on pvals
head(res_D18_A28_D22_A28)
summary(res_D18_A28_D22_A28)

#make a list of which genes in our comparisons of interest are differnetially expressed (list of DEGs)
degs_D18_A28_D22_A28 <- row.names(res_D18_A28_D22_A28[res_D18_A28_D22_A28$padj < .05,])
plotMA(res_D18_A28_D22_A28, ylim=c(-4,4))


#3 compare gene expression between developmental temp treatment groups at A28
res_D18_A33_D22_A33 <- results(dds, contrast=c("group", "D18A33", "D22A33"), alpha = 0.05) #filtering for significance
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[!is.na(res_D18_A33_D22_A33$padj),] #filtering out NA pvals
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[order(res_D18_A33_D22_A33$padj),] #ordering A33d on pvals
head(res_D18_A33_D22_A33)
summary(res_D18_A33_D22_A33)

#make a list of which genes in our comparisons of interest are differnetially expressed (list of DEGs)
degs_D18_A33_D22_A33 <- row.names(res_D18_A33_D22_A33[res_D18_A33_D22_A33$padj < .05,])
plotMA(res_D18_A33_D22_A33, ylim=c(-4,4))

length(degs_D18_BASE_D22_BASE)
#1935 differentially expressed genes
length(degs_D18_A28_D22_A28)
#296
length(degs_D18_A33_D22_A33)
#78

#look at the overlaps in which genes are differentially expressed in multiple contrasts
length(intersect(degs_D18_BASE_D22_BASE, degs_D18_A28_D22_A28))
#107
length(intersect(degs_D18_BASE_D22_BASE, degs_D18_A33_D22_A33))
#44
length(intersect(degs_D18_A33_D22_A33, degs_D18_A28_D22_A28))
#29
length(intersect(degs_D18_BASE_D22_BASE, intersect(degs_D18_A33_D22_A33, degs_D18_A28_D22_A28)))
#23 #~nested interesects~ to compare all 3

#cool, see your venn diagram in your notes. picking up here tuesday.