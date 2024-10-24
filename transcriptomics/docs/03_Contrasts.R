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

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_A28_D22_A28 <- row.names(res_D18_A28_D22_A28[res_D18_A28_D22_A28$padj < .05,])
plotMA(res_D18_A28_D22_A28, ylim=c(-4,4))


#3 compare gene expression between developmental temp treatment groups at A33
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

###### making euler plots cont'd ######

# calculate the # of unique genes in each portion of the Euler plot
1935-107-44+23 #(adding back in the triple overlap so we don't double-count)
#1807 transcripts diff expressed uniquely at BASE btwn 18&22

296-107-29+23 #183 uniquely expressed when exposed to 28

78-44-29+23 #28 genes uniquely expressed whhen exposed to 33

107-23 #84 unique to BASE & 28

44-23 #21 unique to BASE & A33

29-23 #6 unique to A28 & A33

myEuler <- euler(c("BASE"=1807, "A28" = 183, "A33" = 28, "BASE&A28" = 84, 
                   "BASE&A33" = 21, "A28&A33" = 6, "BASE&A28&A33" = 23))

plot(myEuler, lty=1:3, quantities=TRUE)

#cool!

###### make a scatterplot of responses to A28/23 when copepods develop at 18 vs 22 #####

# contrast D18_A28 vs BASE
res_D18_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group","D18BASE","D18A28"), alpha = 0.05))

#contrast D22_A28 vs BASE
res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group","D22BASE","D22A28"), alpha = 0.05))

#merge dfs
res_df28 <- merge(res_D18_BASEvsA28, res_D22_BASEvsA28, by = "row.names", suffixes = c(".18",".22"))

#putting row.names back as true rownames and then deleting
rownames(res_df28) <- res_df28$Row.names
res_df28 <- res_df28[,-1]

library(dplyr)
library(tidyr)

#define colormapping logic with the mutate function

res_df28 <- res_df28 %>%  
  mutate(fill = case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "turquoise2",
    padj.18 < 0.05 & stat.18 > 0 ~ "magenta",
    padj.22 < 0.05 & stat.22 < 0 ~ "blue2",
    padj.22 < 0.05 & stat.22 > 0 ~ "red"
  ))

#plot:

ggplot(res_df28, aes(x = log2FoldChange.18, y=log2FoldChange.22, color = fill)) +
  geom_point(alpha=0.8) +
  scale_color_identity() +
  labs(x= "Log2FoldChange 28 vs. BASE at 18",
       y= "Log2FoldChange 28 vs. BASE 22",
       title= "How does response to 28 C vary by DevTemp?") +
  theme_minimal()

#neat!

#count the # of points per fill color
color_counts <- res_df28 %>% 
  group_by(fill) %>% 
  summarise(count=n())

label_positions <- data.frame(
  fill= c("blue", "magenta", "red", "turquoise2"),
  x_pos = c(1,5,0,-7.5),
  y_pos = c(-5, 0, 9, 3)
)

label_data <- merge(color_counts, label_positions, by = "fill")

plot28 <- ggplot(res_df28, aes(x = log2FoldChange.18, y = log2FoldChange.22, color = fill)) +
       geom_point(alpha=0.8) +
         scale_color_identity() +
  geom_text(data = label_data, aes(x=x_pos, y= y_pos, label = count, color = fill), size=5)+
  geom_abline(intercept = 0, slope=1, linetype = "dashed", color ="black")+
  geom_abline(intercept = 0, slope=-1, linetype= "dashed", color = "grey")+
  xlim(-10,10) + ylim(-10,10)+
  labs(x="Log2FoldChange 28 vs BASE at 18",
       y="Log2FoldCHange 28 vs BASE at 22",
       title = "How does response to 28 C vary by DevTemp?") +
  theme_minimal()

#what is label_data?
# + the label stuff look at your picture of sage's script

#copied & pasted from above + then will modify for 33 (maybe on hw???)


######### Repeating everything for A33 ########

# contrast D18_A33 vs BASE
res_D18_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group","D18BASE","D18A33"), alpha = 0.05))

#contrast D22_A33 vs BASE
res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group","D22BASE","D22A33"), alpha = 0.05))

#merge dfs
res_df33 <- merge(res_D18_BASEvsA33, res_D22_BASEvsA33, by = "row.names", suffixes = c(".18",".22"))

#putting row.names back as true rownames and then deleting
rownames(res_df33) <- res_df33$Row.names
res_df33 <- res_df33[,-1]

library(dplyr)
library(tidyr)

#define colormapping logic with the mutate function

res_df33 <- res_df33 %>%  
  mutate(fill = case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "turquoise2",
    padj.18 < 0.05 & stat.18 > 0 ~ "magenta",
    padj.22 < 0.05 & stat.22 < 0 ~ "blue2",
    padj.22 < 0.05 & stat.22 > 0 ~ "red"
  ))

#plot:

ggplot(res_df33, aes(x = log2FoldChange.18, y=log2FoldChange.22, color = fill)) +
  geom_point(alpha=0.8) +
  scale_color_identity() +
  labs(x= "Log2FoldChange 33 vs. BASE at 18",
       y= "Log2FoldChange 33 vs. BASE 22",
       title= "How does response to 33 C vary by DevTemp?") +
  theme_minimal()




label_data <- merge(color_counts, label_positions, by = "fill")

plot33 <- ggplot(res_df33, aes(x = log2FoldChange.18, y = log2FoldChange.22, color = fill)) +
  geom_point(alpha=0.8) +
  scale_color_identity() +
  geom_text(data = label_data, aes(x=x_pos, y= y_pos, label = count, color = fill), 
            size=5)+
  geom_abline(intercept = 0, slope=1, linetype = "dashed", color ="black")+
  geom_abline(intercept = 0, slope=-1, linetype= "dashed", color = "grey")+
  xlim(-10,10) + ylim(-10,10)+
labs(x="Log2FoldChange 33 vs BASE at 18",
     y="Log2FoldCHange 33 vs BASE at 22",
     title = "How does response to 33 C vary by DevTemp?") +
  theme_minimal()

#idk why the lower count label isn't working
#also the count labels are the same for 33 for some reason :/
library(gridExtra)

combined_plot <- grid.arrange(plot28, plot33, ncol = 2)

ggsave("~/projects/eco_genomics/transcriptomics/figures/combined_scatter_plot.png", 
       combined_plot, width = 12, height = 6)
#ope nvm figured it out

