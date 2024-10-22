#load in day 1 transcriptomics script and run DESeq object
library(pheatmap)
resultsNames(dds)
options(bitmapType = "cairo")

# "Intercept" "DevTemp_D22_vs_D18" "FinalTemp_A33_vs_A28" "FinalTemp_BASE_vs_A28"

#pull out the results for Developmental Temperature 22 vs 18
res_D22vsD18 <- results(dds, name = "DevTemp_D22_vs_D18", alpha = .05)

#order by significance
res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),] #ordering by pval (adjusted (i'm assuming))

head(res_D22vsD18)
#the top of this list has the transcripts that have the most differentiated expression between 22 and 18

#look at counts of a specific top gene that we're interested in to validate that the model is working
#checking out top one
d <- plotCounts(dds, gene ="TRINITY_DN140854_c0_g5_i2", int = (c("DevTemp", "FinalTemp")), returnData = T)

p <- ggplot(d, aes(x=DevTemp, y=count, color=DevTemp, shape=FinalTemp)) +
  theme_minimal() + theme(text=element_text(size=20), panel.grid.major=element_line(color="grey"))

#points not appearing, asking it to add them lol
p <- p + geom_point(position=position_jitter(w=0.2, h=0), size=3)

p

#checking this graph is good to check magnitude, as well as directionality!
#so logfold2 -1.19 is a little over 2 times more in D18.

#MA plot (not manhattan plot!)
#M is another way to put log fold change and A is average

plotMA(res_D22vsD18, ylim=c(-4,4))
#ones on the far right are highly expressed, but relatively few (dna synth, metabolism, etc)
#another way visualize overdispersion!

#volcano plot!

#convert our DESeq results object into a df to plot:
res_df <- as.data.frame(res_D22vsD18)
#add a column to indicate significance:
res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")
#puts another constraint on significance by making magnitude be at least a doubling (log fold change > 1)

#plotting

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=Significant))+
  geom_point(alpha=0.8) +
  scale_color_manual(values = c("slateblue", "tomato")) +
  labs(x = "Log2 Fold Change", y="log10 Adjusted P-value", title = "Volcano Plot") +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange") +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color = "orange")

#cool! more upregulation in 22/downreggulation in 18 (upper right quatdrant thingy)
#remember: each point is a transcript here


#now looking at expression thru a heatmap

vsd <- vst(dds, blind=FALSE)

topgenes <- head(rownames(res_D22vsD18), 20)
mat <- assay(vsd)[topgenes,]
df <- as.data.frame(colData(dds)[,c("DevTemp", "FinalTemp")])
pheatmap(mat, annotation_col=df, show_rownames=FALSE, clusrer_cols=T, cliter_rows=T)
#columns are samples, rows are genes
#trees: phenogram, hierearcical clustering thingy that's trying to group both scripts and samples on a tree
#made by ubmga? or something, doesn't actually imply phylogeny
#transcripts clustering together could be on the same pathway possibly
#notice individuals are clustered by devtemp, not so much by finaltemp!

