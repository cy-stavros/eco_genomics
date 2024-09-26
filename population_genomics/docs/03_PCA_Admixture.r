library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

options(bitmapType = "cairo")

setwd("~/projects/eco_genomics/population_genomics/")

vcf <- read.vcfR("outputs/vcf_final.filtered.vcf.gz")

#we need to thin the SNPs for LD (linkage diequalibirum) before we run PCA and 
# Admixture analyses to satisfy the assumptuions of independence among loci

vcf.thin <- distance_thin(vcf, min.distance = 500) #Making new vcf that eliminates SNPS between 500 bp

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)

meta2 <- meta[meta$id %in% colnames(vcf.thin@gt[, -1]) , ]

dim(meta2)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")

#hide uncompressed vcf file (too big for github!) outside of our repo

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf") 
#system() basically lets you do stuff like you would in command line

geno <- vcf2geno(input.file="/gpfs1/home/c/s/cstavros/vcf_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned.geno")

CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale=TRUE)

#if you've already done PCA preiviously, you can load the results without running it again like so:
CentPCA <- load.pcaProject("vcf_final.filtered.thinned.pcaProject")

show(CentPCA)

plot(CentPCA) #shows scree plot that shows eigen values for each principal component
#units aren't important, but you can tell it levels off very quickly

#plot(CentPCA$projections,
#     col=as.factor(meta2$region))
#     legend("bottomright", legend=as.factor(unique(meta2$region)),
#                                            fill=as.factor(unique(meta2$region)))

#making a prettier version in ggplot

ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent)) +
       geom_point(alpha=.5) +
  labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region",shape="Continent")
#(if you wanna zoom in) + xlim(-10,10) + ylim(-10,10)

#clear signs of structure!!!
# 2 clouds for NE corresponding to C. jaceaea and C. nigra!

#saving!
#also note: the percentage of variation each eigenvector explains is its eigenvalue / sum(all eigenvalues)
ggsave("figures/CentPCA_PC1vPC2.pdf", width=6, height=6, units="in")

