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

plot(CentPCA$projections,
     col=as.factor(meta2$region))
     legend("bottomright", legend=as.factor(unique(meta2$region)),
                                            fill=as.factor(unique(meta2$region)))
