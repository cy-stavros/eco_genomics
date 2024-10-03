library(vcfR)
library(SNPfiltR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")


#not sure if i this is necessary, but bringing in reference genome

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format = "fasta")

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep = "\t", quote = "")

chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)

plot(chr1)

#no snp densities found?

chromoqc(chr1, xlim=c(1e1, 1.1e8))

DP <- extract.gt(vcf, element = "DP", as.numeric = T) #DP is depth!

dim(DP)

#18233 SNPS, 629 individuals

quantile(DP) #not sure why this only has 330 individuals??

DP[DP==0] <- NA

quantile(DP, na.rm = T) #still only 350 some

heatmap.bp(DP[1:1000,], rlabels=F, clabels=F)

vcf.filt <- hard_filter(vcf, depth=3)
#removes all genotypes with less than 3 reads
#5 or 10 would be more stringent

max_depth(vcf.filt)

#put ceiling at 2x the average (so like 60 here)

vcf.filt <- max_depth(vcf.filt, maxdepth = 60)

#filtering missingness!

meta <- read.csv("metadata/meta4vcf.csv", header=T)

meta2 <- meta[,c(1,4)] #subsetting to have all rows, columns 1-4 (format is rows, columns)

names(meta2) <- c("id", "pop") #b/c the function needs the name "pop"

meta2$id = as.factor(meta$id)
meta2$pop = as.factor(meta$pop)

#creating missingness variable i can come back and change later

#indMisscuttoff = 0.75 actually nvm i'll just create 3 separate vcfs

vcf.filt.indMiss.75 <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff = 0.75)

vcf.filt.indMiss.80 <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff = 0.80)

vcf.filt.indMiss.85 <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff = 0.85)
#alrighty saving my stuff now to return to later
#also i just realized i filtered in the wrong direction lol