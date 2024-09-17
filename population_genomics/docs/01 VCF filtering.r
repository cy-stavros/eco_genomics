library(vcfR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics")

list.files()

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")

vcf

head(vcf)

#bringing in reference genome
dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format = "fasta")


gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep = "\t", quote = "")

chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)

plot(chr1)

pdf(file = "~/projects/eco_genomics/population_genomics/figures/ChromoPlot.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))

dev.off()

#9/17/2024

DP <- extract.gt(vcf, element = "DP", as.numeric = T)

dim(DP)

#18233 SNPS, 629 somethings (individuals I think)

DP[1:5, 1:10]

quantile(DP)

DP[DP==0] <- NA

quantile(DP, na.rm = T)

# Visualize the matrix of DP (depth) and missingness in our VCF files

heatmap.bp(DP[1:1000,], rlabels=F, clabels=F)

#big scary graph ^

library(SNPfiltR)

vcf.filt <- hard_filter

#removes all genotypes with less than 3 reads
#5 or 10 would be more stringent

max_depth(vcf.filt)

#max depth rule of thumb: put ceiling at twice the average (so like 60 here)

vcf.filt <- max_depth(vcf.filt, maxdepth = 60)

#missingness time

meta <- read.csv("metadata/meta4vcf.csv", header=T)

meta2 <- meta[,c(1,4)] #subsetting to have all rows, columns 1-4 (format is rows, columns)

names(meta2) <- c("id", "pop") #b/c the function needs the name "pop"

meta2$id = as.factor(meta$id)
meta2$pop = as.factor(meta$pop)


vcf.filt.indMiss <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff = 0.75) #ex, 0.5 would filter out individuals missing 50% or more of its data

vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss) #so that we include only ones where there are two alleles! ig that kinda stuff is only relevant in things like phylogenetics
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)

vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff=0.5)

vcf.filt.indSNPMiss

DP2 <- extract.gt(vcf.filt.indSNPMiss,
                  element = "DP",
                  as.numeric = T)
heatmap.bp(DP2[1:5000,],
           rlabels=F,
           clabels=F)

write.vcf(vcf.filt.indSNPMiss,
          "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")
