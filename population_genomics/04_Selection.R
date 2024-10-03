library(tidyverse)
library(ggplot2)
library(vcfR)
library(qqman)
library(pcadapt)

vcf <- read.pcadapt("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf",
                    type="vcf")
#vcfR <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz") #this one had a different number of rows, so taking from class set

vcfR <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf.gz")

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta2 <- meta[meta$id %in% colnames(vcfR@gt[,-1]),]

pcadapt.pca <- pcadapt(vcf,
                       K=2,
                       method="componentwise",
                       min.maf=0.01,
                       LD.clumping = list(size=500, thr=0.2)) #thinning out for LD
summary(pcadapt.pca)
plot(pcadapt.pca, options="stat.distribution", pop=meta2$region,
     i=1, j=2,
     K=2)
#results in basically a manhattan plot without being divided up into chromosomes

View(head(vcfR@fix))

vcfR.fix <- as.data.frame(vcfR@fix[,1:2])

chr.main <- unique(vcfR.fix$CHROM)[1:8] #filtering for just the chromosomes

chrnum <- as.data.frame(cbind(chr.main, seq(1,8,1)))

Pval <- pcadapt.pca$pvalues #NAs are alleles w/ freq <1% (we filtered those out!)

pcadapt.MHplot <- cbind(vcfR.fix, Pval)
pcadapt.MHplot <- left_join(chrnum, pcadapt.MHplot, join_by(chr.main==CHROM)) #stuck here for a while

pcadapt.MHplot <- pcadapt.MHplot %>% 
  mutate(SNP=paste0(chr.main, "_", POS))

str(pcadapt.MHplot)

pcadapt.MHplot$V2 = as.numeric(pcadapt.MHplot$V2)
pcadapt.MHplot$POS = as.numeric(pcadapt.MHplot$POS)
pcadapt.MHplot$pPC1 = as.numeric(pcadapt.MHplot[,4])
pcadapt.MHplot$pPC2 = as.numeric(pcadapt.MHplot[,5])

pcadapt.MHplot <- pcadapt.MHplot %>%  drop_na(pPC1)

view(pcadapt.MHplot)

manhattan(pcadapt.MHplot,
          chr="V2",
          bp="POS",
          p="pPC1",
          col=c("blue4", "orange3"),
          logp=T,
          ylab=("-log10 p-value"),
          genomewideline = F,
          main = "PCAdapt genome scan for selection (PC1)")

manhattan(pcadapt.MHplot,
          chr="V2",
          bp="POS",
          p="pPC2",
          col=c("blue4", "orange3"),
          logp=T,
          ylab=("-log10 p-value"),
          genomewideline = F,
          main = "PCAdapt genome scan for selection (PC2)")

View(pcadapt.MHplot %>% 
       filter(pPC1<quantile(pcadapt.MHplot$pPC1, 0.001)) %>% 
select(chr.main, POS, pPC1))
#viewing the most extreme outliers

# a peak of several outliers can be evidence that a gene (not just a snp! is under selection)

hist(pcadapt.pca$maf)
