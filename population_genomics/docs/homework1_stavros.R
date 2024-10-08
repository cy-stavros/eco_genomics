library(vcfR)
library(SNPfiltR)
library(tidyverse)
library(qqman)
library(LEA)

#helps solve plotting issues
X11.options(type="cairo")
options(bitmapType = "cairo")

setwd("/gpfs1/cl/pbio3990/PopulationGenomics")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")


#not sure if this is necessary, but bringing in reference genome 
#try running w/o it later

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

#least stringent - 36 inds removed
vcf.filt.indMiss.75 <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff = 0.75)
#more stringent - 58 inds removed
vcf.filt.indMiss.70 <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff = 0.70)
#most stringent - 84 inds removed
vcf.filt.indMiss.65 <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff = 0.65)

#alrighty saving my stuff now to return to later

#picking back up monday
#first, filtering the .75

vcf.filt.indMiss.75 <- filter_biallelic(vcf.filt.indMiss.75) #removing non-biallelic
vcf.filt.indMiss.75 <- min_mac(vcf.filt.indMiss.75, min.mac = 1) #mac = Minor Allele Count, min.mac=1 excludes all minor alleles

vcf.filt.indSNPMiss.75 <- missing_by_snp(vcf.filt.indMiss.75, cutoff=0.5) #filtering missingness by SNP

vcf.filt.indSNPMiss.75 #taking a peek at it:
# ***** Object of Class vcfR *****
# 593 samples
# 19 CHROMs
# 15,456 variants
# Object size: 103.4 Mb
# 29.84 percent missing data
# *****        *****         *****

#creating a file of depth
DP2.75 <- extract.gt(vcf.filt.indSNPMiss.75,
                  element = "DP",
                  as.numeric = T)
#viewing it on the heatmap -note, not necessary, we can prob leave this out for the other two/delete later
heatmap.bp(DP2.75[1:5000,],
           rlabels=F,
           clabels=F)

#writing this so I don't have to rerun this all in the future
write.vcf(vcf.filt.indSNPMiss.75,
          "~/projects/eco_genomics/population_genomics/outputs/vcf.75.filtered.vcf.gz")

#alrighty, filtering done for this level of missingness. moving on to the analyses. 
#we'll reiterate all that's above for the next two levels later

#we need to find:

#b. Genome-wide diversity for each region (Hs) and its standard deviation (SD)
#c. The number of loci with Hs=0 vs. Hs>0 for each region 
#d. Genetic structure using either PCA or Admixture analysis

#finding genome-wide diversity for each region (Hs) and its Std Dev:

#importing metadata again for next steps

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta) # our vcf file has 593 samples, but
dim(meta) # meta has 629 inds

newmeta2 <- meta[meta$id %in% colnames(vcf.filt.indSNPMiss.75@gt[,-1]),] #filtering meta to consist only of the inds remaining post-filtering
dim(newmeta2) #now it matches!

#calculate diversity stats using genetic_diff function in vcfR
vcf.filt.indSNPMiss.75.div <- genetic_diff(vcf.filt.indSNPMiss.75,
                        pops=as.factor(newmeta2$region),
                        method = "nei")

str(vcf.filt.indSNPMiss.75.div) #peeking at structure

chr.main <- unique(vcf.filt.indSNPMiss.75.div$CHROM)[1:8] #taking the first 8 chromosomes (ignoring the tiny scaffolds)

chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1))) #numbering the chromosomes

#don't need to MH plot, so skipping the plotting steps but keeping w/ the same filtering
#in case i forget, .Hs is gonna have the same filtering steps as .MHplot in 02

#joining chromosome num to our vcf
vcf.filt.indSNPMiss.75.div.Hs <- left_join(chrnum, vcf.filt.indSNPMiss.75.div, join_by(chr.main==CHROM)) 

#filtering out Gst(Fst) values <0 and creating "SNP" column w/ both chromo # and bp position
vcf.filt.indSNPMiss.75.div.Hs <- vcf.filt.indSNPMiss.75.div.Hs %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main,"_",POS))

#making these columns numeric, not characters
vcf.filt.indSNPMiss.75.div.Hs$V2 = as.numeric(vcf.filt.indSNPMiss.75.div.Hs$V2)
vcf.filt.indSNPMiss.75.div.Hs$POS = as.numeric(vcf.filt.indSNPMiss.75.div.Hs$POS)


#note to self: all this chromosome stuff might be completely unnecessary. consider removing when cleaning up the script
#this is where 02 moves on to Hs
names(vcf.filt.indSNPMiss.75.div.Hs) #Hs values are in columns 4-9

#taking a peek at what our Hs values look like for each region w/ a histogram
vcf.filt.indSNPMiss.75.div.Hs %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) +
  geom_histogram(position="identity", alpha=0.5, bins=50) +
  labs(title="Genome-wide expected heterozygosity (Hs)",fill="Regions",
       x="Gene diversity (Hs) within Regions", y="Counts of SNPs")

#generating table with the summary statistics we need!!!

loci <- nrow(vcf.filt.indSNPMiss.75.div.Hs) #15403 loci at this stage
#to get number of loci where Hs = 0, we'll do loci - the N we calculate after filtering values that equal 0

Hs_table.75 <- vcf.filt.indSNPMiss.75.div.Hs %>%
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  filter(value!=0) %>% 
  summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n(), N_0Hs=loci-n()) 

# A tibble: 6 x 5
# name   avg_Hs StdDev_Hs  N_Hs N_0Hs
# <chr>   <dbl>     <dbl> <int> <int>
#1 Hs_CEU  0.174     0.158  8132  7271
#2 Hs_NE   0.105     0.139 14041  1362
#3 Hs_NEU  0.198     0.155  7297  8106
#4 Hs_PNW  0.147     0.158  9947  5456
#5 Hs_SEU  0.242     0.153  5041 10362
#6 Hs_WEU  0.192     0.155  7671  7732

write.csv(Hs_table.75, "~/projects/eco_genomics/population_genomics/outputs/Hs_table.75_noZeros.csv",
          quote=F,
          row.names=F)

#alrighty we got our Hs and StdDev by region, as well as the number of Hs=0 vs. Hs>0 for each region!
#taking a break here. saving all my stuff!

#alright i'm back. i now see the value of saving this stuff lol. I don't wanna run this all again!
#saved the vcf.gz this time!
#now to get cracking on the pca
#d. Genetic structure using either PCA or Admixture analysis

#incase i wanna pick back up where i was
#vcf.filt.indSNPMiss.75 <- read.vcfR("/outputs/vcf.75.filtered.vcf.gz")

#thinning out, eliminating SNPs within 500bp
vcf.filt.indSNPMiss.75.thin <- distance_thin(vcf.filt.indSNPMiss.75, min.distance = 500)

#rereading metadata
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)

#including only the filtered individuals
pcameta2 <- meta[meta$id %in% colnames(vcf.filt.indSNPMiss.75@gt[, -1]) , ]

dim(pcameta2)

#sneaking under the door (cuz github doesn't like big files)
setwd("~/projects/eco_genomics/population_genomics/")
write.vcf(vcf.filt.indSNPMiss.75.thin, "outputs/vcf.75.filtered.thinned.vcf.gz")

#stuck here

#hide uncompressed vcf file (too big for github!) outside of our repo

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf.75.filtered.thinned.vcf.gz > ~/vcf.75.filtered.thinned.vcf") 
#system() basically lets you do stuff like you would in command line

geno.75 <- vcf2geno(input.file="/gpfs1/home/c/s/cstavros/vcf.75.filtered.thinned.vcf",
                 output.file = "outputs/vcf.75.filtered.thinned.geno")

PCA.75 <- LEA::pca("outputs/vcf.75.filtered.thinned.geno", scale=TRUE)

ggplot(as.data.frame(PCA.75$projections),
       aes(x=V1, y=V2, color=pcameta2$region, shape=pcameta2$continent)) +
  geom_point(alpha=.5) +
  labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region",shape="Continent")

#yay! saving rn
ggsave("figures/PCA.75_PC1vPC2.pdf", width=6, height=6, units="in")





#####moving on to .70!#####




vcf.filt.indMiss.70 <- filter_biallelic(vcf.filt.indMiss.70) #removing non-biallelic
vcf.filt.indMiss.70 <- min_mac(vcf.filt.indMiss.70, min.mac = 1) #mac = Minor Allele Count, min.mac=1 excludes all minor alleles

vcf.filt.indSNPMiss.70 <- missing_by_snp(vcf.filt.indMiss.70, cutoff=0.5) #filtering missingness by SNP

vcf.filt.indSNPMiss.70 #taking a peek at it:

#***** Object of Class vcfR *****
#  571 samples
#19 CHROMs
#15,802 variants
#Object size: 102.5 Mb
#28.74 percent missing data
#*****        *****         *****


#creating a file of depth
DP2.70 <- extract.gt(vcf.filt.indSNPMiss.70,
                     element = "DP",
                     as.numeric = T)
write.vcf(vcf.filt.indSNPMiss.70,
          "~/projects/eco_genomics/population_genomics/outputs/vcf.70.filtered.vcf.gz")

#finding genome-wide diversity for each region (Hs) and its Std Dev:
  
  #importing metadata again for next steps
  
  meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta) # our vcf file has 571 samples, but
dim(meta) # meta has 629 inds

newmeta2 <- meta[meta$id %in% colnames(vcf.filt.indSNPMiss.70@gt[,-1]),] #filtering meta to consist only of the inds remaining post-filtering
dim(newmeta2) #now it matches!

#calculate diversity stats using genetic_diff function in vcfR
vcf.filt.indSNPMiss.70.div <- genetic_diff(vcf.filt.indSNPMiss.70,
                                           pops=as.factor(newmeta2$region),
                                           method = "nei")

str(vcf.filt.indSNPMiss.70.div) #peeking at structure

chr.main <- unique(vcf.filt.indSNPMiss.70.div$CHROM)[1:8] #taking the first 8 chromosomes (ignoring the tiny scaffolds)

chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1))) #numbering the chromosomes

#joining chromosome num to our vcf
vcf.filt.indSNPMiss.70.div.Hs <- left_join(chrnum, vcf.filt.indSNPMiss.70.div, join_by(chr.main==CHROM)) 

#filtering out Gst(Fst) values <0 and creating "SNP" column w/ both chromo # and bp position
vcf.filt.indSNPMiss.70.div.Hs <- vcf.filt.indSNPMiss.70.div.Hs %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main,"_",POS))

#making these columns numeric, not characters
vcf.filt.indSNPMiss.70.div.Hs$V2 = as.numeric(vcf.filt.indSNPMiss.70.div.Hs$V2)
vcf.filt.indSNPMiss.70.div.Hs$POS = as.numeric(vcf.filt.indSNPMiss.70.div.Hs$POS)

#note to self: all this chromosome stuff might be completely unnecessary. consider removing when cleaning up the script
#this is where 02 moves on to Hs
names(vcf.filt.indSNPMiss.70.div.Hs) #Hs values are in columns 4-9

#taking a peek at what our Hs values look like for each region w/ a histogram
vcf.filt.indSNPMiss.70.div.Hs %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) +
  geom_histogram(position="identity", alpha=0.5, bins=50) +
  labs(title="Genome-wide expected heterozygosity (Hs)",fill="Regions",
       x="Gene diversity (Hs) within Regions", y="Counts of SNPs")

#generating table with the summary statistics we need!!!

loci <- nrow(vcf.filt.indSNPMiss.70.div.Hs) #15748 loci at this stage
#to get number of loci where Hs = 0, we'll do loci - the N we calculate after filtering values that equal 0

Hs_table.70 <- vcf.filt.indSNPMiss.70.div.Hs %>%
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  filter(value!=0) %>% 
  summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n(), N_0Hs=loci-n()) 

## A tibble: 6 x 5
#name   avg_Hs StdDev_Hs  N_Hs N_0Hs
#<chr>   <dbl>     <dbl> <int> <int>
#1 Hs_CEU  0.174     0.158  8326  7422
#2 Hs_NE   0.105     0.139 14344  1404
#3 Hs_NEU  0.199     0.155  7438  8310
#4 Hs_PNW  0.147     0.158 10159  5589
#5 Hs_SEU  0.244     0.152  5144 10604
#6 Hs_WEU  0.195     0.155  7751  7997

write.csv(Hs_table.70, "~/projects/eco_genomics/population_genomics/outputs/Hs_table.70_noZeros.csv",
          quote=F,
          row.names=F)

#incase i wanna pick back up where i was
#vcf.filt.indSNPMiss.70 <- read.vcfR("/outputs/vcf.70.filtered.vcf.gz")

#thinning out, eliminating SNPs within 500bp
vcf.filt.indSNPMiss.70.thin <- distance_thin(vcf.filt.indSNPMiss.70, min.distance = 500)

#rereading metadata
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)

#including only the filtered individuals
pcameta2 <- meta[meta$id %in% colnames(vcf.filt.indSNPMiss.70@gt[, -1]) , ]

dim(pcameta2)

#sneaking under the door (cuz github doesn't like big files)
setwd("~/projects/eco_genomics/population_genomics/")
write.vcf(vcf.filt.indSNPMiss.70.thin, "outputs/vcf.70.filtered.thinned.vcf.gz")

#hide uncompressed vcf file (too big for github!) outside of our repo

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf.70.filtered.thinned.vcf.gz > ~/vcf.70.filtered.thinned.vcf") 
#system() basically lets you do stuff like you would in command line

geno.70 <- vcf2geno(input.file="/gpfs1/home/c/s/cstavros/vcf.70.filtered.thinned.vcf",
                    output.file = "outputs/vcf.70.filtered.thinned.geno")

PCA.70 <- LEA::pca("outputs/vcf.70.filtered.thinned.geno", scale=TRUE)

ggplot(as.data.frame(PCA.70$projections),
       aes(x=V1, y=V2, color=pcameta2$region, shape=pcameta2$continent)) +
  geom_point(alpha=.5) +
  labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region",shape="Continent")

#okay, it's weird and upside down from the last one, but it'll do for the moment.
#apparently john had the same problec, i'll see if he has any wisdom about this

#yay! saving rn
ggsave("figures/PCA.70_PC1vPC2.pdf", width=6, height=6, units="in")




####yipee, moving onto .65!#####




vcf.filt.indMiss.65 <- filter_biallelic(vcf.filt.indMiss.65) #removing non-biallelic
vcf.filt.indMiss.65 <- min_mac(vcf.filt.indMiss.65, min.mac = 1) #mac = Minor Allele Count, min.mac=1 excludes all minor alleles

vcf.filt.indSNPMiss.65 <- missing_by_snp(vcf.filt.indMiss.65, cutoff=0.5) #filtering missingness by SNP

vcf.filt.indSNPMiss.65 #taking a peek at it:

#***** Object of Class vcfR *****
#  545 samples
#20 CHROMs
#16,106 variants
#Object size: 100.8 Mb
#27.4 percent missing data
#*****        *****         *****

#creating a file of depth
DP2.65 <- extract.gt(vcf.filt.indSNPMiss.65,
                     element = "DP",
                     as.numeric = T)
write.vcf(vcf.filt.indSNPMiss.65,
          "~/projects/eco_genomics/population_genomics/outputs/vcf.65.filtered.vcf.gz")

#finding genome-wide diversity for each region (Hs) and its Std Dev:

#importing metadata again for next steps

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta) # our vcf file has 545 samples, but
dim(meta) # meta has 629 inds

newmeta2 <- meta[meta$id %in% colnames(vcf.filt.indSNPMiss.65@gt[,-1]),] #filtering meta to consist only of the inds remaining post-filtering
dim(newmeta2) #now it matches!

#calculate diversity stats using genetic_diff function in vcfR
vcf.filt.indSNPMiss.65.div <- genetic_diff(vcf.filt.indSNPMiss.65,
                                           pops=as.factor(newmeta2$region),
                                           method = "nei")

str(vcf.filt.indSNPMiss.65.div) #peeking at structure

chr.main <- unique(vcf.filt.indSNPMiss.65.div$CHROM)[1:8] #taking the first 8 chromosomes (ignoring the tiny scaffolds)

chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1))) #numbering the chromosomes

#joining chromosome num to our vcf
vcf.filt.indSNPMiss.65.div.Hs <- left_join(chrnum, vcf.filt.indSNPMiss.65.div, join_by(chr.main==CHROM)) 

#filtering out Gst(Fst) values <0 and creating "SNP" column w/ both chromo # and bp position
vcf.filt.indSNPMiss.65.div.Hs <- vcf.filt.indSNPMiss.65.div.Hs %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main,"_",POS))

#making these columns numeric, not characters
vcf.filt.indSNPMiss.65.div.Hs$V2 = as.numeric(vcf.filt.indSNPMiss.65.div.Hs$V2)
vcf.filt.indSNPMiss.65.div.Hs$POS = as.numeric(vcf.filt.indSNPMiss.65.div.Hs$POS)

#note to self: all this chromosome stuff might be completely unnecessary. consider removing when cleaning up the script
#this is where 02 moves on to Hs
names(vcf.filt.indSNPMiss.65.div.Hs) #Hs values are in columns 4-9

#taking a peek at what our Hs values look like for each region w/ a histogram
vcf.filt.indSNPMiss.65.div.Hs %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) +
  geom_histogram(position="identity", alpha=0.5, bins=50) +
  labs(title="Genome-wide expected heterozygosity (Hs)",fill="Regions",
       x="Gene diversity (Hs) within Regions", y="Counts of SNPs")

#generating table with the summary statistics we need!!!

loci <- nrow(vcf.filt.indSNPMiss.65.div.Hs) #15748 loci at this stage
#to get number of loci where Hs = 0, we'll do loci - the N we calculate after filtering values that equal 0

Hs_table.65 <- vcf.filt.indSNPMiss.65.div.Hs %>%
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  filter(value!=0) %>% 
  summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n(), N_0Hs=loci-n()) 

# A tibble: 6 Ã 5
#name   avg_Hs StdDev_Hs  N_Hs N_0Hs
#<chr>   <dbl>     <dbl> <int> <int>
#1 Hs_CEU  0.176     0.158  8416  7634
#2 Hs_NE   0.106     0.139 14595  1455
#3 Hs_NEU  0.201     0.155  7550  8500
#4 Hs_PNW  0.147     0.157 10352  5698
#5 Hs_SEU  0.251     0.150  5118 10932
#6 Hs_WEU  0.196     0.155  7896  8154

write.csv(Hs_table.65, "~/projects/eco_genomics/population_genomics/outputs/Hs_table.65_noZeros.csv",
          quote=F,
          row.names=F)




#putting the bookmark in here, i'll be back after class!



#incase i wanna pick back up where i was
#vcf.filt.indSNPMiss.65 <- read.vcfR("/outputs/vcf.65.filtered.vcf.gz")

#thinning out, eliminating SNPs within 500bp
vcf.filt.indSNPMiss.65.thin <- distance_thin(vcf.filt.indSNPMiss.65, min.distance = 500)

#rereading metadata
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)

#including only the filtered individuals
pcameta2 <- meta[meta$id %in% colnames(vcf.filt.indSNPMiss.65@gt[, -1]) , ]

dim(pcameta2)

#sneaking under the door (cuz github doesn't like big files)
setwd("~/projects/eco_genomics/population_genomics/")
write.vcf(vcf.filt.indSNPMiss.65.thin, "outputs/vcf.65.filtered.thinned.vcf.gz")

#hide uncompressed vcf file (too big for github!) outside of our repo

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf.65.filtered.thinned.vcf.gz > ~/vcf.65.filtered.thinned.vcf") 
#system() basically lets you do stuff like you would in command line

geno.65 <- vcf2geno(input.file="/gpfs1/home/c/s/cstavros/vcf.65.filtered.thinned.vcf",
                    output.file = "outputs/vcf.65.filtered.thinned.geno")

PCA.65 <- LEA::pca("outputs/vcf.65.filtered.thinned.geno", scale=TRUE)

ggplot(as.data.frame(PCA.65$projections),
       aes(x=V1, y=V2, color=pcameta2$region, shape=pcameta2$continent)) +
  geom_point(alpha=.5) +
  labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region",shape="Continent")

#okay, it's weird and upside down from the last one, but it'll do for the moment.
#apparently john had the same problec, i'll see if he has any wisdom about this

#yay! saving rn
ggsave("figures/PCA.65_PC1vPC2.pdf", width=6, height=6, units="in")