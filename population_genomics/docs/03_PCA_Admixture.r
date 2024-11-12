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

meta2 <- meta[meta$id %in% colnames(vcf@gt[, -1]) , ]

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


#10/1/24
#Now, we will run admixture analyses and create plots
#for admixture, we're going to use the LEA R package
#the function inside LEA is called 'snmf'
#snmf is very fast, much faster than bayesian structure analyses, but yields comparable results

CentAdmix <- snmf("outputs/vcf_final.filtered.thinned.geno",
                  K=1:10, #giving it a range of K values
                  entropy = T, #asks it to do cross-validation
                  repetitions = 3,
                  project = "new") # if you're adding to this analysis later, you could choose project="

plot(CentAdmix) #plot of the cross-validatation experiment - if it doesn't start to level off, you should probably increase K
#look for the elbow in the plot! consider a range within the elbow (~3-5 for this data)
#the graph will loop back up as K increases, some sample at the bottom of the U
#basically, you're looking for the least complex model that explains your data reasonably well

par(mfrow=c(2,1)) #creates stacked plots
plot(CentAdmix, col="blue4", main="SNMF")
plot(CentPCA$eigenvalues[1:10], ylab="Eigenvalues", xlab="Number of PCs", col="blue4", main="PCA") #notice the similarity in graphs!
dev.off() #resets to non-stacked plots

myK = 5 #creates little placeholder we can modify later

CE = cross.entropy(CentAdmix, K=myK)
best = which.min(CE)
#not sure why my results aren't printing to screen rn

myKQ = Q(CentAdmix, K=myK, run=best)
#results in a big matrix, where each column is a different group

myKQmeta = cbind(myKQ, meta2) #to use cbind(), you have to have the same number of rows in the same order!

my.colors = c("blue4", "gold", "tomato", "lightblue", "olivedrab") #setting up little color palette

myKQmeta = as_tibble(myKQmeta) %>% 
  group_by(continent) %>% 
  arrange(region, pop, .by_group = T) #says: first group by continent, and then group within contnent by region and pop

pdf("figures/Admixture_K4.pdf",width=10, height=5)
barplot(as.matrix(t(myKQmeta[ , 1:myK])),
        border=NA,
        space=0,
        col=my.colors[1:myK],
        xlab="Geographic regions",
        ylab="Ancestry proportions",
        main=paste0("Ancestry Matrix K=", myK)) #creates a label with whatever myK is in the title

axis(1,
     at=1:length(myKQmeta$region),
     labels=myKQmeta$region,
     tick=F,
     cex.axis=0.5,
     las=3) #turns labels on side (allegedly)
dev.off()
#notice how this compares to the PCA! distinct cluster in PNW, also CEU, decent admixture for everything else

######experimenting for my final project#####

#trying to train k model on european data, and then use it to group new world data

