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

#first trying to see if i can just subset the file for just european data

geno <- vcf2geno(input.file="/gpfs1/home/c/s/cstavros/vcf_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned.geno")
#somehow subset this^

EuCentAdmix <- snmf("outputs/vcf_final.filtered.thinned.geno",
                  K=1:10, #giving it a range of K values
                  entropy = T, #asks it to do cross-validation
                  repetitions = 3,
                  project = "new") # if you're adding to this analysis later, you could choose project="
#use subsetted file in snmf function

#alrighty it looks like we might just have to use structure after all

######### Converting vcf to structre ########

#trying vcf2others

library(vcf2others)
#no return message, but it looks like it loaded

#okay i need to use the vcf2structure() function

thinnedvcf <- read.vcfR("/gpfs1/home/c/s/cstavros/vcf_final.filtered.thinned.vcf")

#as i understand itvcf2structre needs metadata in a specific format, 
#I'm not super sure how to make that happen rn

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)

meta2 <- meta[meta$id %in% colnames(thinnedvcf@gt[, -1]) , ]

dim(meta2)

#pasting stuff I don't think is gonna work from chatgpt:
#(i was right lol)

#checking vcf and meta align
vcf_samples <- colnames(thinnedvcf@gt)  # thinnedvcf is the VCF object

# Check the sample names in your metadata
metadata_samples <- meta2$id  # Adjust this column name if it's different

# Check if sample names match
identical(vcf_samples, metadata_samples)

#alright trying with slightly less dumb gpt answer:

ind_pop <- meta2$region
#keep_pop <- unique(ind_pop)
#keep_pop <- rep(TRUE, length(ind_pop))

#after arguing w/ chatgpt for like 30 minutes and then actually 
#just finding the source code:

# Make sure ind_pop is a factor (you might already have this)
ind_pop <- factor(ind_pop)

# Create a factor for keep_pop with all unique population names (same as ind_pop)
keep_pop <- factor(unique(ind_pop))





vcf2structure(
  vcf = thinnedvcf,           # Your VCF object
  ind_pop = ind_pop,       # Vector of population info from metadata
  keep_pop = keep_pop,     # Whether to keep all populations or filter
  inc_missing = TRUE,      # Whether to include missing genotypes
  out_file = "thinned.str",  # Output file name
  method = "S"             # Structure format method
)


#thank god that finally worked! checking the formatting now

# Read the .str file into R as a data frame
#my_structure_file <- read.table("thinned.str", header = FALSE, 
                                #sep = " ")
#yields a df with one column
#my_structure_file <- read.delim("thinned.str", header = FALSE)
#yields the V1, V2 thingy

#my_structure_file <- read.fwf("thinned.str", widths = rep(10, 100))
#doesn't even work

my_structure_file <- read.table("thinned.str", header = FALSE, sep = "\t")


# Preview the first few rows
head(my_structure_file)

#okay none of this crap is working, but I downloaded the file it and it looks 
#pretty normal in my notepad. wrapping up for the night.
#hold on i just looked at it in the environment tab and it looks super normal, 
#i think the head function is just buggin bc there are so many dang columns

#duplicating the metadata and stitching IDs, also flag info to the structure file
#[ignore all, just found out vcf2structure reorders it all]
# #meta2_repeated <- meta2 %>%
#   slice(rep(1:n(), each = 2))
# 
# region_repeated <- meta2_repeated %>%
#   select(region)

regions <- c("PNW", "NE", "WEU", "CEU", "NEU", "SEU")
#in the order they appear in .str file!

regiondf <- regions[my_structure_file$V1]

popflagdf <- ifelse(my_structure_file$V1 %in% 3:6, 1, 0)

structure_with_meta <- bind_cols(regiondf, popflagdf, my_structure_file)
#sucess! overwriting old file
write.table(structure_with_meta, file = "thinnedwithmeta3.str", append = TRUE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


#i fear i have to put the region info after the population index in the structure file
#the order of columns should be:
#pop index (column 1 of structure file)
#popflag
#region (characters)
#everything else

structure_with_meta <- bind_cols(my_structure_file[,1], popflagdf, regiondf, my_structure_file[,-1])
#sucess! overwriting old file
write.table(structure_with_meta, file = "thinnedwithmeta4.str", append = TRUE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


structure_and_meta <- read.table("thinnedwithmeta.str", header = TRUE, sep = "\t")


#okay i thought i was good, but the popflag isn't correct for WEU, also idk 
#what the ascending 1, 2, 3, 4, 5, 6 column is after popflag?
#okay apparently these are population assignments (which group pops together),
#so my meta is junk. row order is not preserved at all.

#okay here are what the indices mean:
as_tibble(keep_pop)

# 1 PNW  
# 2 NE   
# 3 WEU  
# 4 CEU  
# 5 NEU  
# 6 SEU 

#looks like i won't be able to salvage individual id, but i'll go in and 
#rebind with correct pop labels and popflags
#okay that's all done, but i wanna use LEA to estimate K just within the european subpops


####### estimating european k ########

setwd("~/projects/eco_genomics/population_genomics/")
fullvcf <- read.vcfR("outputs/vcf_final.filtered.vcf.gz")

fullvcf.thin <- distance_thin(fullvcf, min.distance = 500)

#okay i have to somehow only get the european samples out of here, but the 
#metadata is in the meta file :/

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta2 <- meta[meta$id %in% colnames(fullvcf.thin@gt[, -1]) , ]


europe <- c("WEU", "CEU", "NEU", "SEU")

europemeta <- meta2 %>% filter(region %in% europe)
europeids <- europemeta$id

europevcf.thin <- fullvcf.thin[, c("CHROM", "POS", europeids)] #that didn't work for some reason

colnames(fullvcf.thin@gt)[-1] #these all look good

europe_indices  <- which(colnames(fullvcf.thin@gt) %in% europeids)

europevcf.thin <- fullvcf.thin[, c(1, europe_indices)]

x <- colnames(europevcf.thin@gt)[-1]
#okay looks good (i think) 

europeids == x
#okay not all the same, but maybe that just has to do with ordering?

setequal(x, europeids)
#okay they are the same!!! HUGE

#okay now to determine K.

write.vcf(europevcf.thin, "outputs/europevcf.thin.vcf.gz")


#hide uncompressed vcf file (too big for github!) outside of our repo

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/europevcf.thin.vcf.gz > ~/europevcf.thin.vcf") 
#system() basically lets you do stuff like you would in command line



eurogeno <- vcf2geno(input.file="/gpfs1/home/c/s/cstavros/europevcf.thin.vcf",
                 output.file = "outputs/europevcf.thin.geno")
#"3 line(s) were removed because these are not SNPs."
#idk if i need to worry about that, probably not

EuroCentAdmix <- snmf("outputs/europevcf.thin.geno",
                  K=1:10, #giving it a range of K values
                  entropy = T, #asks it to do cross-validation
                  repetitions = 3,
                  project = "new")

plot(EuroCentAdmix)
#it looks like 4 might be a good number, it's lowest at 5 but there are
#4 sampled pops (see screenshot) (11/19/24)
#making barplot!

myK = 6 #creates little placeholder we can modify later

CE = cross.entropy(EuroCentAdmix, K=myK)
best = which.min(CE)
#not sure why my results aren't printing to screen rn

myKQ = Q(EuroCentAdmix, K=myK, run=best)
#results in a big matrix, where each column is a different group

myKQmeta = cbind(myKQ, europemeta) #to use cbind(), you have to have the same number of rows in the same order!

my.colors = c("blue4", "gold", "tomato", "lightblue", "olivedrab", "purple") #setting up little color palette

myKQmeta = as_tibble(myKQmeta) %>% 
  group_by(continent) %>% 
  arrange(region, pop, .by_group = T) #says: first group by continent, and then group within contnent by region and pop

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
