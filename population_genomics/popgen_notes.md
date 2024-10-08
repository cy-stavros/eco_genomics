# Coding and Data notes for the population genomics module

## Author: Cy Stavros

### 09-10-2024 - Intro to *Centaurea* GBS data and working with VCF files

We'll be analyzing the GBS data from 3 regions (EU, NE, PNW) starting today with variant call format files (VCFs)

Notes:

-   Lots to take note of!

-   Paired-end data is expensive and somewhat redundant due to loci linkage, depending on the study it might be better to spend money on getting more fragments than reading from the other end.

-   Don't forget to pull from github at the start of your R sesh!

-   Q scores are encoded in fastq files using characters, if you see a letter (I is 40!!), that's good, other characters (#@%) not so much

-   zcat command lets you take a peek into a these huge files using something like head

-   If you are in the correct directory w/ a specific file or directory, you can use tab to auto-complete the file name or directory

### 09-12-2024

- Worked on viewing sequence data and Q-scores in command line

- Most of this was in command line, so I don't have inline notes for it unfortunately

- Command line allows us to take a peak into these HUGE files to see what the actual sequence reads look like

-unsurprisingly, different areas have different numbers of reads (depth), and sometimes these reads are conflicting

- conflicting reads can be a result of a genuine mistake in the sequencer, or the result of an individual being heterozygous at this site. (Q-scores and proportion of the reads can help parse this)

### 09-17-2024

- worked on filtering the VCF file from 3 directions: depth, missingness, and low-frequency alleles

- low depth (hereforth DP) results in inaccurate calls, while excessively high DP can be the result of assembly errors or paralogy (reads from duplicated genes mapping to the same position)

- Missingness is filtered in two aspects, individual-level excludes individuals that don't have data for many of the observed snps, while snp-level excludes snps that only apply to a few individuals (not really sure why analyses can't account for this missingness though, cuz it seems to throw out a lot of data (unless that data is inaccurate trash?))

- (putative) low-frequency alleles (like less than 1%) are often just fake articfacts of bad reads, so those are filtered out.

- also note, you typically want to only handle biallelic snps (b/c from my understanding, more than one allele at any one position is exceedingly rare within a species (b/c you'd have to get a mutation at the same position twice!)) this might apply for larger, phylogenetic analyses but we don't have to worry about that right now. 

### 09-19-2024

- Taking filtered vcf from 01 VCF filtering and generating a manhattan plot in 02_Diversity_Differentiation

- first, wrangling data to join in the metadata with info like pop, region (as always, see script for more notes)

-breaking down data into chromosomes, removing scaffolds

- using genetic_diff to get Fst (here as Gst, but the two are the same for biallelic SNPs)

- then, generating manahattan plot!!

### 09-24-2024

- starting new script (03) to generate a PCA

- the line " options(bitmapType = "cairo") " does something with the way the vacc generates figures, I'm not quite sure what but the pca wasn't generating properly without it.

- PCAs need to be filtered for linkage disequlibrium b/c close snps will not likely interact as separate variables

- "thinned" file (w/ no loci within 500 bp)

- (as always, more notes on this stuff in the script)

- We've been working with a compressed version of our file (.gz) but we need to open it outside of our repo to generate the pca file

- To "hide" this from github, we're using the system() command which from my understanding basically lets you treat rstudio like command line (refer back to line 29)

- Once we do the proper data wrangling to create a pca file, it can be graphed using baseR plot()

- Ran out of time today, but we'll come back thurs and make a nicer looking one with ggplot


### 09-26-2024

-Picking up where we left off last time with the PCA plot.

- Note - if returning to a pca you ran last time, you can load the results without running it again by using load.pcaProject()

- This time we're making in in ggplot which gives us more customizability

- you can also change it from V1 and V2 to 2 and 3

- cropping the graph and adjusting alpha (opacitiy) helps you see things you might not see otherwise

- also learned that plot() of a pca gives you a scree plot

- then discussed Admixture Analysis (structure), but didn't get to coding it yet. Looking forward to doing that soon!

### 10-1-24

- Starting with admixture analysis (03) and then moving on to scanning for selection (04)

- snmf takes much less time, computational power than bayesean analyses but yields similar results

- (03 script is heavily annotated, see it for more detailed notes)

- You provide a range of possible K values (we did 1:10)

- the plot of CentAdmix helps us determine an optimal K value (in the "elbow" of the graph)

- the shape of this graph is similar to the amount of variation explained by each subsecquent PC

- generated figures showinf the results of this AA.

- for 04, first regenerating the vcfr and metadata files

- creating manahattan plot

- cont



