# Coding and Data notes for the population genomics module

## Author: Cy Stavros

### 09-10-2024 - Intro to *Centaurea* GBS data and working with VCF files

We'll be analyzing the GBS data from 3 regions (EU, NE, PNW) starting today with variant call format files (VCFs)

Notes:

-   Lots to take note of!

-   Paired-end data is expensive and somewhat redundant due to loci linkage, depending on the study it might be better to spend money on getting more fragments than reading from the other end.

-   Don't forget to pull at the start of your R sesh!

-   Q scores are encoded in fastq files using characters, if you see a letter (I is 40!!), that's good, other characters not so much

-   zcat command lets you take a peek into a these huge files using something like head

-   If you are in the correct directory w/ a specific file or directory, you can use tab to auto-complete the file name or directory

### 09-12-2024

- Worked on viewing sequence data and Q-scores in command line

### 09-17-2024

- worked on filtering the VCF file from 3 directions: depth, missingness, and low-frequency alleles

- low depth (hereforth DP) results in inaccurate calls, while excessively high DP can be the result of assembly errors or paralogy (reads from duplicated genes mapping to the same position)

- Missingness is filtered in two aspects, individual-level excludes individuals that don't have data for many of the observed snps, while snp-level excludes snps that only apply to a few individuals (not really sure why analyses can't account for this missingness though, cuz it seems to throw out a lot of data (unless that data is inaccurate trash?))

- (putative) low-frequency alleles (like less than 1%) are often just fake articfacts of bad reads, so those are filtered out.

- also note, you typically want to only handle biallelic snps (b/c from my understanding, more than one allele at any one position is exceedingly rare within a species (b/c you'd have to get a mutation at the same position twice!)) this might apply for larger, phylogenetic analyses but we don't have to worry about that right now. 
