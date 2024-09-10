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

-   If you are in the correct directory w/ a specific file or directory, you can use tab to auto-complete the file name or directory/
