# Coding and Data notes for the transcriptomics module

## Author: Cy Stavros

### 10-10-2024 - Cleaning up Copepod data and generating a PCA!

- After some trouble getting DESeq2 to load tuesday, we're back today

- loading in DESeq and ggplot

- bringing in the table with out counts

- note that transcript counts are to transcriptomics what SNPs are to population genomics (maybe they should call it SNPomics lol)

- to begin with, we averaged around 18mil transcripts per sample group 

- note that this data has a serious right skew. most trascripts are lowly expressed, a few have a lottttt of expression (lotsa transcripts)

- filtered to exclude transcripts with less than