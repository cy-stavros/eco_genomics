# Coding and Data notes for the transcriptomics module

## Author: Cy Stavros

### 10-10-2024 - Cleaning up Copepod data and generating a PCA!

- After some trouble getting DESeq2 to load tuesday, we're back today

- loading in DESeq and ggplot

- bringing in the table with out counts

- note that transcript counts are to transcriptomics what SNPs are to population genomics (maybe they should call it SNPomics lol)

- to begin with, we averaged around 18mil transcripts per sample group 

- note that this data has a serious right skew. most trascripts are lowly expressed, a few have a lottttt of expression (lotsa transcripts)

- filtered to exclude transcripts with less than 10 reads

- filtered transcripts that appear in less than 15 samples

- running the DESeq program!

- generating a quick pca, looks like clustering is mostly happening between finaltemp

### 10-10-2024 - Taking DESeq data and making MAplot, volcano plot, and heatmap

- starting by rerunning previous script (01) to get dds back

- from that, getting a smaller df with just the results we're interested in

- picking out the highest transcript to graph + compare between dev temps

- this helps orient us, now we know a positive logfold2 means upreg in d22/downreg in d18

-next making MA plot

- this helps us visualize the overdispersion, also the amount of upreg/downreg on each side.

### 10-17-2024 - Doing contrast stuff
