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

- as always, check out the inline notes!

- next, transfering results and p-vals only into a new df to make a volcano plot

- cool! in addition to the .05 p-value threshold, a threshold of a log fold 2 change of >1 is added (meaning diff exp must be at least 2x in magnitude!!)

- in this case, it looks like there is more upreg in 22 / downreg in 18

- next, looking at things through a heatmap

- columns are samples, rows are genes

- note that the phenograms on the side do not imply evolutionary relatedness of the transcripts,
but could be indicative transcripts working in the same pathway

### 10-17-2024 - Doing contrast stuff

- within the DESeq file (dds), we're not setting up groups by each unique comparison of factors (e.g. d18a33 vs d18a28)

- by grouping them like this we can compare within and between factors

- first, we're comparing gene expression at baseline between developmental temps ("18BASE" vs "22BASE")

- to do so, making results file where NA's are filtered out, pvals must be over .05, and ordering based on pvals

- making a list of transcripts (and MA plotting them) that a differentially expressed

- rinse & repeating for A28 and A33

- looking at the length of the dfs we made (ie the number of transcripts in each comparison)

- using nested intersect()s to find the overlap in all 3 comparisons

- check venn diagram in physical notes if this isn't clicking!

### 10-22-2024 - continuing w/ euler plots and then making plots that compare log fold change (idk what they're called)
### 10-24-2024 - continuing with yesterday's plots and starting wgcna network analysis