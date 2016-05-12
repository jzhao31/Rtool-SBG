# Rtool-SBG
This tool can be used in Seven Bridges Platform, in order to compare the differential expression results from STAR -- Cufflinks workflow and STAR - HTSeq-count - DESeq2 workflow. 
  
  
##### DE_comparison.R
R script performing pairwise comparisons for 3 groups of samples (HER2, TNBC, and NonTNBC).

##### Dockerfile
Dockerfile to build the docker image for this tool. The docker image can also be pulled from jzhao31/de-comparison:v1

##### DE_comparison_Tool.json
Json file for the tool. Users can import this Json file to Seven Bridges SDK to build the tool.



### test files
##### CuffDiff_gene_exp_small.diff
A test file of differential gene expression result from STAR -- Cufflinks workflow

##### DESeq2_TNBC_HER2_small.csv
A test file of differential gene expression result from STAR - HTSeq-count - DESeq2 workflow. It compares TNBC tumor samples with HER2 tumor samples.

##### DESeq2_TNBC_NonTNBC_small.csv
A test file of differential gene expression result from STAR - HTSeq-count - DESeq2 workflow. It compares TNBC tumor samples with NonTNBC tumor samples.

##### DESeq2_NonTNBC_HER2_small.csv
A test of differential gene expression result from STAR - HTSeq-count - DESeq2 workflow. It compares HER2 tumor samples with NonTNBC tumor samples.
