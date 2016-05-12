# Rtool-SBG
This tool can be used in Seven Bridges Platform, in order to compare the differential expression results from STAR -- Cufflinks workflow and STAR - HTSeq-count - DESeq2 workflow. 


##### DE_comparison.R
R script performing pairwise comparisons for 3 groups of samples (HER2, TNBC, and NonTNBC).

##### Dockerfile
Dockerfile to build the docker image for this tool. The docker image can also be pulled from jzhao31/de-comparison:v1

##### DE_comparison_Tool.json
Json file for the tool. Users can import this Json file to Seven Bridges SDK to build the tool.



### test files
##### gene_exp.diff
Differential gene expression result from STAR -- Cufflinks workflow

##### DESeq2_results_TNBC_Tumor_vs_HER2_Tumor.csv
Differential gene expression result from STAR - HTSeq-count - DESeq2 workflow. It compares TNBC tumor samples with HER2 tumor samples.

##### DESeq2_results_TNBC_Tumor_vs_NonTNBC_Tumor.csv
Differential gene expression result from STAR - HTSeq-count - DESeq2 workflow. It compares TNBC tumor samples with NonTNBC tumor samples.

##### DESeq2_results_NonTNBC_Tumor_vs_HER2_Tumor.csv
Differential gene expression result from STAR - HTSeq-count - DESeq2 workflow. It compares HER2 tumor samples with NonTNBC tumor samples.
