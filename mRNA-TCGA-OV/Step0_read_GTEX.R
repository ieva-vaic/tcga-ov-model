#STEP0 read GTEX file
#this script will read the gtex to a manageable data format
#this is my home wsl because it can tar gz files.
#data downloaded form https://www.gtexportal.org/home/datasets; GTEx Analysis V8 release; counts by tissue
library(phantasus)
library(SummarizedExperiment)
setwd("/home/ieva/rprojects/TCGA-OV-data/GTEX")
gtex <- read.gct("gene_reads_2017-06-05_v8_ovary.gct")
gtex_counts <- exprs(gtex)
saveRDS(gtex_counts, "gtex_counts.RDS")
