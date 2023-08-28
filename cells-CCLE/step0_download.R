#this script will read the gtex to a managable data format
#this is my home wls because it can tar gz files and also because conda!
library(phantasus)
library(SummarizedExperiment)
setwd("home/ieva/rprojects/TCGA-OV-data/GTEX")
CLLE <- read.gct("CCLE_RNAseq_genes_counts_20180929.gct")
CLLE_COUNTS <- exprs(CLLE)
saveRDS(CLLE_COUNTS, "CLLE_COUNTS.RDS")
