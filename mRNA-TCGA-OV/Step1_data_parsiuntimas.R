#This is the script no.1 required for downloading the TCGA data with TCGAbiolinks
# Load packages
library("TCGAbiolinks")
library("SummarizedExperiment")
setwd("~/rprojects/TCGA-OV-data/") #wsl

## Build your query
query_TCGA = GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts", 
  access = "open" )

##Download your data (do not forget to setwd to data folder, creates a big folder)
GDCdownload(query = query_TCGA, files.per.chunk = 200) 
##Prepare your data
tcga_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)

## Save clinical data to separate file 
pheno <- as.data.frame(colData(tcga_data)) 
#I have since found a more full clinical data on XENA

#############################################################################
#final saved TCGA-OC data: full file RangedSummarizedExperiment tcga_data and clinical data (colData)
#might need to remove all of the downloaded files, because it takes a tone of space
saveRDS(object = tcga_data,
        file = "tcga_data.RDS",
        compress = FALSE)

saveRDS(object = pheno,
        file = "pheno.RDS",
        compress = FALSE)

