#This is the script no.1 required for downloading the TCGA data with TCGAbiolinks
BiocManager::install("edgeR")
BiocManager::install("genefilter")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data", force = TRUE)
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks", force = TRUE)
#devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
# Load packages
library("TCGAbiolinks")
library("SummarizedExperiment")

setwd("/home/ieva/tcga-ov-data/")
GDCprojects = getGDCprojects()
GDCprojects[c("project_id", "name")] #I choose ovarian cancer project TCGA-OV

TCGAbiolinks:::getProjectSummary("TCGA-OV")

query_TCGA1 = GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling")
#Find out about all transcriptomic data in TCGA project
View(getResults(query_TCGA1)) 
#Build your query
query_TCGA = GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts", 
  access = "open" )
#Chek your query
View(getResults(query_TCGA))
# Make results as table
lihc_res = getResults(query_TCGA) 
#Chek our your data
head(lihc_res) # data of the first 6 patients.
colnames(lihc_res) # columns present in the table
unique(lihc_res$sample_type)
summary(factor(lihc_res$sample_type))
#Download your data (do not forget to setwd to data folder, creates a big folder)
GDCdownload(query = query_TCGA, files.per.chunk = 200) 
#Prepare your data
tcga_data <- GDCprepare(query_TCGA,
                        summarizedExperiment = TRUE)
#Chek out final downloaded data
dim(tcga_data)
colnames(colData(tcga_data))
# Save clinical data to separate file
pheno <- as.data.frame(colData(tcga_data)) #impotant to save clinical data
# this is on original col data, but could be done on pheno:
#count some most important clinical data
table(tcga_data@colData$figo_stage)
table(tcga_data@colData$vital_status)
#Look at count data:
dim(assay(tcga_data))
head(assay(tcga_data)[,1:10])
head(rowData(tcga_data))     # ensembl id and gene id of the first 6 genes.
#############################################################################
#final saved TCGA-OC data: full file RangedSummarizedExperiment tcga_data and clinical data (colData)
# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again
saveRDS(object = tcga_data,
        file = "tcga_data.RDS",
        compress = FALSE)

saveRDS(object = pheno,
        file = "pheno.RDS",
        compress = FALSE)

