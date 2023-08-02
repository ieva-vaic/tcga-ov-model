BiocManager::install("edgeR")
BiocManager::install("genefilter")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data", force = TRUE)
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks", force = TRUE)
# BiocManager::install("TCGAbiolinks") # and abobe doesnt work, installing related packages
#install.packages(c("boot", "class", "cluster", "codetools", "foreign", "KernSmooth", "lattice", "MASS", "Matrix", "mgcv", "nlme", "nnet", "rpart", "spatial", "survival")) # did not work instaling packages thus installing dependaciens ‘biomaRt’, ‘GenomicRanges’, ‘XML’, ‘rvest’, ‘SummarizedExperiment’, ‘xml2’, ‘httr’
#BiocManager::install("biomaRt")
# ‘GenomicRanges’, ‘XML’, ‘rvest’, ‘SummarizedExperiment’, ‘xml2’, ‘httr’ 
#BiocManager::install(c("XML", "rvest", "SummarizedExperiment", "xml2", "httr" ))
#withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true"), remotes::install_github('BioinformaticsFMRP/TCGAbiolinks'))
devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")

# Load packages
library("TCGAbiolinks")
# library("limma")
# library("edgeR")
# library("glmnet")
# library("factoextra")
# library("FactoMineR")
# library("caret")
library("SummarizedExperiment")
# library("gplots")
# library("survival")
# library("survminer")
# library("RColorBrewer")
# library("gProfileR")
# library("genefilter")

setwd("/home/ieva/tcga-ov-data/")
GDCprojects = getGDCprojects()
GDCprojects[c("project_id", "name")] #aš noriu TCGA-OV

TCGAbiolinks:::getProjectSummary("TCGA-OV")

query_TCGA1 = GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling")
View(getResults(query_TCGA1)) #čia sužinau apie visus OV mėginius

query_TCGA = GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts", 
  access = "open" )

View(getResults(query_TCGA))

lihc_res = getResults(query_TCGA) # make results as table
head(lihc_res) # data of the first 6 patients.
colnames(lihc_res) # columns present in the table
unique(lihc_res$sample_type)
summary(factor(lihc_res$sample_type))

GDCdownload(query = query_TCGA, files.per.chunk = 200) #downloads GDCdata files, I moved it out afterwards

tcga_data <- GDCprepare(query_TCGA,
                        summarizedExperiment = TRUE)

dim(tcga_data)

colnames(colData(tcga_data))
pheno <- as.data.frame(colData(tcga_data))

table(tcga_data@colData$figo_stage)

table(tcga_data@colData$vital_status)

dim(assay(tcga_data))
head(assay(tcga_data)[,1:10])
head(rowData(tcga_data))     # ensembl id and gene id of the first 6 genes.

#final saved TCGA-OC data

# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again
saveRDS(object = tcga_data,
        file = "tcga_data.RDS",
        compress = FALSE)

saveRDS(object = pheno,
        file = "pheno.RDS",
        compress = FALSE)

