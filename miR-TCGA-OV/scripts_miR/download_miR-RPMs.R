#this is a STEP1: a script for downloading TCGA-OV miRNA-seq data in RPMs (reads per million mapped). 
#There is no gene lenght normalization performed as microRNA is allways ~22nt long
setwd("/home/ieva/rprojects/TCGA-OV-data/") #home wsl - data folder 
#load libraries
library("tidyverse")
library("SummarizedExperiment")
library("TCGAbiolinks")
#make TCGA query
query_TCGA_mi = GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "miRNA Expression Quantification",
  experimental.strategy = "miRNA-Seq",
  workflow.type = "BCGSC miRNA Profiling", 
  access = "open" )
#Download your data (do not forget to setwd to data folder, creates a big folder)
GDCdownload(query = query_TCGA_mi, files.per.chunk = 200) 
#Prepare your data
tcga_mir_data <- GDCprepare(query_TCGA_mi)
#chek how much data was preped
dim(tcga_mir_data) #way too many colums per sample
rownames(tcga_mir_data) <- tcga_mir_data$miRNA_ID #this is needed to keep mir names
#count data 
read_real_countData <-  colnames(tcga_mir_data)[grep("read_count", colnames(tcga_mir_data))]
tcga_mir_countdata<- tcga_mir_data[,read_real_countData]
colnames(tcga_mir_countdata) <- gsub("read_count_","", colnames(tcga_mir_countdata)) 
dim(tcga_mir_countdata) #now only 499
tcga_mir_countdata <- as.matrix(tcga_mir_countdata)
str(tcga_mir_countdata)
saveRDS(object = tcga_mir_countdata,
        snakemake@output[[1]],
        compress = FALSE)
########################################        
# I'll be using RPMS's data, then log-transform 
read_countData <-  colnames(tcga_mir_data)[grep("reads", colnames(tcga_mir_data))]
tcga_mir_data <- tcga_mir_data[,read_countData]
colnames(tcga_mir_data) <- gsub("reads_per_million_miRNA_mapped_","", colnames(tcga_mir_data)) 
dim(tcga_mir_data) #now only 499
saveRDS(object = tcga_mir_data,
        snakemake@output[[2]],
        compress = FALSE)
##################################################################################
#This is step 1b. LOG-trasnform/Normalize RPMs, resulting i a smaller countmatrix
#filter low expressed mirs
table(rowSums(tcga_mir_data > 0.5))
#at least 2 rpm larger that 0.5
summary(rowSums(tcga_mir_data > 0.5) >= 2) #1338 miRs
#filter low counts
tcga_mir_data.keep <- tcga_mir_data[rowSums(tcga_mir_data > 0.5) >= 2,]
dim(tcga_mir_data.keep)
#normalize
tcga_mir_data.keep <- as.matrix(tcga_mir_data.keep)
logcounts <- as.matrix(log2(tcga_mir_data.keep+0.5))
#save!
saveRDS(object = logcounts,
        snakemake@output[[3]],
        compress = FALSE)
   
