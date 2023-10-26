#Step2.3 add biomart -> gene names and types
library(tidyverse)
library(biomaRt)
setwd("~/rprojects/TCGA-OV-data")
tcga_data <- readRDS("tcga_no_weird_counts.RDS")
#read from step 2.2
mRNA_counts <- t(tcga_data)
##assay data converted to dataframe
mRNA_counts <- as.data.frame(mRNA_counts) 
dim(mRNA_counts)#60660 transcripts (rows), 416 zmones (cols) 
mRNA_counts[,1:416] <- lapply(mRNA_counts[,1:416], as.numeric)
str(mRNA_counts) #60660 genes

#get gene names in the count df
mRNA_counts$ensembl <- rownames(mRNA_counts)

mRNA_counts$ensembl_gene_id <- gsub("\\..*", "",mRNA_counts$ensembl)
length(mRNA_counts$ensembl_gene_id) # 60660 

#biomart
listEnsembl() #shows available databases
ensembl <- useEnsembl(biomart = "genes") #unsuportive gali but
datasets <- listDatasets(ensembl) #issirinkti human genes
ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl') #name of data base ir data set
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

tcga_genes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', "gene_biotype"), #values to retreve
                     filters = "ensembl_gene_id", #input on quary
                     values = mRNA_counts$ensembl_gene_id,
                     mart = ensembl.con) #sukurtas conection object
dim(tcga_genes) #60419 = 241 pamesta

#prisijoinsiu dfs
counts_tcga_with_gene_names <- left_join(mRNA_counts, tcga_genes, by= "ensembl_gene_id")
rownames(counts_tcga_with_gene_names) <- counts_tcga_with_gene_names$ensembl
dim(counts_tcga_with_gene_names) #60660   420
saveRDS(counts_tcga_with_gene_names, "tcga_with_names_all.RDS")

