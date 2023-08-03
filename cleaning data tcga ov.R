#valymas, quality control
setwd("~/tcga-ov-data")
library("TCGAbiolinks")
library("caret")
library("dplyr")
library("tidyr")
library("tidyverse")
library("EDASeq")
#library("GDCRNATools")#neveikia nes reik sudo apt install cmake 
library(edgeR)
train_counts <- read.csv( "train_counts.csv", header = T) #ant readu 
saveRDS(train_counts, "train_counts_cleaned.RDS")
# as noriu kad butu countuose pats pirmas stulpelis barkodai ir kartu kad griztu i rownamus jei neber
train_counts[1:10, 1:10]
train_counts <- rownames(train_counts$barcode)
train_counts <- train_counts %>% relocate(barcode, .before = ENSG00000000003.15)
#################################################################################
#oulier ieskot


###########################################################################################
#normalizavima
# for most normal analyses I want a matrix as it is given by assay function
#option 1 go back and make filtered dataframe as a matrix: not doing
#option 2: go back and filter matrix
train_counts_matrix <- readRDS("train_count_matrix.RDS")

# option 1: pradziai normalizacijai pabandysiu paprasciausiai panaudoti TCGAbiolinks funkcija
#kadangi reik gene info pirma reikai tcga original object uzkraut ir ji istraukt
tcga_data = readRDS(file = "tcga_data.RDS") #pakeiciau pavadinima
#normalization su TCGAanalyze
train_dataNorm <- TCGAanalyze_Normalization(
  tabDF = train_counts_matrix, 
  geneInfo =  geneInfoHT
)

#2nd option for norm is limavoom
#trainrnaExpr <- gdcVoomNormalization(counts = train_counts_matrix, filter = FALSE)

#3rd option edgeR limma bet jie nori is karto kad duociau modeli kazkoki

######################
#pasitikrint ar nera double transkriptu, gal genus
#remove the duplicate genes and make gene names the matrix row names
train_data_filtered <- dataFilt[which(!duplicated(rownames(dataFilt))),] # nieko nenufiltravo


#pasalint low counts
# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(
  tabDF = train_dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)
train_clinical <- read.csv( "train_data.csv", header = T) #ant readu 
write.csv(train_clinical, "train_data_cleaned.csv", quote = T, row.names = F)


