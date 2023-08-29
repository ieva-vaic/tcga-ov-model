#Cleaning, quality control 
library("TCGAbiolinks")
library("caret")
library("dplyr")
library("tidyr")
library("tidyverse")
library("edgeR")
library("EDASeq")
library("GDCRNATools")#neveikia nes reik sudo apt install cmake 
setwd("~/tcga-ov-data")
# train_counts <- read.csv( "train_counts.csv", header = T) #ant readu 
# saveRDS(train_counts, "train_counts_cleaned.RDS")
# #if using count data 
# # as noriu kad butu countuose pats pirmas stulpelis barkodai ir kartu kad griztu i rownamus jei neber
# train_counts[1:10, 1:10]
# train_counts <- rownames(train_counts$barcode)
# train_counts <- train_counts %>% relocate(barcode, .before = ENSG00000000003.15)

###########################################################################################
#Normalization
# for most normal analyses I want a matrix as it is given by assay function
train_counts_matrix <- readRDS("train_count_matrix.RDS") #load matrix #60660 genes
# option 1: pradziai normalizacijai pabandysiu paprasciausiai panaudoti TCGAbiolinks funkcija
#kadangi reik gene info pirma reikai tcga original object uzkraut ir ji istraukt just in case
tcga_data = readRDS(file = "tcga_data.RDS") #pakeiciau pavadinima
#normalization su TCGAanalyze #46357 genes
train_dataNorm <- TCGAanalyze_Normalization(
  tabDF = train_counts_matrix, 
  geneInfo =  geneInfoHT
)
#clean low and high counts
# quantile filter of genes #34766 genes
train_dataNorm_filtered <- TCGAanalyze_Filtering(
  tabDF = train_dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)
#2nd option for norm is limavoom (does not fiflter off any genes)
trainrnaExpr <- gdcVoomNormalization(counts = train_counts_matrix, filter = FALSE) #60660 genes
#3rd option edgeR limma bet jie nori is karto kad duociau modeli kazkoki
##############################################################################
#pasitikrint ar nera double transkriptu, gal genus
#remove the duplicate genes and make gene names the matrix row names
train_dataNorm_filtered <- dataFilt[which(!duplicated(rownames(dataFilt))),] # 34766 genes left, nothings filtered
trainrnaExpr <- trainrnaExpr[which(!duplicated(rownames(trainrnaExpr))),] # 60660 genes left, nothings filtered
#################################################################################
#oulier ieskot 
#kol kas praleidziu nes tutoriale nedare 
###############################################################################
#save final normalized data
write.csv(train_dataNorm_filtered, "train_count_data_cleaned.csv", quote = T) #nuimu rownames false optiona nes man reiks ju po to
saveRDS(trainrnaExpr, "limawoom_tcgaov.RDS")

