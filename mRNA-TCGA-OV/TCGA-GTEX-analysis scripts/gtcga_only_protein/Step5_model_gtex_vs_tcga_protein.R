#Step 5, model gtex vs tcga
setwd("~/rprojects/TCGA-OV-data") #wsl
library(glmnet)
library(biomaRt)
library(tidyverse)
#because the normalized hplot was bad this might not be the most appropirate?
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_counts_train <- data.matrix(gtex_counts_train)
#
snames = rownames(gtex_counts_train);
group = substr(snames, 1, 4); #Sets up level information for samples.
group = as.factor(group)
#Model selection, lasso (no weak values left)

#using norm.counts
#clinical feature: gtex or TCGA data
res_gtex = cv.glmnet(
  x = gtex_counts_train,
  y = group,
  alpha = 1,
  family = "binomial"
)
res_gtex #atrenka 13
# Getting genes that contribute for the prediction
res_coef_gtex = coef(res_gtex, s="lambda.min") # the "coef" function returns a sparse matrix
head(res_coef_gtex) # in a sparse matrix the "." represents the value of zero
# get coefficients with non-zero values
res_coef_gtex = res_coef_gtex[res_coef_gtex[,1] != 0,] 
# note how performing this operation changed the type of the variable
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef_gtex = res_coef_gtex[-1]
res_coef_gtex_names = names(res_coef_gtex) # get names of the (non-zero) variables.
length(res_coef_gtex_names) # number of selected genes #13 lambda min, ir #13 genes with 1SE
# TTC4"    "SLC39A1" "TMEM110" "RAD50"   "ZBTB9"   "NUDT3"   "RPS10"   "CLDN4"   "PFDN5"   "PAGR1"  
# "RNASEK"  "GPS2"    "RTEL1"  
GDC_filt_protein_only <- res_coef_gtex_names
################################################################################
#visu ju palyginimas

#  GDC RNA tools, no filter
GDCnofilter <- c("EVA1B" ,  "COQ8A" ,  "NICN1" ,  "ACAD11",  "CCDC39"  ,"ZBTB9"  , 
                 "TOMM6" ,  "ILK"  ,   "EEF1G"   ,"TAX1BP3",
                 "GABARAP", "WDR83OS" ,"MT-TP")

# GDC RNA tools, with filter
GDCfilter <- c( "TTC4" , "SLC39A1","RP5-1061H20.4","TMEM110", "RAD50","ZBTB9",     
                "NUDT3","RPS10","CLDN4","RNASEK", "GPS2","MT-TP")

# XENA tut, with filter, model added
GTEXdata_TutorialXENA <- c("EVA1B","NICN1","ACAD11","CCDC39","ZBTB9","TOMM6", 
                           "EEF1G","TAX1BP3","WDR83OS","MT-TP" )

#GDC RNA TOOLS, with filter, model added

GDC_filt_model<- c("TTC4", "SLC39A1","RP5-1061H20.4", "TMEM110" , "RAD50", "ZBTB9",        
                   "NUDT3","RPS10","CLDN4","RNASEK" ,"GPS2","MT-TP")

#XENA
XENA <- c("TPX2", "MISP", "FAM83D", "NEK2", "EPCAM", "KLK8","GLUL", "ABCA10", "FOXQ1", "CHMP4C", "TUBA1C",
          "KLK7", "KSR2", "PNLIP", "FAM83H", "RABIF", "TCEAL3", "TMEM185B")
gene_list <- list(GDC_filt_model = GDC_filt_model,GTEXdata_TutorialXENA= GTEXdata_TutorialXENA, 
                  GDCfilter = GDCfilter,GDCnofilter= GDCnofilter, GDC_filt_protein_only= GDC_filt_protein_only)
venn(gene_list)

intersect(GDCfilter,GDC_filt_protein_only )
