#Step 5, model gtex vs tcga
setwd("~/rprojects/TCGA-OV-data") #wsl
library(glmnet)
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
res_coef_gtex = res_coef_gtex[-1]
res_coef_gtex_names = names(res_coef_gtex) # get names of the (non-zero) variables.
res_coef_gtex_names 
# "TTC4"    "SLC39A1" "TMEM110" "RAD50"   "ANKHD1"  "ZBTB9"   "RPS10"   "CLDN4"   "PFDN5"   "PAGR1"  
# "RNASEK"  "GPS2"    "RTEL1" 
################################################################################
