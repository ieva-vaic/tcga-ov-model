#Step 5, model gtex vs tcga
setwd("~/rprojects/TCGA-OV-data") #wsl
library(glmnet)
library(biomaRt)
library(tidyverse)
#because the normalized hplot was bad this might not be the most appropirate?
gtex_counts_train <- readRDS("train_gtcga_normcounts_np.RDS")
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
res_gtex #atrenka 9
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
#"EVA1B"   "COQ8A"   "NICN1"   "ACAD11"  "CCDC39"  "ZBTB9"   "TOMM6"   "ILK"     "EEF1G"   "TAX1BP3"
# "GABARAP" "WDR83OS" "MT-TP" 
bad_apjungimas <- c("EVA1B","ARL6IP4","TSPYL2","CLDN4", "MT-TP","UTP14C","RPL17","UPK3BL1","TMEM35B")
intersect(res_coef_gtex_names, bad_apjungimas)
################################################################################
