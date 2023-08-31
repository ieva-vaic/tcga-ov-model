#Step4 Train/test split
setwd("~/rprojects/TCGA-OV-data") #wsl
library(tidyverse)
gtcga_counts <- readRDS("mrna_voom_no_pseudo.RDS")
gtcga_counts <- t(gtcga_counts)
# Set (random-number-generator) seed so that results are consistent between runs
set.seed(18) #choose favorite number
train_ids <- rbinom(nrow(gtcga_counts), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
gtex_counts_train = gtcga_counts[train_ids, ] #new train data 107 35117 -> 494 10925
gtex_counts_test  = gtcga_counts[!train_ids, ] #new test data 107 35117 -> 107 10925

#add clinical data -> no change to clinical data 
saveRDS(gtex_counts_train, "train_gtcga_normcounts_np.RDS")
saveRDS(gtex_counts_test, "test_gtcga_normcounts_np.RDS")
