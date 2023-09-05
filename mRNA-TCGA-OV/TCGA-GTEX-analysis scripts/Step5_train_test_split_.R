#Step4 Train/test split
setwd("~/rprojects/TCGA-OV-data") #wsl
library(tidyverse)
gtcga_counts <- readRDS("mrna_voom_protein.RDS")
gtcga_counts <- t(gtcga_counts)
gtcga_counts <- as.data.frame(gtcga_counts)
# Set (random-number-generator) seed so that results are consistent between runs
set.seed(18) #choose favorite number
train_ids <- rbinom(nrow(gtcga_counts), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
gtex_counts_train = gtcga_counts[train_ids, ] #494 13680
gtex_counts_test  = gtcga_counts[!train_ids, ] #107 13680

#add clinical data -> no change to clinical data 
saveRDS(gtex_counts_train, "train_gtcga_normcounts_prot.RDS")
saveRDS(gtex_counts_test, "test_gtcga_normcounts_prot.RDS")
