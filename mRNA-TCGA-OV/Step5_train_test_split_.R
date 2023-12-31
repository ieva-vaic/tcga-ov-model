#Step5 Train/test split
setwd("~/rprojects/TCGA-OV-data") #wsl
library(tidyverse)
gtcga_counts <- readRDS("mrna_voom_protein.RDS")
gtcga_counts <- t(gtcga_counts)
gtcga_counts <- as.data.frame(gtcga_counts)
# Set (random-number-generator) seed so that results are consistent between runs
set.seed(18) #choose favorite number
train_ids <- rbinom(nrow(gtcga_counts), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
gtex_counts_train = gtcga_counts[train_ids, ] 
gtex_counts_test  = gtcga_counts[!train_ids, ] 

dim(gtex_counts_train) #489 13674
dim(gtex_counts_test) #106 13674

snames_train = rownames(gtex_counts_train);
group_train = as.factor(substr(snames_train, 1, 4))
summary(group_train) #153gtex  336tcga

snames_test = rownames(gtex_counts_test);
group_train = as.factor(substr(snames_test, 1, 4))
summary(group_train) #27gtex   79 tcga

#save
saveRDS(gtex_counts_train, "train_gtcga_normcounts_prot.RDS")
saveRDS(gtex_counts_test, "test_gtcga_normcounts_prot.RDS")

#split clinical data as well
pheno <- readRDS("tcga_no_weird_pheno_XENA_TCGA.RDS") #full clinical
dim(pheno) #416  69
#the pheno only has the tcga ids, thus the need to filter
tcga_train_ <- grep("TCGA", snames_train)
tcga_train_ids <- snames_train[tcga_train_]
pheno_train  = pheno[tcga_train_ids, ]  #336 samples
dim(pheno_train)


tcga_test_ <- grep("TCGA", snames_test)
tcga_test_ids <- snames_test[tcga_test_]
pheno_test = pheno[tcga_test_ids, ]  #79 samples
dim(pheno_test)

#save
saveRDS(pheno_train, "tcga_pheno_train.RDS")
saveRDS(pheno_test, "tcga_pheno_test.RDS")