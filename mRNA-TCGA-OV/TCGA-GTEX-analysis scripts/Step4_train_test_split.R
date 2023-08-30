#Step4 Train/test split
library(tidyverse)
gtcga_counts <- readRDS("C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mrna_voom.RDS")
gtcga_counts <- t(gtcga_counts)
# Set (random-number-generator) seed so that results are consistent between runs
set.seed(18) #choose favorite number
train_ids <- rbinom(nrow(gtcga_counts), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
gtex_counts_train = gtcga_counts[train_ids, ] #new train data 107 35117
gtex_counts_test  = gtcga_counts[!train_ids, ] #new test data 107 35117

#add clinical data 
coldata_gtcga <- readRDS("C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/coldata_gtcga_full.RDS")

gtex_cpheno_train = coldata_gtcga[train_ids, ] #new train data 494  48
gtex_pheno_test  = coldata_gtcga[!train_ids, ] #new test data 107  48

#pasitikrinu kiek kokokiu zmoniu atskyre
table(gtex_cpheno_train$gtex)  
table(gtex_pheno_test$gtex) 
#145 control vs 349 ovarian train sete
#35 control vs 72 ovarian test sete


saveRDS(gtex_counts_train, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/train_gtcga_normcounts.RDS")
saveRDS(gtex_counts_test, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/test_gtcga_normcounts.RDS")

saveRDS(gtex_cpheno_train, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/train_gtcga_coldata.RDS")
saveRDS(gtex_pheno_test, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/test_gtcga_coldata.RDS")