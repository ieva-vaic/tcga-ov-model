#STEP7 only TCGA analysis:
## grade
## stage
## dead-alive (vital status)

setwd("~/rprojects/TCGA-OV-data") #wsl
library(glmnet)
library(tidyverse)
#get counts dataframe
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_counts_train <- as.data.frame(t(gtex_counts_train))
#get only TCGA counts
tcga_counts_train = gtex_counts_train %>% select(starts_with("TCGA")) #336 samples
#get clinical
pheno <- readRDS("tcga_no_weird_pheno_XENA_TCGA.RDS") #full clinical
#split off the tcga train samples
train_ids <- colnames(tcga_counts_train)
pheno_train  = pheno[train_ids, ]  #336 samples

#clean-up
rm(pheno)
rm(gtex_counts_train)
rm(train_ids)

###grade
table(pheno_train$neoplasmhistologicgrade, useNA="a") #lyginisu G2 vs G3
grade <- pheno_train
grade <- grade[grade$neoplasmhistologicgrade %in% c("G2", "G3"), ]
table(grade$neoplasmhistologicgrade, useNA = "a") #bus 36 G2 vs 291 G3

grade_counts <- tcga_counts_train[, (colnames(tcga_counts_train) %in% rownames(grade))]
grade_counts <- data.matrix(t(grade_counts))
grade_factor <- grade$neoplasmhistologicgrade

grade_glm = cv.glmnet(
  x = grade_counts,
  y = grade_factor,
  alpha = 1, 
  family = "binomial"
)
grade_glm # 37
grade_coef= coef(grade_glm, s="lambda.min") # the "coef" function returns a sparse matrix
grade_coef = grade_coef[grade_coef[,1] != 0,] 
grade_coef = grade_coef[-1]
relevant_genes_grade= names(grade_coef) # get names of the (non-zero) variables.
relevant_genes_grade 

###stage
table(pheno_train$clinicalstage_num, useNA="a") 
#will remove -2 (STAGE1), ir NA
stage <- pheno_train[pheno_train$clinicalstage_num %in% c("-1", "0", "1"), ]
stage_counts <- tcga_counts_train[, (colnames(tcga_counts_train) %in% rownames(stage))]
stage_counts <- data.matrix(t(stage_counts))
grade_factor <- stage$clinicalstage_num

stage_glm = cv.glmnet(
  x = stage_counts,
  y = grade_factor,
  alpha = 1, 
  family = "gaussian"
)
stage_glm #13
stage_coef= coef(stage_glm, s="lambda.min") # the "coef" function returns a sparse matrix
stage_coef = stage_coef[stage_coef[,1] != 0,] 
stage_coef = stage_coef[-1]
relevant_genes_stage= names(stage_coef) # get names of the (non-zero) variables.
relevant_genes_stage 

intersect(relevant_genes_grade, relevant_genes_stage) #0 witch makes sense

###vital_status
table(pheno_train$vital_status, useNA="a") #lyginisu 208 dead 127 alive 
tcga_counts <- data.matrix(t(tcga_counts_train))
vital_factor <- pheno_train$vital_status

vital_glm = cv.glmnet(
  x = tcga_counts,
  y = vital_factor,
  alpha = 1, 
  family = "binomial"
)
vital_glm # 1? really?
vital_coef= coef(vital_glm, s="lambda.min") # the "coef" function returns a sparse matrix
vital_coef = vital_coef[vital_coef[,1] != 0,] 
vital_coef = vital_coef[-1]
relevant_genes_vital= names(vital_coef) # get names of the (non-zero) variables.
relevant_genes_vital #PPL
