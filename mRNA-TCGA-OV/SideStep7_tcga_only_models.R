#STEP7 only TCGA analysis
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

saveRDS(pheno_train, "tcga_pheno_train.RDS")
saveRDS(tcga_counts_train, "tcga_counts_train.RDS")

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
grade_glm #45
grade_coef= coef(grade_glm, s="lambda.min") # the "coef" function returns a sparse matrix
grade_coef = grade_coef[grade_coef[,1] != 0,] 
grade_coef = grade_coef[-1]
relevant_genes_grade= names(grade_coef) # get names of the (non-zero) variables.
relevant_genes_grade 
saveRDS(relevant_genes_grade, "tcga_grade_45.RDS")

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
  alpha = 0.5, 
  family = "gaussian"
)
stage_glm #1
stage_coef= coef(stage_glm, s="lambda.min") # the "coef" function returns a sparse matrix
stage_coef = stage_coef[stage_coef[,1] != 0,] 
stage_coef = stage_coef[-1]
relevant_genes_stage= names(stage_coef) # get names of the (non-zero) variables.
relevant_genes_stage  #1

#Lymphovascular invasion
table(pheno_train$lymphaticinvasion)
lumph <- pheno_train[pheno_train$lymphaticinvasion %in% c("0", "1"), ] #lieka 135  69
lumph_counts <- tcga_counts_train[, (colnames(tcga_counts_train) %in% rownames(lumph))] #lieka 135 zmones
lumph_counts <- data.matrix(t(lumph_counts))
lumph_factor <- lumph$lymphaticinvasion
lumph_glm = cv.glmnet(
  x = lumph_counts,
  y = lumph_factor,
  alpha = 1, 
  family = "binomial"
)
lumph_glm
lumph_coef= coef(lumph_glm, s="lambda.min") # the "coef" function returns a sparse matrix
lumph_coef = lumph_coef[lumph_coef[,1] != 0,] 
lumph_coef = lumph_coef[-1]
relevant_genes_lumph= names(lumph_coef) # get names of the (non-zero) variables.
relevant_genes_lumph 
saveRDS(relevant_genes_lumph, "tcga_lumph_4.RDS")

#Lymphovascular invasion
table(pheno_train$tumorresidualdisease, useNA = "a")
residualdis <- pheno_train %>% drop_na(tumorresidualdisease) #lieka 303 zmones
residualdis_counts <- tcga_counts_train[, (colnames(tcga_counts_train) %in% rownames(residualdis))]
residualdis_counts <- data.matrix(t(residualdis_counts))
resdis_factor <- residualdis$tumorresidualdisease
resdis_factor_bin <- recode(resdis_factor, "1-10 mm" = "residual disease",
                            "11-20 mm" = "residual disease", 
                            ">20 mm" = "residual disease")
resdis_glm = cv.glmnet(
  x = residualdis_counts,
  y = resdis_factor_bin,
  alpha = 1, 
  family = "binomial"
)
resdis_glm 
resdis_glmcoef= coef(resdis_glm, s="lambda.min") # the "coef" function returns a sparse matrix
resdis_glmcoef = resdis_glmcoef[resdis_glmcoef[,1] != 0,] 
resdis_glmcoef = resdis_glmcoef[-1]
relevant_genes_resdis= names(resdis_glmcoef) # get names of the (non-zero) variables.
relevant_genes_resdis 
saveRDS(relevant_genes_resdis, "tcga_resdis_13.RDS")

