#this is a STEP3: a script for spliting data sets
#setwd("/home/ieva/rprojects/TCGA-OV-data") #home wsl - data folder 
#load libraries

library(tidyverse)
library("glmnet")
library("rockchalk")
# parse data
train_LOGfull_mir <- readRDS(snakemake@input[[1]])

#already normalized in 1st step
###############################################################################
#model training
# step 1: set up data
train_LOGfull_mir <- readRDS("y_train.RDS")

#I need a count matrix  where persons are in rows and genes in columns
#I need response factor
train_response <- as.factor(train_LOGfull_mir$vital_status)
# now separate normalized counts
mir_train_LOGnorms <- train_LOGfull_mir[1:1338]
mir_train_LOGnorms <- as.matrix(mir_train_LOGnorms)
str(mir_train_LOGnorms)#how is it suddenly numeric? fiiiiine
############model
#parameter alpha: Elastic Net will behave more like LASSO (alpha = 1) or like Ridge Regression (alpha = 0)
# Train model on training dataset using cross-validation
res_mir = cv.glmnet(
  x = mir_train_LOGnorms,
  y = train_response,
  alpha = 0.5, 
  family = "binomial")# 43 seclected with min lambda!
res_mir
res_coef_mirnos_norm_RPM= coef(res_mir, s="lambda.min") # the "coef" function returns a sparse matrix
# get coefficients with non-zero values
res_coef_mirnos_norm_RPM = res_coef_mirnos_norm_RPM[res_coef_mirnos_norm_RPM[,1] != 0,] 
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef_mirnos_norm_RPM = res_coef_mirnos_norm_RPM[-1]
relevant_mirs = names(res_coef_mirnos_norm_RPM) # get names of the (non-zero) variables.
length(relevant_mirs) # number of selected genes #43 lambda min, ir #6 genes with 1SE
relevant_mirs #final list
#write(relevant_mirs, xargs$output1)
#saveRDS(relevant_mirs, xargs$output1)
write.csv(relevant_mirs, snakemake@output[[1]], quote = F, row.names = T)
##############################################################################
#try other selection variables (stage maybe)
table(train_LOGfull_mir$figo_stage, useNA = "a") #404 zmones, 3 na
#remove na's
train_LOGfull_mir_stage <- filter(train_LOGfull_mir, !is.na(train_LOGfull_mir$figo_stage)) #401
#remove stage I because only one case is present
train_LOGfull_mir_stage <- filter(train_LOGfull_mir_stage, train_LOGfull_mir_stage$figo_stage != "Stage IC") #400
#make stage II, III and iv binomial
train_LOGfull_mir_stage$figo_stage_f = gsub("[ABC]$", "", train_LOGfull_mir_stage$figo_stage)
train_LOGfull_mir_stage$figo_stage_f <- as.factor(train_LOGfull_mir_stage$figo_stage_f)
train_LOGfull_mir_stage$figo_stage_f <- combineLevels(train_LOGfull_mir_stage$figo_stage_f,
                                                      levs = c("Stage II", "Stage III"), newLabel = c("Stage II & III"))

#dabar pasisalino 3 na atvejai, todel reikia suvienodinti ir pakartoti nuskelima
mir_train_LOGnorms_stage <- train_LOGfull_mir_stage[1:1338]
mir_train_LOGnorms_stage <- as.matrix(mir_train_LOGnorms_stage)
str(mir_train_LOGnorms_stage)
#mir_train_LOGnorms_stage <- t(mir_train_LOGnorms_stage) #how is it suddenly numeric? fiiiiine
# Train model on training dataset using cross-validation
res_mir_stage = cv.glmnet(
  x = mir_train_LOGnorms_stage,
  y = train_LOGfull_mir_stage$figo_stage_f,
  alpha = 0.5, 
  family = "binomial"
)# 24 seclected with min lambda!
res_mir_stage
res_coef_mirnos_norm_RPM_stage= coef(res_mir_stage, s="lambda.min") # the "coef" function returns a sparse matrix
# get coefficients with non-zero values
res_coef_mirnos_norm_RPM_stage = res_coef_mirnos_norm_RPM_stage[res_coef_mirnos_norm_RPM_stage[,1] != 0,] 
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef_mirnos_norm_RPM_stage = res_coef_mirnos_norm_RPM_stage[-1]
relevant_mirs2 = names(res_coef_mirnos_norm_RPM_stage) # get names of the (non-zero) variables.
length(relevant_mirs2) # number of selected genes #26 lambda min, ir #0 genes with 1SE
relevant_mirs2 #final list of 24
#write(relevant_mirs2, xargs$output2)
#saveRDS(relevant_mirs2, xargs$output2)
write.csv(relevant_mirs2, snakemake@output[[2]], quote = F, row.names = T)
###############################################################################
#multinomial: stage 3 groups
train_LOGfull_mir_stage$figo_stage_f3 = gsub("[ABC]$", "", train_LOGfull_mir_stage$figo_stage)
train_LOGfull_mir_stage$figo_stage_f3 <- as.factor(train_LOGfull_mir_stage$figo_stage_f3)
#normalizuoti counts jau nuskelti
res_mir_stage3 = cv.glmnet(
  x = mir_train_LOGnorms_stage,
  y = train_LOGfull_mir_stage$figo_stage_f3,
  alpha = 0.5, 
  family = "multinomial")# 24 seclected with min lambda!
res_coef_mirnos_norm_RPM_stage3= coef(res_mir_stage3, s="lambda.min") # the "coef" function returns a sparse matrix
beta <- Reduce(cbind, res_coef_mirnos_norm_RPM_stage3)
beta <- beta[apply(beta != 0, 1, any),]
colnames(beta) <- names(res_coef_mirnos_norm_RPM_stage3)
beta
relevant_mirs_3figo = rownames(beta) # get names of the (non-zero) variables.
# remove first coefficient as this is the intercept, a variable of the model itself
relevant_mirs_3figo = relevant_mirs_3figo[-1]
length(relevant_mirs_3figo) # number of selected genes #26 lambda min, ir #0 genes with 1SE
relevant_mirs_3figo #final list of 24
#write(relevant_mirs_3figo, xargs$output3)
#saveRDS(relevant_mirs_3figo, xargs$output3)
write.csv(relevant_mirs_3figo, snakemake@output[[3]], quote = F, row.names = T)