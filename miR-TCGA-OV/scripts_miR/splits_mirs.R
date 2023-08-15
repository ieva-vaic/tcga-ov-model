#this is a STEP3: a script for spliting data sets
#setwd("/home/ieva/rprojects/TCGA-OV-data") #home wsl - data folder 
#load libraries
library(tidyverse)
#read input
full_LOGdata_mirs <- readRDS(snakemake@input[[1]])

#full_LOGdata_mirs <- readRDS("MERGED.RDS")
colnames(full_LOGdata_mirs[1337:1405])
rownames(full_LOGdata_mirs) <- full_LOGdata_mirs$full_barcode

#split off the weird and test sets
#first split weird data by definitions/colnames
table(full_LOGdata_mirs$prior_treatment, useNA = "a") #noriu nusplitint kaip ir mrna datoj, bet cia 4 su na!
split_prior_treatment <- split(full_LOGdata_mirs, f = full_LOGdata_mirs$prior_treatment, drop = T)
full_LOGdata_mirs <- split_prior_treatment$No #nuo siol man lieka tiek cases splitui
weird_cases_mir <- split_prior_treatment$Yes #ir tie 3 na kaip ir prarasti lieka....

# Set (random-number-generator) seed so that results are consistent between runs
set.seed(18) #choose favorite number
train_ids <-rbinom(nrow(full_LOGdata_mirs), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
y_train = full_LOGdata_mirs[train_ids, ] #new train data
y_test  = full_LOGdata_mirs[!train_ids, ] #new test data 
#save 3 files in RDS
saveRDS(y_train, snakemake@output[[1]])
saveRDS(y_test, snakemake@output[[2]])
saveRDS(weird_cases_mir, snakemake@output[[3]])