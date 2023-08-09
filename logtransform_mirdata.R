#logtranform the rpms 

setwd("~/tcga-ov-data") #home computer
# Load packages
library("TCGAbiolinks")
library("caret")
library("dplyr")
library("tidyr")
library("tidyverse")
library("edgeR")
library("EDASeq")
library("GDCRNATools")
library("glmnet")
library("gplots")
library("SummarizedExperiment")
library("rockchalk")
#library(multiMiR) #neatidaro

#colums are samples 
tcga_mir_data = readRDS(file = "tcga_mir_data.RDS") 
#filter low expressed mirs
table(rowSums(tcga_mir_data > 0.5))
#at least 2 rpm larger that 0.5
summary(rowSums(tcga_mir_data > 0.5) >= 2) #1338 mirnos all good
#filter low counts
tcga_mir_data.keep <- tcga_mir_data[rowSums(tcga_mir_data > 0.5) >= 2,]
dim(tcga_mir_data.keep)
#normalize
tcga_mir_data.keep <- as.matrix(tcga_mir_data.keep)
logcounts <- as.matrix(log2(tcga_mir_data.keep+0.5))

saveRDS(object = logcounts,
        file = "tcga_mir_LOGdata.RDS",
        compress = FALSE)


# Read mirdata and full clinical data 
tcga_mir_LOGdata = readRDS(file = "tcga_mir_LOGdata.RDS") 
full_clinical <- readRDS(file = "full_clinical.RDS")
#i will first transpose the count data 
tcga_mir_LOGdata_t <- as.data.frame(t(tcga_mir_LOGdata))
dim(tcga_mir_LOGdata_t) #499 cases ir 1881 mirs
#as there are no full barcodes in clinical data need to fix
tcga_mir_LOGdata_t$bcr_patient_barcode <- rownames(tcga_mir_LOGdata_t)
tcga_mir_LOGdata_t$bcr_patient_barcode <- gsub("-01.*","", tcga_mir_LOGdata_t$bcr_patient_barcode)
#these 2A samples should be removed as they are "Recurrent Solid Tumor"
tcga_mir_LOGdata_t <- tcga_mir_LOGdata_t[!grepl("-02A", tcga_mir_LOGdata_t$bcr_patient_barcode),] #491 cases
str(tcga_mir_LOGdata_t$bcr_patient_barcode) #491 names liko
#still 2 non unique zmones yra:
x <- data.frame(table(tcga_mir_LOGdata_t$bcr_patient_barcode))
x <- x[x$Freq == 2, ]
xnames <- c("TCGA-09-0366",
            "TCGA-23-1023")
x <- filter(tcga_mir_LOGdata_t, tcga_mir_LOGdata_t$bcr_patient_barcode %in% xnames)
rownames(x)
#pašalinu nes su 1R žymėm pavadinimas "TCGA-23-1023-01R-01R-1564-13"
#"TCGA-09-0366-01A-01R-1986-13" is in the "removed samples list in firehose:
#https://gdac.broadinstitute.org/runs/tmp/sample_report__2018_01_17/Replicate_Samples.html
x <- c("TCGA-09-0366-01A-01R-1986-13", "TCGA-23-1023-01R-01R-1564-13" )
tcga_mir_LOGdata_t <- tcga_mir_LOGdata_t[!(row.names(tcga_mir_LOGdata_t) %in% x),]
dim(tcga_mir_LOGdata_t) #489 liko
#filter clinical data set (so that the extra clinical data cases are not present)
sum(full_clinical$bcr_patient_barcode %in% tcga_mir_LOGdata_t$bcr_patient_barcode) #bus 489 zmones
mir_data_LOGdf <- semi_join(full_clinical, tcga_mir_LOGdata_t, by='bcr_patient_barcode')
dim(mir_data_LOGdf) #489 cases
#first I wanna add the full barcodes to the phenodata
#best way is to join clinical data (mir_data_df) with count data (tcga_mir_data_t)
#first add back the full barcodes to 
tcga_mir_LOGdata_t$full_barcode <- rownames(tcga_mir_LOGdata_t)
full_LOGdata_mirs <- full_join(tcga_mir_LOGdata_t, mir_data_df, "bcr_patient_barcode" )
dim(full_LOGdata_mirs) # 1340 + 66 = 1405 nes liko tik vienas stulpelis per kuri junge!
colnames(full_LOGdata_mirs[1337:1405]) #see the join at the end
saveRDS(full_LOGdata_mirs, "full_LOGdata_mirs.RDS")

#split off the weird and test sets
#first split weird data by definitions/colnames
colnames(mir_data_df)
table(mir_data_df$prior_treatment, useNA = "a") #noriu nusplitint kaip ir mrna datoj, bet cia 4 su na!
split_prior_treatment <- split(mir_data_df, f = mir_data_df$prior_treatment, drop = T)
main_mir_clinical <- split_prior_treatment$No #nuo siol man lieka tiek cases splitui
weird_cases_mir <- split_prior_treatment$Yes #ir tie 3 na kaip ir prarasti lieka....

# Set (random-number-generator) seed so that results are consistent between runs
set.seed(18) #choose favorite number
train_ids <-rbinom(nrow(main_mir_clinical), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
y_train = main_mir_clinical[train_ids, ] #new train data
y_test  = main_mir_clinical[!train_ids, ] #new test data 
#save 3 files in csv
write.csv(weird_cases_mir, "weird_cases_mir.csv", quote = T, row.names = F) #normal quote F uzdedam bet kadangi kai kur gali but kableliu neveiks
write.csv(y_train, "train_data.csv_mir", quote = T, row.names = F)
write.csv(y_test, "test_data.csv_mir", quote = T, row.names = F )

#spilt the microRNA counts as well 
test_barcodes_mir <- y_test[['bcr_patient_barcode']]
test_LOGcounts_mir <- tcga_mir_LOGdata_t[tcga_mir_LOGdata_t$bcr_patient_barcode %in% test_barcodes_mir, ]
#nepamirst kad countuose yra ir bcr_patient_barcode

train_barcodes_mir <- y_train[['bcr_patient_barcode']]
train_LOGcounts_mir <- tcga_mir_LOGdata_t[tcga_mir_LOGdata_t$bcr_patient_barcode %in% train_barcodes_mir, ]

weird_barcodes_mir <- weird_cases_mir[['bcr_patient_barcode']]
weird_LOGcounts_mir <- tcga_mir_LOGdata_t[tcga_mir_LOGdata_t$bcr_patient_barcode %in% weird_barcodes_mir, ]

#save count data
write.csv(test_LOGcounts_mir, "test_LOGcounts_mir.csv", quote = F, row.names = T)
write.csv(train_LOGcounts_mir, "train_LOGcounts_mir.csv", quote = F, row.names = T)
write.csv(weird_LOGcounts_mir, "weird_LOGcounts_mir.csv", quote = F, row.names = T)

#save normalized counts trains to rds for the ease too 
write_rds(train_LOGcounts_mir, "train_LOGcounts_mir.RDS")
train_LOGcounts_mir = read_rds("train_LOGcounts_mir.RDS") #new train counts
#train_counts_mir <- read.csv("train_counts_mir.csv", header = T)

#split full data (full_data_mirs)
test_LOGfull_mir <- full_LOGdata_mirs[full_LOGdata_mirs$bcr_patient_barcode %in% test_barcodes_mir, ]
train_LOGfull_mir <- full_LOGdata_mirs[full_LOGdata_mirs$bcr_patient_barcode %in% train_barcodes_mir, ]
weird_LOGfull_mir <- full_LOGdata_mirs[full_LOGdata_mirs$bcr_patient_barcode %in% weird_barcodes_mir, ]
saveRDS(test_LOGfull_mir, "test_LOGfull_mir.RDS")
saveRDS(train_LOGfull_mir, "train_LOGfull_mir.RDS")

#already normalized in 1st step
###############################################################################
#model training
# step 1: set up data
train_LOGfull_mir <- read_rds("train_LOGfull_mir.RDS")
rownames(train_LOGfull_mir) <- train_LOGfull_mir$full_barcode
#I need a count matrix  where persons are in rows and genes in columns
#use same train_full_mir 
#I need response factor
train_response <- as.factor(train_LOGfull_mir$vital_status)
# now separate normalized counts
mir_train_LOGnorms <- train_LOGfull_mir[1:1338]
mir_train_LOGnorms <- t(mir_train_LOGnorms) #how is it suddenly numeric? fiiiiine
############model
#parameter alpha: Elastic Net will behave more like LASSO (alpha = 1) or like Ridge Regression (alpha = 0)
# Train model on training dataset using cross-validation
res_mir = cv.glmnet(
  x = mir_train_LOGnorms,
  y = train_response,
  alpha = 0.5, 
  family = "binomial")# 43 seclected with min lambda!
res_coef_mirnos_norm_RPM= coef(res_mir, s="lambda.1se") # the "coef" function returns a sparse matrix
# get coefficients with non-zero values
res_coef_mirnos_norm_RPM = res_coef_mirnos_norm_RPM[res_coef_mirnos_norm_RPM[,1] != 0,] 
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef_mirnos_norm_RPM = res_coef_mirnos_norm_RPM[-1]
relevant_mirs = names(res_coef_mirnos_norm_RPM) # get names of the (non-zero) variables.
length(relevant_mirs) # number of selected genes #43 lambda min, ir #6 genes with 1SE
relevant_mirs #final list

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
mir_train_LOGnorms_stage <- t(mir_train_LOGnorms_stage) #how is it suddenly numeric? fiiiiine
# Train model on training dataset using cross-validation
res_mir_stage = cv.glmnet(
  x = mir_train_LOGnorms_stage,
  y = train_LOGfull_mir_stage$figo_stage_f,
  alpha = 0.5, 
  family = "binomial"
)# 24 seclected with min lambda!
res_coef_mirnos_norm_RPM_stage= coef(res_mir_stage, s="lambda.min") # the "coef" function returns a sparse matrix
# get coefficients with non-zero values
res_coef_mirnos_norm_RPM_stage = res_coef_mirnos_norm_RPM_stage[res_coef_mirnos_norm_RPM_stage[,1] != 0,] 
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef_mirnos_norm_RPM_stage = res_coef_mirnos_norm_RPM_stage[-1]
relevant_mirs = names(res_coef_mirnos_norm_RPM_stage) # get names of the (non-zero) variables.
length(relevant_mirs) # number of selected genes #26 lambda min, ir #0 genes with 1SE
relevant_mirs #final list of 24

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
