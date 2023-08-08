#this script is for microRNA data analysis step 2 (after downloading)
#in this scrip I change counts to normalized data
#set wd, libraries
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
#################################################################################
#Save microRNA data
query_TCGA_mi = GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "miRNA Expression Quantification",
  experimental.strategy = "miRNA-Seq",
  workflow.type = "BCGSC miRNA Profiling", 
  access = "open" )

#Download your data (do not forget to setwd to data folder, creates a big folder)
GDCdownload(query = query_TCGA_mi, files.per.chunk = 200) 
#Prepare your data
tcga_mir_data <- GDCprepare(query_TCGA_mi)
dim(tcga_mir_data) #way too many colums per sample
rownames(tcga_mir_data) <- tcga_mir_data$miRNA_ID #this is needed to keep mir names
# I'll be using read_count's data, will need to normalize
read_countData <-  colnames(tcga_mir_data)[grep("reads", colnames(tcga_mir_data))]
tcga_mir_data <- tcga_mir_data[,read_countData]
colnames(tcga_mir_data) <- gsub("reads_per_million_miRNA_mapped_","", colnames(tcga_mir_data)) #now only 499
############################################################################
#clinical data: I think the mir dataset does not have clincial data conveniently atached
#I need to download the clinial stuff seperatly. 
query_TCGA_clinical = GDCquery_clinic(
  project = "TCGA-OV" )
colnames(query_TCGA_clinical)
dim(query_TCGA_clinical) #too many samples
#chek clinical data
rownames(query_TCGA_clinical) <- query_TCGA_clinical$bcr_patient_barcode
#I will narrow down 608 samples later when combining data set with microRNA counts
#Let's save!
saveRDS(object = tcga_mir_data,
        file = "tcga_mir_data.RDS",
        compress = FALSE)

saveRDS(object = query_TCGA_clinical,
        file = "full_clinical.RDS",
        compress = FALSE)
################################################################################
# Read mirdata and full clinical data 
tcga_mir_data = readRDS(file = "tcga_mir_data.RDS") 
full_clinical <- readRDS(file = "full_clinical.RDS")
#i will first transpose the count data 
tcga_mir_data_t <- as.data.frame(t(tcga_mir_data))
dim(tcga_mir_data_t) #499 cases ir 1881 mirs
#as there are no full barcodes in clinical data need to fix
tcga_mir_data_t$bcr_patient_barcode <- rownames(tcga_mir_data_t)
tcga_mir_data_t$bcr_patient_barcode <- gsub("-01.*","", tcga_mir_data_t$bcr_patient_barcode)
#these 2A samples should be removed as they are "Recurrent Solid Tumor"
tcga_mir_data_t <- tcga_mir_data_t[!grepl("-02A", tcga_mir_data_t$bcr_patient_barcode),] #491 cases
str(tcga_mir_data_t$bcr_patient_barcode) #491 names liko
#still 2 non unique zmones yra:
x <- data.frame(table(tcga_mir_data_t$bcr_patient_barcode))
x <- x[x$Freq == 2, ]
xnames <- c("TCGA-09-0366",
            "TCGA-23-1023")
x <- filter(tcga_mir_data_t, tcga_mir_data_t$bcr_patient_barcode %in% xnames)
rownames(x)
#pašalinu nes su 1R žymėm pavadinimas "TCGA-23-1023-01R-01R-1564-13"
#"TCGA-09-0366-01A-01R-1986-13" is in the "removed samples list in firehose:
#https://gdac.broadinstitute.org/runs/tmp/sample_report__2018_01_17/Replicate_Samples.html
x <- c("TCGA-09-0366-01A-01R-1986-13", "TCGA-23-1023-01R-01R-1564-13" )
tcga_mir_data_t <- tcga_mir_data_t[!(row.names(tcga_mir_data_t) %in% x),]
dim(tcga_mir_data_t) #489 liko
#filter clinical data set (so that the extra clinical data cases are not present)
sum(full_clinical$bcr_patient_barcode %in% tcga_mir_data_t$bcr_patient_barcode) #bus 489 zmones
mir_data_df <- semi_join(full_clinical, tcga_mir_data_t, by='bcr_patient_barcode')
dim(mir_data_df) #489 cases

#first I wanna add the full barcodes to the phenodata
#best way is to join clinical data (mir_data_df) with count data (tcga_mir_data_t)
#first add back the full barcodes to 
tcga_mir_data_t$full_barcode <- rownames(tcga_mir_data_t)
full_data_mirs <- full_join(tcga_mir_data_t, mir_data_df, "bcr_patient_barcode" )
dim(full_data_mirs) # 1883 + 66 = 1948 nes liko tik vienas stulpelis per kuri junge!
colnames(full_data_mirs[1880:1948]) #see the join at the end
saveRDS(full_data_mirs, "full_data_mirs.RDS")
#########################################################################
#split off the weird and test sets
#first split weird data by definitions/colnames
colnames(mir_data_df)
table(mir_data_df$prior_treatment, useNA = "a") #noriu nusplitint kaip ir mrna datoj, bet cia 4 su na!
split_prior_treatment <- split(mir_data_df, f = mir_data_df$prior_treatment, drop = T)
main_mir_clinical <- split_prior_treatment$No #nuo siol man lieka tiek cases splitui
weird_cases_mir <- split_prior_treatment$Yes #ir tie 3 na kaip ir prarasti lieka....

###################################################################################
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
test_counts_mir <- tcga_mir_data_t[tcga_mir_data_t$bcr_patient_barcode %in% test_barcodes_mir, ]
#nepamirst kad countuose yra ir bcr_patient_barcode

train_barcodes_mir <- y_train[['bcr_patient_barcode']]
train_counts_mir <- tcga_mir_data_t[tcga_mir_data_t$bcr_patient_barcode %in% train_barcodes_mir, ]

weird_barcodes_mir <- weird_cases_mir[['bcr_patient_barcode']]
weird_counts_mir <- tcga_mir_data_t[tcga_mir_data_t$bcr_patient_barcode %in% weird_barcodes_mir, ]


#save count data
write.csv(test_counts_mir, "test_counts_mir.csv", quote = F, row.names = T)
write.csv(train_counts_mir, "train_counts_mir.csv", quote = F, row.names = T)
write.csv(weird_counts_mir, "weird_counts_mir.csv", quote = F, row.names = T)

#save normalized counts trains to rds for the ease too 
write_rds(train_counts_mir, "train_counts_mir.RDS")
train_counts_mir = read_rds("train_counts_mir.RDS") #new train counts
#train_counts_mir <- read.csv("train_counts_mir.csv", header = T)

#split full data (full_data_mirs)
test_full_mir <- full_data_mirs[full_data_mirs$bcr_patient_barcode %in% test_barcodes_mir, ]
train_full_mir <- full_data_mirs[full_data_mirs$bcr_patient_barcode %in% train_barcodes_mir, ]
weird_full_mir <- full_data_mirs[full_data_mirs$bcr_patient_barcode %in% weird_barcodes_mir, ]
saveRDS(test_full_mir, "test_full_mir.RDS")
saveRDS(train_full_mir, "train_full_mir.RDS")
################################################################################
#no need to normalize this time as there are no raw counts
###############################################################################
#model training
# step 1: set up data
train_full_mir <- read_rds("train_full_mir.RDS")
#I need a count matrix  where persons are in rows and genes in columns
#use same train_full_mir 
#I need response factor
train_response <- as.factor(train_full_mir$vital_status)
# now separate normalized counts
mir_train_norms <- train_full_mir[1:1881]
mir_train_norms <- t(mir_train_norms) #how is it suddenly numeric? fiiiiine
############model
#parameter alpha: Elastic Net will behave more like LASSO (alpha = 1) or like Ridge Regression (alpha = 0)
# Train model on training dataset using cross-validation
res_mir = cv.glmnet(
  x = mir_train_norms,
  y = train_response,
  alpha = 0.5, 
  family = "binomial"
)# 60 seclected with min lambda!
res_coef_mirnos_norm_RPM= coef(res_mir, s="lambda.min") # the "coef" function returns a sparse matrix
# get coefficients with non-zero values
res_coef_mirnos_norm_RPM = res_coef_mirnos_norm_RPM[res_coef_mirnos_norm_RPM[,1] != 0,] 
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef_mirnos_norm_RPM = res_coef_mirnos_norm_RPM[-1]
relevant_mirs = names(res_coef_mirnos_norm_RPM) # get names of the (non-zero) variables.
length(relevant_mirs) # number of selected genes #26 lambda min, ir #0 genes with 1SE
relevant_mirs #final list

# find targets maybe?
#heatmap?
# define the color palette for the plot
hmcol = colorRampPalette(rev(brewer.pal(9, "RdBu")))(256)
# perform complete linkage clustering
clust = function(x) hclust(x, method="complete")
# use the inverse of correlation as distance.
dist = function(x) as.dist((1-cor(t(x)))/2)
# As you've seen a good looking heatmap involves a lot of parameters: does not work
gene_heatmap = heatmap.2(
  t(mir_train_norms[,relevant_mirs]),
  scale="row",          # scale the values for each gene (row)
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  col=hmcol,            # define the color map
  labCol=FALSE,         # Not showing column labels
  ColSideColors=as.character(as.numeric(train_response)), # Show colors for each response class
  dendrogram="both",    # Show dendrograms for both axis
  hclust = clust,       # Define hierarchical clustering method
  distfun = dist,       # Using correlation coefficient for distance function
  cexRow=.6,            # Resize row labels
  margins=c(2,6)        # Define margin spaces
) #red vs black is vital status


#try other selection variables (stage maybe)
table(train_full_mir$figo_stage, useNA = "a") #404 zmones, 3 na
#remove na's
train_full_mir_stage <- filter(train_full_mir, !is.na(train_full_mir$figo_stage)) #401
#remove stage I because only one case is present
train_full_mir_stage <- filter(train_full_mir_stage, train_full_mir_stage$figo_stage != "Stage IC") #400
#make stage II, III and iv binomial
train_full_mir_stage$figo_stage_f = gsub("[ABC]$", "", train_full_mir_stage$figo_stage)
train_full_mir_stage$figo_stage_f <- as.factor(train_full_mir_stage$figo_stage_f)
train_full_mir_stage$figo_stage_f <- combineLevels(train_full_mir_stage$figo_stage_f,
                                                levs = c("Stage II", "Stage III"), newLabel = c("Stage II & III"))

#dabar pasisalino 3 na atvejai, todel reikia suvienodinti ir pakartoti nuskelima
mir_train_norms_stage <- train_full_mir_stage[1:1881]
mir_train_norms_stage <- t(mir_train_norms_stage) #how is it suddenly numeric? fiiiiine
# Train model on training dataset using cross-validation
res_mir_stage = cv.glmnet(
  x = mir_train_norms_stage,
  y = train_full_mir_stage$figo_stage_f,
  alpha = 0.5, 
  family = "binomial"
)# 3 seclected with min lambda!
res_coef_mirnos_norm_RPM_stage= coef(res_mir_stage, s="lambda.min") # the "coef" function returns a sparse matrix
# get coefficients with non-zero values
res_coef_mirnos_norm_RPM_stage = res_coef_mirnos_norm_RPM_stage[res_coef_mirnos_norm_RPM_stage[,1] != 0,] 
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef_mirnos_norm_RPM_stage = res_coef_mirnos_norm_RPM_stage[-1]
relevant_mirs = names(res_coef_mirnos_norm_RPM_stage) # get names of the (non-zero) variables.
length(relevant_mirs) # number of selected genes #26 lambda min, ir #0 genes with 1SE
relevant_mirs #final list of 3

