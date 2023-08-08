#this script is for microRNA data analysis step 2 (after downloading)

#set wd, libraries
setwd("C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data") #home computer
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
library("SummarizedExperiment")
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
read_countData <-  colnames(tcga_mir_data)[grep("count", colnames(tcga_mir_data))]
tcga_mir_data <- tcga_mir_data[,read_countData]
colnames(tcga_mir_data) <- gsub("read_count_","", colnames(tcga_mir_data)) #now only 499
############################################################################
#clinical data: I think the mir dataset does not have clincial data conveniently atached
#I need to download the linial stuff seperatly. 
query_TCGA_clinical = GDCquery_clinic(
  project = "TCGA-OV" )
colnames(query_TCGA_clinical)
dim(query_TCGA_clinical) #too many samples
#chek clinical data
sum(query_TCGA_clinical$bcr_patient_barcode %in% mir_samples) #489 matching
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
tcga_mir_data_t <- tcga_mir_data_t[!grepl("-02A", tcga_mir_data_t$bcr_patient_barcode),]
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
dim(tcga_mir_data_t)
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
colnames(full_data_mirs[1900:1948]) #see the join at the end
saveRDS(full_data_mirs, "full_data_mirs.RDS")
#########################################################################
#split off the weird and test sets
#first split weird data by definitions/colnames
colnames(mir_data_df)
table(mir_data_df$prior_treatment) #noriu nusplitint kaip ir mrna datoj
split_prior_treatment <- split(mir_data_df, f = mir_data_df$prior_treatment, drop = T)
main_mir_clinical <- split_prior_treatment$No
weird_cases_mir <- rbind(weird_cases, split_prior_treatment$Yes) 

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

train_counts_mir <- read.csv("train_counts_mir.csv", header = T)

#split full data (full_data_mirs)
test_full_mir <- full_data_mirs[full_data_mirs$bcr_patient_barcode %in% test_barcodes_mir, ]
train_full_mir <- full_data_mirs[full_data_mirs$bcr_patient_barcode %in% train_barcodes_mir, ]
weird_full_mir <- full_data_mirs[full_data_mirs$bcr_patient_barcode %in% weird_barcodes_mir, ]
saveRDS(test_full_mir, "test_full_mir.RDS")
saveRDS(train_full_mir, "train_full_mir.RDS")
################################################################################
#apply global normalization 
#step one get count matrix with genes in rows and names in colums
train_counts_mir_matrix <- t(train_counts_mir)
colnames(train_counts_mir_matrix) <- train_counts_mir_matrix[1,]
train_counts_mir_matrix <- train_counts_mir_matrix[-1,]
train_counts_mir_matrix <- train_counts_mir_matrix[-1882, ]
train_counts_mir_matrix <- data.frame(train_counts_mir_matrix)

train_counts_mir_matrix = data.frame(lapply(train_counts_mir_matrix, function(x) as.numeric(as.character(x))),
                check.names=F, row.names = rownames(train_counts_mir_matrix))
train_counts_mir_matrix <- as.matrix(train_counts_mir_matrix)
#train_counts_mir_matrix <- lapply(train_counts_mir_matrix,as.numeric) #neveikia man reikia kad butu skai2iais viskas
trainmirnaExpr <- gdcVoomNormalization(counts = train_counts_mir_matrix, filter = FALSE) #1881 mirs
saveRDS(trainmirnaExpr, "trainmirnaExpr.RDS") #save

###########################################################################
#model training
# step 1: set up data
trainmirnaExpr <- readRDS("trainmirnaExpr.RDS")
#I need a count matrix  where persons are in rows and genes in colums
train_mir_mat <- t(trainmirnaExpr)
#I need response factor
train_response <- as.factor(y_train$vital_status)
############model
#parameter alpha: Elastic Net will behave more like LASSO (alpha = 1) or like Ridge Regression (alpha = 0)
# Train model on training dataset using cross-validation
res_mir = cv.glmnet(
  x = train_mir_mat,
  y = train_response,
  alpha = 0.5, 
  family = "binomial"
)# noting is selected 
####################
#try deg for normalization
train_full_mir <- readRDS("train_full_mir.RDS")
rownames(train_full_mir) <- train_full_mir$full_barcode
#i need two filtered count matrixes
table(train_full_mir$vital_status)#count how many should be in the groups
mir_alive <- filter(train_full_mir, train_full_mir$vital_status == "Alive" )
mir_dead <- filter(train_full_mir, train_full_mir$vital_status == "Dead" )
#for this i need numeric matrixes where rows represent genes, so let's transpose
mir_alive <- t(mir_alive)
mir_alive<- mir_alive[1:1881, ] #only counts left
mir_dead <- t(mir_dead)
mir_dead<- mir_dead[1:1881, ] #only counts left
#it's not numeric yet though
mir_alive <- as.data.frame(mir_alive)
mir_alive = data.frame(lapply(mir_alive, function(x) as.numeric(as.character(x))),
                                     check.names=F, row.names = rownames(mir_alive))
mir_alive <- as.matrix(mir_alive)#now it's numeric matrix with rows as genes
mir_dead <- as.data.frame(mir_dead)
mir_dead = data.frame(lapply(mir_dead, function(x) as.numeric(as.character(x))),
                       check.names=F, row.names = rownames(mir_dead))
mir_dead <- as.matrix(mir_dead)#now it's numeric matrix with rows as genes
#the deg from tcgabiolinks
mir_dataDEGs <- TCGAanalyze_DEA(
  mat1 = mir_alive,
  mat2 = mir_dead,
  Cond1type = "alive",
  Cond2type = "dead",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT" #Fit a negative binomial generalized log-linear model 
)  #gives 18 mirs, but not suitable for the model, so back to costallabs tutorial

###############################################################################

res_mir = cv.glmnet(
  x = train_mir_mat,
  y = train_response,
  alpha = 0.5, 
  family = "binomial"
)# noting is selected 



#try other selection variables (stage maybe)