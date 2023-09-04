#XENA TCGA
# as esu pastebejusi kad parsisiunciant data is XENA bazes yra tvarkingiau ir daugiau duomenu
#noriu palyginti su TCGA siaip
#padaryti glm su grade
#padaryti glm su stage
library(tidyverse)
library(WGCNA)
v <- readRDS("v.RDS")
v_counts <- v$E
v_counts <- as.data.frame(v_counts)

#get only TCGA expression data
dataExpr0 = v_counts
#rownames(dataExpr0) <- dataExpr0$X
#head(dataExpr0[1:6]); #contains gtex
dataExpr1 = dataExpr0 %>% select(starts_with("TCGA"));
head(dataExpr1[1:6]); #only tcga left



#Add a new variable "Count" that counts the number of samples that have normalized expression value <0 for gene x
dataExpr1$Count = rowSums(dataExpr1 < 0);
table(dataExpr1$Count);

#Keep genes that have "Count" = 0 (i.e., Keep genes that have non-negative normalized expression value across ALL samples).
dataExpr2 = dataExpr1 %>% filter(Count == 0);
#Remove the "Count" variable from the gene expression matrix.
dataExpr2$Count = NULL
dataExpr3 = as.data.frame(t(dataExpr2))

#praleidziu clusterinima viskas plius minus kartu
dataExpr = dataExpr3


###############################################################################
#prep clinical, cia atrinkti tik clinical yra
dataTraitALL = read.csv("00_ClinTraits.csv");
head(dataTraitALL);

dataTraitALL$sampleID = gsub("-", ".", dataTraitALL$sampleID);
head(dataTraitALL)
#Form a trait data frame analogous to the expression data frame.
exprRows = rownames(dataExpr);
traitRows = match(exprRows, dataTraitALL$sampleID);
dataTrait = dataTraitALL[traitRows, -1];
head(dataTrait);

rownames(dataTrait) = dataTrait$sampleID;
head(dataTrait);

dataTrait = select(dataTrait, -sampleID);
head(dataTrait)

############################################################################
## split
xena_counts <- dataExpr3
xena_counts <-as.data.frame(xena_counts)
set.seed(18) #choose favorite number
train_ids <- rbinom(nrow(xena_counts), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
xena_counts_train = xena_counts[train_ids, ] #new train data 346 2710
xena_counts_test  = xena_counts[!train_ids, ] #new test data 72 2710


####################################################################
#glm for grade II vs III
table(dataTrait$neoplasmhistologicgrade, useNA = "a") #reik padaryt g2 vs g3

#grade
grade <- dataTrait
grade <- filter(grade, !is.na(grade$neoplasmhistologicgrade))#nebutina turbut
grade <- grade[grade$neoplasmhistologicgrade %in% c("G2", "G3"), ]
table(grade$neoplasmhistologicgrade, useNA = "a") #406 lieka 47 G2 +359 G3

#susivienodinam 
grade_samples <- rownames(grade)
xena_counts_train_grade <- filter(xena_counts_train, rownames(xena_counts_train) %in% grade_samples)
dim(xena_counts_train_grade) #xena_counts_train_grade
grade_train <- filter(grade, rownames(grade) %in% rownames(xena_counts_train_grade))
dim(grade_train) #335   7

#prep for glm
xena_counts_train_grade <- data.matrix(xena_counts_train_grade)
grade_factor <- as.factor(grade_train$neoplasmhistologicgrade)


res_factor = cv.glmnet(
  x = xena_counts_train_grade,
  y = grade_factor,
  alpha = 1, 
  family = "binomial"
)
res_factor # 18


# Getting genes that contribute for the prediction
res_coef= coef(res_factor, s="lambda.min") # the "coef" function returns a sparse matrix
head(res_coef) # in a sparse matrix the "." represents the value of zero
# get coefficients with non-zero values
res_coef = res_coef[res_coef[,1] != 0,] 
# note how performing this operation changed the type of the variable
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef = res_coef[-1]
relevant_genes_xena_grade= names(res_coef) # get names of the (non-zero) variables.
relevant_genes_xena_grade 

xena_grade <- c("MTMR11", "POLQ","SDCCAG8","XRCC1","KLHL20", "PES1", "TTLL1" ,
"TMEM30A","PPP2CA","DESI2", "DYM","PTPRN2","SGCB","TMEM182","STARD5","SLC39A10",
"LGR4","TSTD3")  
XENA_gtex_tcga <- c("TPX2", "MISP", "FAM83D", "NEK2", "EPCAM", "KLK8","GLUL",
                    "ABCA10", "FOXQ1", "CHMP4C", "TUBA1C",
          "KLK7", "KSR2", "PNLIP", "FAM83H", "RABIF", "TCEAL3", "TMEM185B")
intersect(xena_grade, XENA_gtex_tcga) #nesutampa niekas

############################################################################
#glm for STAGE II vs III
table(dataTrait$clinicalstage2, useNA = "a") #reik nusalint stage I, na
#stage
stage <- dataTrait
stage <- filter(stage, !is.na(stage$clinicalstage2))#nebutina turbut
stage <- stage[stage$clinicalstage2 %in% c("Stage II", "Stage III"), ]
table(stage$clinicalstage2, useNA = "a") #451 lieka 24 G2 +327 G3

#susivienodinam 
stage_samples <- rownames(stage)
xena_counts_train_stage <- filter(xena_counts_train, rownames(xena_counts_train) %in% stage_samples)
dim(xena_counts_train_stage) #xena_counts_train_grade
stage_train <- filter(stage, rownames(stage) %in% rownames(xena_counts_train_stage))
dim(stage_train) #290   7

#prep for glm
xena_counts_train_stage <- data.matrix(xena_counts_train_stage)
stage_factor <- as.factor(stage_train$clinicalstage2)

res_factor_Stage = cv.glmnet(
  x = xena_counts_train_stage,
  y = stage_factor,
  alpha = 1, 
  family = "binomial"
)
res_factor_Stage #18

# Getting genes that contribute for the prediction
res_coef2= coef(res_factor_Stage, s="lambda.min") # the "coef" function returns a sparse matrix
head(res_coef2) # in a sparse matrix the "." represents the value of zero
# get coefficients with non-zero values
res_coef2 = res_coef2[res_coef2[,1] != 0,] 
# note how performing this operation changed the type of the variable
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef2 = res_coef2[-1]
relevant_genes_xena_stage= names(res_coef2) # get names of the (non-zero) variables.
relevant_genes_xena_stage 

intersect(xena_grade, relevant_genes_xena_stage) #no intersect

XENA_stage <- c( "CAMKK1","ZDHHC6","ZNF76","MYL6","SRRD","MYL9","PSMD3","STAT1","PGD","INTS4" , 
 "PDCD4","PDK1", "CRTAP", "RARG", "ZNF852","METTL7A" ,"TMEM17", "GDAP2","UBE2J1","C1QTNF5")

tcga_stage <- c("MIR4660", "NDUFB8P1", "HSPA8P5", "LINC02384", "SPIN2P1", "SPIN2P1", 
                "TMEM97P1", "EIF3JP1", "IGHV7-81", "FAM180A", "BACE1", "SPPL2C",
                "CCL13", "MTLN", "SUGCT", "KLK15", "FRZB", "TMEM123", 
                "MPZL2", "CREG1", "MMP7", "PRRG1", "SLC44A2", "NDUFA8", "CCL7", 
                "CNFN", "DMAC2", "TRIB3", "SHD", "TRIB3", "NRDC", "REEP1")
intersect(XENA_stage, tcga_stage)

############################################################################
#stage gausian
#susivienodinam -> beda kad yr na

train_traits <- dataTrait
train_traits <- filter(train_traits, !is.na(train_traits$clinicalstage_num))#nebutina turbut
train_traits <- train_traits[train_traits$clinicalstage2 %in% c("Stage II", "Stage III", "Stage IV"), ]
table(train_traits$clinicalstage2, useNA = "a") #451 lieka 24 G2 +327 G3

xena_counts_train_stage_num <- xena_counts_train[rownames(xena_counts_train) %in% rownames(train_traits), ]
dim(xena_counts_train_stage_num)#  343 2710
xena_counts_train_stage_num <- data.matrix(xena_counts_train_stage_num)
stage_num <- train_traits$clinicalstage_num
res_Stage = cv.glmnet(
  x = xena_counts_train_stage_num,
  y = stage_num,
  alpha = 1, 
  family = "gaussian"
)
res_Stage #116

# Getting genes that contribute for the prediction
res_coef3= coef(res_Stage, s="lambda.min") # the "coef" function returns a sparse matrix
res_coef3 = res_coef3[res_coef3[,1] != 0,] 
res_coef3 = res_coef3[-1]
relevant_genes_xena_stage_num= names(res_coef3) # get names of the (non-zero) variables.
relevant_genes_xena_stage_num 

intersect(relevant_genes_xena_stage_num, relevant_genes_xena_stage)
intersect(relevant_genes_xena_stage_num, tcga_stage)
