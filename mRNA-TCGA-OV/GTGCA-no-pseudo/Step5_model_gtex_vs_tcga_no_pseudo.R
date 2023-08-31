#Step 5, model gtex vs tcga
setwd("~/rprojects/TCGA-OV-data") #wsl
library(glmnet)
library(biomaRt)
library(tidyverse)
#because the normalized hplot was bad this might not be the most appropirate?
gtex_counts_train <- readRDS( "train_gtcga_normcounts_np.RDS")
gtex_pheno_train <- readRDS( "train_gtcga_coldata.RDS")

#Model selection, lasso (no weak values left)

#using norm.counts
#clinical feature: gtex or TCGA data
dim(gtex_counts_train) #buvo tiek:494 35117 ->494 10925
train_response <- as.factor(gtex_pheno_train$gtex)
res_gtex = cv.glmnet(
  x = gtex_counts_train,
  y = train_response,
  alpha = 1,
  family = "binomial"
)
res_gtex #atrenka 9
# Getting genes that contribute for the prediction
res_coef_gtex = coef(res_gtex, s="lambda.min") # the "coef" function returns a sparse matrix
head(res_coef_gtex) # in a sparse matrix the "." represents the value of zero
# get coefficients with non-zero values
res_coef_gtex = res_coef_gtex[res_coef_gtex[,1] != 0,] 
# note how performing this operation changed the type of the variable
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef_gtex = res_coef_gtex[-1]
res_coef_gtex_names = names(res_coef_gtex) # get names of the (non-zero) variables.
length(res_coef_gtex_names) # number of selected genes #13 lambda min, ir #13 genes with 1SE

################################################################################
#nebereikia vardu kankintis perdarinet
# sutampa "EVA1B"   "ARL6IP4" "TSPYL2"  "CLDN4"   "MT-TP"    "UTP14C"  "RPL17"   "UPK3BL1"
# nesutampa TMEM35B - Transmembrane Protein 35B - zinomas smegenu navikuose
# ensembl_gene_id external_gene_name: neprartynus vardu
# 1  ENSG00000142694              EVA1B Eva-1 homolog B, maybe immune function
# 2  ENSG00000182196            ARL6IP4 ADP Ribosylation Factor Like GTPase 6 Interacting Protein 4, maybe a splicing regulator
# 3  ENSG00000184205             TSPYL2 Testis-specific Y-encoded-like protein -> in damage repair, G1 chekpoint
# 4  ENSG00000189143              CLDN4 claudinas4 yra tight junction narys
# 6  ENSG00000210196              MT-TP mitochondrial transfer RNA (tRNA) proline
# 9  ENSG00000253797             UTP14C a ribosome processome component
# 11 ENSG00000265681              RPL17 from the large ribosomal subunit
# 13 ENSG00000267368            UPK3BL1 Uroplakin -urothelium-specific
#pasisalino visi neaiskus pseudogenai realiai, kurie buvo is visu genu darant.
# 5  ENSG00000207165            SNORA70 RNU70, small nucleolar RNA, H/ACA box 70
# 7  ENSG00000226232           NPIPB14P #Nuclear Pore Complex Interacting Protein Family Member B14 Pseudogene.
# 8  ENSG00000231503             PTMAP4 #pseudogene that encodes a prothymosin alpha pseudogene 4.
# 10 ENSG00000256393            RPL41P5 #ribosomal protein L41 pseudogene 5 
# 12 ENSG00000266820            KPNA2P3 #pseaudogene, unknown function