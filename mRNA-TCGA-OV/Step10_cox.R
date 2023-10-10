#Step10 survival:
#esmė bus ta, kad čia pasiliksiu tik TCGA žmones, nes jiems turiu survival
#ir tik 214 geną, kurį atsirenka elatic net tarp sveikų ir ne

setwd("~/rprojects/TCGA-OV-data") #wsl
library(glmnet)
library(tidyverse)
library(survival)
library(gplots)
library(survminer)
library(survivalROC)
library(gridExtra)
library(timeROC)
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_genes <- readRDS("gtcga_elastic.RDS")
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 

gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train))
gtex_filtered_counts_train2 <- gtex_filtered_counts_train2 %>%  dplyr::select(starts_with("TCGA")) 
dim(gtex_filtered_counts_train2)
#336 samples #214 genes

#clinical data
pheno_train <- readRDS("tcga_pheno_train.RDS") #FORM STEP7
#336 samples, 69 variables

# survival df
clin_df = pheno_train[,
                   c("barcode",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "neoplasmhistologicgrade",
                     "figo_stage")]

# create a new boolean variable that has TRUE for dead patients and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)
# remove any of the letters "a", "b" or "c", but only if they are at the end
# of the name, eg "stage iiia" would become simply "stage iii"
clin_df$tumor_stage = gsub("[ABC]$", "", clin_df$figo_stage)
clin_df[which(clin_df$tumor_stage == "Stage I"), "tumor_stage"] = NA
table(clin_df$tumor_stage, useNA = "a")

saveRDS(clin_df, "pheno_survival_only.RDS")
###############################################################################
gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train2))

#add all genes to clin_df
colnames(clin_df)
gtex_filtered_counts_train2$barcode <- rownames(gtex_filtered_counts_train2)
clin_df_joined <- left_join(clin_df, gtex_filtered_counts_train2, by = "barcode")
colnames(clin_df_joined) <- gsub("[.]", "_", colnames(clin_df_joined)) 
colnames(clin_df_joined) <- gsub("[-]", "_", colnames(clin_df_joined)) 
##############################################################################
#COXNET!
clin_df_joined2 <- clin_df_joined %>% drop_na(overall_survival) %>% drop_na(deceased) 
clin_df_joined2$decesed2 <- as.integer(as.logical(clin_df_joined2$deceased))
time <- clin_df_joined2$overall_survival
status <- clin_df_joined2$decesed2
y2 <- clin_df_joined2[ ,c(8,7)]
names(y2) <- c("time", "status")
y2 <- as.matrix(y2)
head(y2)
surv_counts <- clin_df_joined2[, 10:223]

cox_fitx <- glmnet(surv_counts, y2, family="cox")
cox_fitx
png(file="figures/coxnet_cof.png", height=350, width=350)
plot(cox_fitx)
dev.off()

surv_counts <- as.matrix(surv_counts)
set.seed(5)
cvfit <- cv.glmnet(surv_counts, y2, family = "cox", type.measure = "C") #takes some time
png(file="figures/coxnet_cof_cvfit.png", height=350, width=350)
plot(cvfit)
dev.off()
cvfit$lambda.min
cvfit$lambda.1se

coef_x <- coef(cox_fitx, s = 0.088)
head(coef_x)
coef_x = coef_x[coef_x[,1] != 0,] 
res_coef_cox_names = names(coef_x) # get names of the (non-zero) variables.
res_coef_cox_names #10
write.csv(res_coef_cox_names, "res_coef_coxnet_names.csv")

###############################################################################
saveRDS(cox_fitx, "coxnet_fit.RDS")
saveRDS(cvfit, "coxnet_cvfit.RDS")

res_coef_cox_names <- read.csv("res_coef_coxnet_names.csv")
res_coef_cox_names <- res_coef_cox_names$x
cox_fitx <- readRDS("coxnet_fit.RDS")
cvfit <- readRDS("coxnet_cvfit.RDS")
###############################################################################
#test my 10 genes by survroc
#fist I need a matrix with overall survival, cencorship, and genes form elastic net
rownames(clin_df_joined2) <- clin_df_joined2$barcode
surv_df <- clin_df_joined2[, c(8, 224, 10:223)]

surv_df <- surv_df %>% dplyr::rename(censor = decesed2, surv_time = overall_survival) 
write.csv(surv_df, "top_glm_for_surv.csv")
#need these for survrock
nobs <- NROW(surv_df)
cutoff <- 365

coxnet.df <- surv_df[, (colnames(surv_df) %in% res_coef_cox_names)]
dim(coxnet.df) #334  10 #kazkur pametu 2 zmones, kol kas neieskosiu
time <- surv_df$surv_time
event <- surv_df$censor
rez_list <- apply(coxnet.df, 2, survivalROC, Stime = time, status = event, predict.time = cutoff, method="KM")
for (i in seq_along(rez_list)) {
png(paste("figures/plot_", i, ".png", sep = ""), width=600, height=500, res=120) # start export
p =  plot(rez_list[[i]]$FP, rez_list[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab=paste( "FP", "\n", "AUC = ", round(rez_list[[i]]$AUC,3)),
       ylab="TP",
       main=paste(names(rez_list)[i],", Method = KM, Year = 1"), )
  abline(0,1)
print(p) 
dev.off()
}

################################################################################
#20231010: veiksmas train data
#1. coxnet mažiau genų = noriu 8
coef_7 <- coef(cox_fitx, s = 0.097) 
head(coef_7)
coef_7 = coef_7[coef_7[,1] != 0,] 
res_coef_cox7_names = names(coef_7) # get names of the (non-zero) variables.
res_coef_cox7_names #7
write.csv(res_coef_cox7_names, "res_coef_coxnet_7names.csv")
#2. coxnet daugiau genų = noriu ~30
coef_33 <- coef(cox_fitx, s = 0.051) 
head(coef_33)
coef_33 = coef_33[coef_33[,1] != 0,] 
res_coef_cox33_names = names(coef_33) # get names of the (non-zero) variables.
res_coef_cox33_names #33
write.csv(res_coef_cox33_names, "res_coef_coxnet_33names.csv")
#3. nne metodas survroc 
#PAGAL PVZ
nobs <- NROW(surv_df)
cutoff <- 365
## METHOD = NNE,  Nearest Neighbor Estimation (NNE)
#3.1. nne 7 genes 
coxnet.df7 <- surv_df[, (colnames(surv_df) %in% res_coef_cox7_names)]
dim(coxnet.df7)
rez_list_7 <- apply(coxnet.df7, 2, survivalROC, Stime = time,
 status = event, predict.time = cutoff, span = 0.25*nobs^(-0.20), method="NNE")
for (i in seq_along(rez_list_7)) {
  png(paste("figures/plot_nne7_", i, ".png", sep = ""), width=600, height=500, res=120) # start export
p =  plot(rez_list_7[[i]]$FP, rez_list_7[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab=paste( "FP", "\n", "AUC = ", round(rez_list_7[[i]]$AUC,3)),
       ylab="TP",
       main=paste(names(rez_list_7)[i],", Method = NNE, Year = 1"), )
  abline(0,1)
print(p) 
dev.off()
}
#3.2. nne 33 genes 
coxnet.df33 <- surv_df[, (colnames(surv_df) %in% res_coef_cox33_names)]
dim(coxnet.df33)
rez_list_33 <- apply(coxnet.df33, 2, survivalROC, Stime = time,
 status = event, predict.time = cutoff, span = 0.25*nobs^(-0.20), method="NNE")
for (i in seq_along(rez_list_33)) {
  png(paste("figures/plot_nne33_", i, ".png", sep = ""), width=600, height=500, res=120) # start export
p =  plot(rez_list_33[[i]]$FP, rez_list_33[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab=paste( "FP", "\n", "AUC = ", round(rez_list_33[[i]]$AUC,3)),
       ylab="TP",
       main=paste(names(rez_list_33)[i],", Method = NNE, Year = 1"), )
  abline(0,1)
print(p) 
dev.off()
}
#4. pusės metų cutoff, 2 metų cutoff
#4.1. nne 7 genes, cutoff puse metu
coxnet.df7 <- surv_df[, (colnames(surv_df) %in% res_coef_cox7_names)]
dim(coxnet.df7)
cutoff0_5 <- 182
rez_list_7 <- apply(coxnet.df7, 2, survivalROC, Stime = time,
 status = event, predict.time = cutoff0_5, span = 0.25*nobs^(-0.20), method="NNE")
for (i in seq_along(rez_list_7)) {
  png(paste("figures/plot_nne7_0.5_", i, ".png", sep = ""), width=600, height=500, res=120) # start export
p =  plot(rez_list_7[[i]]$FP, rez_list_7[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab=paste( "FP", "\n", "AUC = ", round(rez_list_7[[i]]$AUC,3)),
       ylab="TP",
       main=paste(names(rez_list_7)[i],", Method = NNE, Year = 0.5"), )
  abline(0,1)
print(p) 
dev.off()
}
#4.2. nne 7 genes, cutoff 2 metai
coxnet.df7 <- surv_df[, (colnames(surv_df) %in% res_coef_cox7_names)]
dim(coxnet.df7)
cutoff2 <- 730
rez_list_7 <- apply(coxnet.df7, 2, survivalROC, Stime = time,
 status = event, predict.time = cutoff2, span = 0.25*nobs^(-0.20), method="NNE")
for (i in seq_along(rez_list_7)) {
  png(paste("figures/plot_nne7_2_", i, ".png", sep = ""), width=600, height=500, res=120) # start export
p =  plot(rez_list_7[[i]]$FP, rez_list_7[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab=paste( "FP", "\n", "AUC = ", round(rez_list_7[[i]]$AUC,3)),
       ylab="TP",
       main=paste(names(rez_list_7)[i],", Method = NNE, Year = 2"), )
  abline(0,1)
print(p) 
dev.off()
}
#4.3. nne 33 genes, cutoff puse metu
#1. nne 33 genes 
coxnet.df33 <- surv_df[, (colnames(surv_df) %in% res_coef_cox33_names)]
dim(coxnet.df33)
rez_list_33 <- apply(coxnet.df33, 2, survivalROC, Stime = time,
 status = event, predict.time = cutoff0_5, span = 0.25*nobs^(-0.20), method="NNE")
for (i in seq_along(rez_list_33)) {
  png(paste("figures/plot_nne33_0.5_", i, ".png", sep = ""), width=600, height=500, res=120) # start export
p =  plot(rez_list_33[[i]]$FP, rez_list_33[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab=paste( "FP", "\n", "AUC = ", round(rez_list_33[[i]]$AUC,3)),
       ylab="TP",
       main=paste(names(rez_list_33)[i],", Method = NNE, Year = 0.5"), )
  abline(0,1)
print(p) 
dev.off()
}
#5. time rock for CI
#pirmiausia vienam genui pabandyt
timeROC_EXO1 <- timeROC(
T=time,
delta=event,marker=coxnet.df33$EXO1,
cause=1,weighting="marginal", #cia gali būti "aalen"
times=c(182, 365, 730, 1095, 1460, 1825, 3650), #half-year, year, 3, 5, 10 year survivals
iid=TRUE)
timeROC_EXO1 #estimates
plot(timeROC_EXO1,time = 365 )
#visiems 33 genams
rez_list_survroc33 <- apply(coxnet.df33, 2, timeROC, T = time,
 delta = event, cause=1,weighting="marginal", #cia gali būti "aalen"
times=c(182, 365, 730, 1095, 1460, 1825, 3650), #half-year, year, 3, 5, 10 year survivals
iid=TRUE)
rez_list_survroc33

for (i in seq_along(rez_list_survroc33)) {
  png(paste("figures/plot_survroc", i, ".png", sep = ""), width=600, height=500, res=120) 
  plot(rez_list_survroc33[[i]],time = 365 )
dev.off()
}