#test final
setwd("~/rprojects/TCGA-OV-data") #wsl
library(glmnet)
library(tidyverse)
library(yardstick)
library(survivalROC)
#modeliai
elastic_net_model <- readRDS("elastic_net_model_gtex.RDS")
coxnet_model <- readRDS("coxnet_fit.RDS")
#genai
gtex_genes <- readRDS("gtcga_elastic.RDS")
res_coef_cox_names <- read.csv("res_coef_coxnet_names.csv")
res_coef_cox_names <- res_coef_cox_names$x
#duomenys
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_counts_test <- readRDS("test_gtcga_normcounts_prot.RDS")
gtex_counts_train <- data.matrix(gtex_counts_train)
snames = rownames(gtex_counts_train)
group = substr(snames, 1, 4) #Sets up level information for samples.
group = as.factor(group)

gtex_counts_test <- data.matrix(gtex_counts_test)
# Test/Make prediction on test datasetcfor elastic net
y_pred = predict(elastic_net_model, newx=gtex_counts_test, type="class", s="lambda.min")
length(y_pred) #106
snamest <- rownames(gtex_counts_test)
groupt <-substr(snamest, 1, 4)
y_test <- data.frame(snamest, groupt)
rownames(y_test) <- y_test$snamest
y_test <- as.matrix(y_test)
y_test <- y_test[,-1]
y_test <- as.matrix(y_test)
confusion_matrix = table(y_pred, y_test) 
confusion_matrix #niekas nieko neconfusina.

#CIA KAIP IŠ PAT PRADŽIŲ TESTINOM KAI "BOXES" BUVO ON WEEK 1 IN SWEDEN
print(paste0("Precision: ",precision(confusion_matrix)))
y_pred = predict(elastic_net_model, type="response", newx =gtex_counts_test )
pred_class <- y_pred >0.5
test_class <- y_test == "TCGA"
v <- pred_class == test_class
sum(v) / nrow(y_pred) #loss accuracy?
##########################################################################
pheno <- readRDS("tcga_no_weird_pheno_XENA_TCGA.RDS") #full clinical
#split off the tcga test samples
test_ids <- rownames(gtex_counts_test)
pheno_test  <- pheno[rownames(pheno) %in% test_ids, ]  #106 samples
# susitvarkom survival df for test data
clin_df = pheno_test[,
                   c("barcode",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "neoplasmhistologicgrade",
                     "figo_stage")]
# create a new boolean variable that has TRUE for dead patients and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)
clin_df$tumor_stage = gsub("[ABC]$", "", clin_df$figo_stage)
clin_df[which(clin_df$tumor_stage == "Stage I"), "tumor_stage"] = NA
table(clin_df$tumor_stage, useNA = "a")                                  
colnames(clin_df)

#join with test data
tcga_test_counts <- gtex_counts_test[,(colnames(gtex_counts_test) %in% gtex_genes)]  
tcga_test_counts <- as.data.frame(t(tcga_test_counts))
tcga_test_counts <- tcga_test_counts %>% dplyr::select(starts_with("TCGA")) 
dim(tcga_test_counts) #214 genu ir 79 meginiai 
tcga_test_counts <- as.data.frame(t(tcga_test_counts))
tcga_test_counts$barcode <- rownames(tcga_test_counts)
clin_df_joined <- left_join(clin_df, tcga_test_counts, by = "barcode")
rownames(clin_df_joined) <- clin_df_joined$barcode
head(clin_df_joined)
#be šito coxnetui nepatinka vardai genų, todėl čia darau kad sutaptų kaip reik
colnames(clin_df_joined) <- gsub("[.]", "_", colnames(clin_df_joined)) 
colnames(clin_df_joined) <- gsub("[-]", "_", colnames(clin_df_joined)) 

clin_df_joined2 <- clin_df_joined %>% drop_na(overall_survival) %>% drop_na(deceased) 
clin_df_joined2$decesed2 <- as.integer(as.logical(clin_df_joined2$deceased))
time <- clin_df_joined2$overall_survival
status <- clin_df_joined2$decesed2
y2 <- clin_df_joined2[ ,c(8,7)]
names(y2) <- c("time", "status")
y2 <- as.matrix(y2)
head(y2)
surv_counts <- clin_df_joined2[, 10:223]

surv_df <- clin_df_joined2[, c(8, 224, 10:223)]
surv_df <- surv_df %>% rename(censor = decesed2, surv_time = overall_survival) #
#need these for survrock
nobs <- NROW(surv_df)
cutoff <- 365

coxnet.df <- surv_df[, (colnames(surv_df) %in% res_coef_cox_names)]
dim(coxnet.df) #79 10 
time <- surv_df$surv_time
event <- surv_df$censor

rez_list <- apply(coxnet.df, 2, survivalROC, Stime = time, status = event, predict.time = cutoff, method="KM")
for (i in seq_along(rez_list)) {
png(paste("figures/plot_test", i, ".png", sep = ""), width=600, height=500, res=120) # start export
p =  plot(rez_list[[i]]$FP, rez_list[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab=paste( "FP", "\n", "AUC = ", round(rez_list[[i]]$AUC,3)),
       ylab="TP",
       main=paste(names(rez_list)[i],", Method = KM, Year = 1"), )
  abline(0,1)
print(p) 
dev.off()
}

################################################################
#boxplotas grupes ar atskiria
res_coef_cox_names2 <- c(res_coef_cox_names, "grupe")
gtex_counts_test2 <- as.data.frame(gtex_counts_test)
gtex_counts_test2$grupe <-  rownames(gtex_counts_test)
gtex_counts_test2$grupe <-  substr(gtex_counts_test2$grupe, 1, 4)

png("figures/test_grupes_box.png", width=600, height=500, res=120)
gtex_counts_test2 %>% select(res_coef_cox_names2) %>%
  pivot_longer(., cols = c(res_coef_cox_names[1:10]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = grupe)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle=-90))+
  guides(fill=guide_legend(title="Study"))
dev.off()

##########################################################################
#20231010: veiksmas test data
#1. coxnet mažiau genų = noriu 8
res_coef_cox7_names <- read.csv("res_coef_coxnet_7names.csv")
res_coef_cox7_names <- res_coef_cox7_names$x #"EXO1"   "PPT2"   "LUC7L2" "PKP3"   "CDCA5"  "VPS33B" "GRB7"

coxnet.df7_test <- surv_df[, (colnames(surv_df) %in% res_coef_cox7_names)]
dim(coxnet.df7_test) #79 7
rez_list_7_test <- apply(coxnet.df7_test, 2, survivalROC, Stime = time,
 status = event, predict.time = cutoff, span = 0.25*nobs^(-0.20), method="NNE")
for (i in seq_along(rez_list_7_test)) {
  png(paste("figures/plot_nne7_test_", i, ".png", sep = ""), width=600, height=500, res=120) # start export
p =  plot(rez_list_7_test[[i]]$FP, rez_list_7_test[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab=paste( "FP", "\n", "AUC = ", round(rez_list_7_test[[i]]$AUC,3)),
       ylab="TP",
       main=paste(names(rez_list_7_test)[i],", Method = NNE, Year = 1"), )
  abline(0,1)
print(p) 
dev.off()
}
#2. nne 33 genes 
res_coef_cox33_names <- read.csv("res_coef_coxnet_33names.csv")
res_coef_cox33_names <- res_coef_cox33_names$x 
# "KDF1"     "C1orf210" "MUTYH"    "TTC4"     "DENND2D"  "PSEN2"   
# "EXO1"     "ITPR1"    "STIMATE"  "AMACR"    "RAD50"    "PPT2"    
# "CLDN4"    "ARPC1B"   "LUC7L2"   "IQANK1"   "PKP3"     "ZFPL1"   
# "FAM181A"  "VPS33B"   "RNASEK"   "GRB7"     "ABCA8"    "KLHL14"  
# "C18orf32" "LSR"      "SNX21"    "PARD6B"   "SAMD10"   "IFNAR2"  
# "ADSL"     "GPRASP2"  "TCEAL4"  

coxnet.df33_test <- surv_df[, (colnames(surv_df) %in% res_coef_cox33_names)]
dim(coxnet.df33_test) #79 33
rez_list_33_test <- apply(coxnet.df33_test, 2, survivalROC, Stime = time,
 status = event, predict.time = cutoff, span = 0.25*nobs^(-0.20), method="NNE")
for (i in seq_along(rez_list_33_test)) {
  png(paste("figures/plot_nne33_test_", i, ".png", sep = ""), width=600, height=500, res=120) # start export
p =  plot(rez_list_33_test[[i]]$FP, rez_list_33_test[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab=paste( "FP", "\n", "AUC = ", round(rez_list_33_test[[i]]$AUC,3)),
       ylab="TP",
       main=paste(names(rez_list_33_test)[i],", Method = NNE, Year = 1"), )
  abline(0,1)
print(p) 
dev.off()
}
#4.1. nne 7 genes, cutoff puse metu
cutoff0_5 <- 182
rez_list_7_test_0.5 <- apply(coxnet.df7_test, 2, survivalROC, Stime = time,
 status = event, predict.time = cutoff0_5, span = 0.25*nobs^(-0.20), method="NNE")
for (i in seq_along(rez_list_7_test_0.5)) {
  png(paste("figures/plot_nne7_0.5_", i, ".png", sep = ""), width=600, height=500, res=120) # start export
p =  plot(rez_list_7_test_0.5[[i]]$FP, rez_list_7_test_0.5[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab=paste( "FP", "\n", "AUC = ", round(rez_list_7_test_0.5[[i]]$AUC,3)),
       ylab="TP",
       main=paste(names(rez_list_7_test_0.5)[i],", Method = NNE, Year = 0.5"), )
  abline(0,1)
print(p) 
dev.off()
} 
#4.2. nne 7 genes, cutoff 2 metai
cutoff2 <- 730
rez_list_7_2yr_test <- apply(coxnet.df7_test, 2, survivalROC, Stime = time,
 status = event, predict.time = cutoff2, span = 0.25*nobs^(-0.20), method="NNE")
for (i in seq_along(rez_list_7_2yr_test)) {
  png(paste("figures/plot_nne7_2_test", i, ".png", sep = ""), width=600, height=500, res=120) # start export
p =  plot(rez_list_7_2yr_test[[i]]$FP, rez_list_7_2yr_test[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab=paste( "FP", "\n", "AUC = ", round(rez_list_7_2yr_test[[i]]$AUC,3)),
       ylab="TP",
       main=paste(names(rez_list_7_2yr_test)[i],", Method = NNE, Year = 2"), )
  abline(0,1)
print(p) 
dev.off()
}
#weirdly visi beveik geriau gaunasi in test data
#4.3. nne 33 genes, cutoff puse metu, test
rez_list_33_test <- apply(coxnet.df33_test, 2, survivalROC, Stime = time,
 status = event, predict.time = cutoff0_5, span = 0.25*nobs^(-0.20), method="NNE")
for (i in seq_along(rez_list_33_test)) {
  png(paste("figures/plot_nne33_0.5_test", i, ".png", sep = ""), width=600, height=500, res=120) # start export
p =  plot(rez_list_33_test[[i]]$FP, rez_list_33_test[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
       xlab=paste( "FP", "\n", "AUC = ", round(rez_list_33_test[[i]]$AUC,3)),
       ylab="TP",
       main=paste(names(rez_list_33_test)[i],", Method = NNE, Year = 0.5"), )
  abline(0,1)
print(p) 
dev.off()
}
#5. time rock for CI
#visiems 33 genams
time <- surv_df$surv_time
event <- surv_df$censor

rez_list_survroc33_test <- apply(coxnet.df33_test, 2, timeROC, T = time,
 delta = event, cause=1,weighting="marginal", #cia gali būti "aalen"
times=c(182, 365, 1095, 1825), #half-year, year, 3, 5, year survivals, daugiau neleido, errors
iid=TRUE)
rez_list_survroc33_test

for (i in seq_along(rez_list_survroc33_test)) {
  png(paste("figures/plot_survroc_test", i, ".png", sep = ""), width=600, height=500, res=120) 
  plot(rez_list_survroc33_test[[i]],time = 365 )
dev.off()
}