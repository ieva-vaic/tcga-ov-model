#test final
setwd("~/rprojects/TCGA-OV-data") #wsl
library(glmnet)
library(tidyverse)
library(yardstick)
library(survivalROC)
library(survival)
library(timeROC)
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
colnames(clin_df_joined2)
surv_counts <- clin_df_joined2[, 10:223]

surv_df <- clin_df_joined2[, c(8, 224, 10:223)]
surv_df <- surv_df %>% dplyr::rename(deceased = decesed2, surv_time = overall_survival) #
saveRDS(surv_df, "surv_df_test.RDS")
#need these for survrock
nobs <- NROW(surv_df)
cutoff <- 365

coxnet.df <- surv_df[, (colnames(surv_df) %in% res_coef_cox_names)]
dim(coxnet.df) #79 10 
time <- surv_df$surv_time
event <- surv_df$deceased

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

###############################################################################
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

###############################################################################
#nne METHOD, test data
rez_list <- apply(coxnet.df, 2, survivalROC, Stime = time, status = event, 
                  predict.time = cutoff, method="NNE", span = 0.25*nobs^(-0.20))
for (i in seq_along(rez_list)) {
  png(paste("figures/plot_test_NNE", i, ".png", sep = ""), width=600, height=500, res=120) # start export
  p =  plot(rez_list[[i]]$FP, rez_list[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=paste( "FP", "\n", "AUC = ", round(rez_list[[i]]$AUC,3)),
            ylab="TP",
            main=paste(names(rez_list)[i],", Method = NNE, Year = 1"), )
  abline(0,1)
  print(p) 
  dev.off()
}

################################################################################
# kaplan mayer plots
#separate each gene into high and low expr
surv_df_test <- readRDS("surv_df_test.RDS")
res_coef_cox_names <- read.csv("res_coef_coxnet_names.csv")
res_coef_cox_names <- res_coef_cox_names$x
coxnet.df_test <- surv_df_test[, (colnames(surv_df_test) %in% res_coef_cox_names)]

median_exprs_test <- apply(coxnet.df_test, 2, median)
median_exprs_test

fit = survfit(Surv(surv_time, deceased) ~ EXO1 > 3.52, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/exo1KM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="EXO1")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ RAD50 > 1.22, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/RAD50KM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="RAD50")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ PPT2 > 3.21, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/PPT2KM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="PPT2")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ LUC7L2 > 5.26, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/LUC7L2KM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="LUC7L2")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ PKP3 > 6.43, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/PKP3KM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="PKP3")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ CDCA5 >5.2, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/CDCA5KM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="CDCA5")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ ZFPL1> 2.82, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/ZFPL1KM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="ZFPL1")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ VPS33B >1.38, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/VPS33BKM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="VPS33B")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ GRB7 >6.87, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/GRB7KM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="GRB7")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ TCEAL4 > 6.10, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/TCEAL4KM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="TCEAL4")
dev.off()

#############################################################################
# TESTING OTHER COX MODELS
###############################################################################
#33 genai, test
# if starting from scratch
res_coef_cox_names33 <- read.csv("res_coef_coxnet_names30.csv")
res_coef_cox_names33 <- res_coef_cox_names33$x
surv_df_test <- readRDS("surv_df_test.RDS")
coxnet.test_df33 <- surv_df_test[, (colnames(surv_df_test) %in% res_coef_cox_names33)] #####
dim(coxnet.test_df33)
#need these for survrock
nobs <- NROW(surv_df_test)
cutoff <- 365
time <- surv_df_test$surv_time
event <- surv_df_test$deceased

rez_test_list33 <- apply(coxnet.test_df33, 2, survivalROC, Stime = time, status = event, predict.time = cutoff, method="KM")

for (i in seq_along(rez_test_list33)) {
  filename <- paste("figures/plot_survroc33_test", i, ".png", sep = "") 
  png(filename, width=600, height=500, res=120) # start export
  p =  plot(rez_test_list33[[i]]$FP, rez_test_list33[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=paste( "FP", "\n", "AUC = ", round(rez_test_list33[[i]]$AUC,3)),
            ylab="TP",
            main=paste(names(rez_test_list33)[i],", TEST, Method = KM, Year = 1"), )
  abline(0,1)
  print(p) 
  dev.off()
}

#17 GENU IS 50, NESUTAMPANCIU SU TOP 33 
diff_genes17 <- c("SMPDL3B",  "POLR3GL" , "NUAK2" ,   "MARVELD2" ,
                  "HMGA1",    "RPS10"  ,  "MAP7"  ,   "BAIAP2L1",
                  "FANCG" ,   "MSANTD7",  "SPTBN2" ,  "PRR11" ,
                  "TRAPPC5",  "EGLN2" ,  "GNL3L"  ,  "SYTL4"  ,  "RPL36A" )
coxnet.test_df17 <- surv_df[, (colnames(surv_df_test) %in% diff_genes17)]
dim(coxnet.test_df17)

rez_test_list17 <- apply(coxnet.test_df17, 2, survivalROC, Stime = time, status = event, predict.time = cutoff, method="KM")

for (i in seq_along(rez_test_list17)) {
  filename <- paste("figures/plot_survroc97_49_test", i, ".png", sep = "") 
  png(filename, width=600, height=500, res=120) # start export
  p =  plot(rez_test_list17[[i]]$FP, rez_test_list17[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=paste( "FP", "\n", "AUC = ", round(rez_test_list17[[i]]$AUC,3)),
            ylab="TP",
            main=paste(names(rez_test_list17)[i],", TEST, Method = KM, Year = 1"), )
  abline(0,1)
  print(p) 
  dev.off()
}

median(surv_df_test$PRR11)
fit = survfit(Surv(surv_time, deceased) ~ PRR11 > 4.78, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/PRR11KM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="PRR11")
dev.off()

median(surv_df_test$BAIAP2L1)
fit = survfit(Surv(surv_time, deceased) ~ BAIAP2L1 > 5.81, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/BAIAP2L1KM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="BAIAP2L1")
dev.off()

median(surv_df_test$POLR3GL)
fit = survfit(Surv(surv_time, deceased) ~ POLR3GL > 2.65, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/POLR3GLKM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="POLR3GL")
dev.off()

#49 GENU Is 97, NESUTAMPANCIU SU TOP 50 
diff_genes49 <- read.csv("diff_genes97.csv")
diff_genes49 <- diff_genes49$x
coxnet.test_df49 <- surv_df_test[, (colnames(surv_df_test) %in% diff_genes49)]
dim(coxnet.test_df49)

rez_test_list49 <- apply(coxnet.test_df49, 2, survivalROC, Stime = time, status = event, predict.time = cutoff, method="KM")
for (i in seq_along(rez_test_list49)) {
  filename <- paste("figures/plot_survroc50_17_test", i, ".png", sep = "") 
  png(filename, width=600, height=500, res=120) # start export
  p =  plot(rez_test_list49[[i]]$FP, rez_test_list49[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=paste( "FP", "\n", "AUC = ", round(rez_test_list49[[i]]$AUC,3)),
            ylab="TP",
            main=paste(names(rez_test_list49)[i],", TEST, Method = KM, Year = 1"), )
  abline(0,1)
  print(p) 
  dev.off()
}
#####################################
#time rock for CI visiems 33 genams
time <- surv_df$surv_time
event <- surv_df$deceased

rez_list_survroc_test <- apply(coxnet.df, 2, timeROC, T = time,
                               delta = event, cause=1,weighting="marginal", #cia gali būti "aalen"
                               times=c(182, 365, 1095, 1825), #half-year, year, 3, 5, year survivals, daugiau neleido, errors
                               iid=TRUE)
rez_list_survroc_test

#############################
#papildomi KM, iš tų 33
median(surv_df_test$SNX21)
fit = survfit(Surv(surv_time, deceased) ~ SNX21 > 3.638, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/SNX21KM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="SNX21")
dev.off()


median(surv_df_test$ARPC1B)
fit = survfit(Surv(surv_time, deceased) ~ ARPC1B >  1.25, data=surv_df_test)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)
png("figures/ARPC1BKM_test.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="ARPC1B")
dev.off()
