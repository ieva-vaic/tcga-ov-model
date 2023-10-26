#Step10 survival:
#ONLY TCGA, as only they have the survival data
#only 214 genes form elastic net

setwd("~/rprojects/TCGA-OV-data") 
library(glmnet)
library(tidyverse)
library(survival)
library(gplots)
library(survminer)
library(survivalROC)
library(data.table)
library(timeROC)
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_genes <- readRDS("gtcga_elastic.RDS")
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 

gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train))
gtex_filtered_counts_train2 <- gtex_filtered_counts_train2 %>%  dplyr::select(starts_with("TCGA")) 
dim(gtex_filtered_counts_train2)
#336 samples #214 genes

#clinical data
pheno_train <- readRDS("tcga_pheno_train.RDS") 
#336 samples, 69 variables

# survival df
clin_df = pheno_train[,
                   c("barcode",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "neoplasmhistologicgrade",
                     "STAGE")]

# create a new boolean variable that has TRUE for dead patients and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)
# remove 1 patient with STAGE I
clin_df[which(clin_df$STAGE == "Stage I"), "STAGE"] = NA
table(clin_df$STAGE, useNA = "a")

saveRDS(clin_df, "pheno_survival_only.RDS")
###############################################################################
#add all genes to clin_df
gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train2))
colnames(clin_df)
gtex_filtered_counts_train2$barcode <- rownames(gtex_filtered_counts_train2)
clin_df_joined <- left_join(clin_df, gtex_filtered_counts_train2, by = "barcode")
#gene names cannot have weird symbols for coxnet
colnames(clin_df_joined) <- gsub("[.]", "_", colnames(clin_df_joined)) 
colnames(clin_df_joined) <- gsub("[-]", "_", colnames(clin_df_joined)) 
##############################################################################
#COXNET!
clin_df_joined2 <- clin_df_joined %>% drop_na(overall_survival) %>% drop_na(deceased) #lost 2 cases
clin_df_joined2$decesed2 <- as.integer(as.logical(clin_df_joined2$deceased))
time <- clin_df_joined2$overall_survival
status <- clin_df_joined2$decesed2
y2 <- clin_df_joined2[ ,c(8,7)]
names(y2) <- c("time", "status")
y2 <- as.matrix(y2)
head(y2)
surv_counts <- clin_df_joined2[, 9:222] #buvo 10:223?

# cox_fit - visas lambdas ismeta
cox_fitx <- glmnet(surv_counts, y2, family="cox")
cox_fitx
png(file="figures/coxnet_cof.png", height=350, width=350)
plot(cox_fitx)
dev.off()

#cv.glmnet - min ir max + paveikslas - takes more time
surv_counts <- as.matrix(surv_counts)
set.seed(5)
cvfit <- cv.glmnet(surv_counts, y2, family = "cox", type.measure = "C") #takes some time
png(file="figures/coxnet_cof_cvfit.png", height=350, width=350)
plot(cvfit)
dev.off()
cvfit #lambda min and max

coef_x <- coef(cox_fitx, s = 0.081) # best 10
head(coef_x)
coef_x = coef_x[coef_x[,1] != 0,] 
res_coef_cox_names = names(coef_x) # get names of the (non-zero) variables.
res_coef_cox_names #10
write.csv(res_coef_cox_names, "res_coef_coxnet_names.csv")


###############################################################################
#save it using later:
saveRDS(cox_fitx, "coxnet_fit.RDS")
saveRDS(cvfit, "coxnet_cvfit.RDS")

#if starting with saved model
res_coef_cox_names <- read.csv("res_coef_coxnet_names.csv")
res_coef_cox_names <- res_coef_cox_names$x
cox_fitx <- readRDS("coxnet_fit.RDS")
cvfit <- readRDS("coxnet_cvfit.RDS")
###############################################################################
#test my 10 genes by survroc
#fist I need a matrix with overall survival, cencorship, and genes form elastic net
rownames(clin_df_joined2) <- clin_df_joined2$barcode
colnames(clin_df_joined2)
surv_df <- clin_df_joined2[, c(8, 223, 9:222)] #buvo 8, 224, 10:223

surv_df <- surv_df %>% dplyr::rename(deceased = decesed2, surv_time = overall_survival) 
write_rds(surv_df, "top_glm_for_surv.RDS")

#need these for survrock
nobs <- NROW(surv_df)
cutoff <- 365

coxnet.df <- surv_df[, (colnames(surv_df) %in% res_coef_cox_names)]
dim(coxnet.df) #334  10 

time <- surv_df$surv_time
event <- surv_df$deceased


rez_list <- apply(coxnet.df, 2, survivalROC, Stime = time, status = event, predict.time = cutoff, method="KM")

for (i in seq_along(rez_list)) {
  filename <- paste("figures/plot_", i, ".png", sep = "") 
  png(filename, width=600, height=500, res=120) # start export
  p =  plot(rez_list[[i]]$FP, rez_list[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=paste( "FP", "\n", "AUC = ", round(rez_list[[i]]$AUC,3)),
            ylab="TP",
            main=paste(names(rez_list)[i],", Method = KM, Year = 1"), )
  abline(0,1)
  print(p) 
  dev.off()
}
###############################################################################
#other survrocs
## METHOD = NNE,  Nearest Neighbor Estimation (NNE)
surv_df <- read
coxnet.df10 <- surv_df[, (colnames(surv_df) %in% res_coef_cox_names)]
dim(coxnet.df10)
rez_list_10 <- apply(coxnet.df10, 2, survivalROC, Stime = time,
                    status = event, predict.time = cutoff, span = 0.25*nobs^(-0.20), method="NNE")
for (i in seq_along(rez_list_10)) {
  png(paste("figures/plot_nne10_", i, ".png", sep = ""), width=600, height=500, res=120) # start export
  p =  plot(rez_list_10[[i]]$FP, rez_list_10[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=paste( "FP", "\n", "AUC = ", round(rez_list_10[[i]]$AUC,3)),
            ylab="TP",
            main=paste(names(rez_list_10)[i],", Method = NNE, Year = 1"), )
  abline(0,1)
  print(p) 
  dev.off()
}
## METHOD =  KM, 0,5 year
cutoff_half_year <- 182
rez_list_half_year <- apply(coxnet.df, 2, survivalROC, Stime = time,
                            status = event, predict.time = cutoff_half_year, method="KM")

for (i in seq_along(rez_list_half_year)) {
  filename <- paste("figures/plot_half_yr", i, ".png", sep = "") 
  png(filename, width=600, height=500, res=120) # start export
  p =  plot(rez_list_half_year[[i]]$FP, rez_list_half_year[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=paste( "FP", "\n", "AUC = ", round(rez_list_half_year[[i]]$AUC,3)),
            ylab="TP",
            main=paste(names(rez_list_half_year)[i],", Method = KM, Year = 0,5"), )
  abline(0,1)
  print(p) 
  dev.off()
}
## METHOD =  KM, 5 year
cutoff_5_year <- 1825
rez_list_5_year <- apply(coxnet.df, 2, survivalROC, Stime = time, 
                         status = event, predict.time = cutoff_5_year, method="KM")

for (i in seq_along(rez_list_5_year)) {
  filename <- paste("figures/plot_5_yr", i, ".png", sep = "") 
  png(filename, width=600, height=500, res=120) # start export
  p =  plot(rez_list_5_year[[i]]$FP, rez_list_5_year[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=paste( "FP", "\n", "AUC = ", round(rez_list_5_year[[i]]$AUC,3)),
            ylab="TP",
            main=paste(names(rez_list_5_year)[i],", Method = KM, Year = 0,5"), )
  abline(0,1)
  print(p) 
  dev.off()
}
## METHOD =  KM, 3 year
cutoff_3_year <- 1095
rez_list_3_year <- apply(coxnet.df, 2, survivalROC, Stime = time, 
                         status = event, predict.time = cutoff_3_year, method="KM")

for (i in seq_along(rez_list_3_year)) {
  filename <- paste("figures/plot_3_yr", i, ".png", sep = "") 
  png(filename, width=600, height=500, res=120) # start export
  p =  plot(rez_list_3_year[[i]]$FP, rez_list_3_year[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=paste( "FP", "\n", "AUC = ", round(rez_list_3_year[[i]]$AUC,3)),
            ylab="TP",
            main=paste(names(rez_list_3_year)[i],", Method = KM, Year = 0,5"), )
  abline(0,1)
  print(p) 
  dev.off()
}
###############################################################################
# time rock for CI
#pirmiausia vienam genui pabandyt
timeROC_EXO1 <- timeROC(
  T=time,
  delta=event,marker=coxnet.df$EXO1,
  cause=1,weighting="marginal", #cia gali būti "aalen"
  times=c(182, 365, 1095, 1825, 3650), #half-year, year, 3, 5, 10 year survivals
  iid=TRUE)
timeROC_EXO1 #estimates
plot(timeROC_EXO1,time = 365 )
#visiems genams timeroc
rez_list_survroc<- apply(coxnet.df, 2, timeROC, T = time,
                            delta = event, cause=1,weighting="marginal", #cia gali būti "aalen"
                            times=c(182, 365, 1095, 1825, 3650), #half-year, year, 3, 5, 10 year survivals
                            iid=TRUE)
rez_list_survroc

for (i in seq_along(rez_list_survroc)) {
  png(paste("figures/plot_survroc", i, ".png", sep = ""), width=600, height=500, res=120) 
  plot(rez_list_survroc[[i]],time = 365 )
  dev.off()
}

rez_table <- data.frame(rez_list_survroc)
#################################################################################
# kaplan mayer plots
#separate each gene into high and low expr
surv_df <- readRDS("top_glm_for_surv.RDS")
res_coef_cox_names <- read.csv("res_coef_coxnet_names.csv")
res_coef_cox_names <- res_coef_cox_names$x
coxnet.df <- surv_df[, (colnames(surv_df) %in% res_coef_cox_names)]
median_exprs <- apply(coxnet.df, 2, median)
median_exprs

#fit survival model: all 10 genes: i could not make all the apply and loops work
res_coef_cox_names
#"EXO1"   "RAD50"  "PPT2"   "LUC7L2" "PKP3"   "CDCA5"  "ZFPL1"  "VPS33B" "GRB7"   "TCEAL4"

fit = survfit(Surv(surv_time, deceased) ~ EXO1 > 3.45, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/exo1KM.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="EXO1")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ RAD50 > 1.28, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/RAD50KM.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="RAD50")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ PPT2 > 3.17, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/PPT2KM.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="PPT2")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ LUC7L2 > 5.21, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/LUC7L2KM.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="LUC7L2")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ PKP3 > 6.5, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/PKP3KM.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="PKP3")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ CDCA5 >5.2, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/CDCA5KM.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="CDCA5")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ ZFPL1> 2.88, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/ZFPL1KM.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="ZFPL1")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ VPS33B >1.23, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/VPS33BKM.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="VPS33B")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ GRB7 >6.78, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/GRB7KM.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="GRB7")
dev.off()

fit = survfit(Surv(surv_time, deceased) ~ TCEAL4 > 6.03, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/TCEAL4KM.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="TCEAL4")
dev.off()



###########################################################################
#2023-10-17 su julium
ggplot(surv_df, aes(log(surv_time), CDCA5, color = deceased)) +
  geom_point() +
  geom_smooth()

#IS TEST FAILO
ggplot(surv_df_test, aes(log(surv_time), CDCA5, color = deceased)) +
  geom_point() +
  geom_smooth() 

mod <- coxph(Surv(surv_time, deceased) ~ CDCA5, data = surv_df_test )
summary(mod)  

mod <- coxph(Surv(surv_time, deceased) ~ RAD50, data = surv_df_test )
summary(mod)  


fit = survfit(Surv(surv_time, deceased) ~ CDCA5 >5, data=surv_df_test)
summary(fit)
pval = surv_pvalue(fit, data=surv_df_test)$pval
print(pval)

ggsurvplot(fit, data=surv_df_test, pval=T, risk.table=T, title="CDCA5")

##############################################################################
#OTHER COX MODELS
#best pagal 1se: 7 genai
coef_x7 <- coef(cox_fitx, s = 0.089) # best 7
head(coef_x7)
coef_x7 = coef_x7[coef_x7[,1] != 0,] 
res_coef_cox_names7 = names(coef_x7) # get names of the (non-zero) variables.
res_coef_cox_names7 #10
write.csv(res_coef_cox_names7, "res_coef_coxnet_names7.csv")

#norint daugiau pasirinkimo: 30 genai
coef_x30 <- coef(cox_fitx, s = 0.051) # best 30
head(coef_x30)
coef_x30 = coef_x30[coef_x30[,1] != 0,] 
res_coef_cox_names30 = names(coef_x30) # get names of the (non-zero) variables.
res_coef_cox_names30 #10
write.csv(res_coef_cox_names30, "res_coef_coxnet_names30.csv")

#kaip sutampa
intersect(res_coef_cox_names30, res_coef_cox_names7)
intersect(res_coef_cox_names30, res_coef_cox_names)
intersect(res_coef_cox_names7, res_coef_cox_names)

#norint daugiau pasirinkimo: 50 genai
coef_x50 <- coef(cox_fitx, s = 0.039) # best 50
head(coef_x50)
coef_x50 = coef_x50[coef_x50[,1] != 0,] 
res_coef_cox_names50 = names(coef_x50) # get names of the (non-zero) variables.
res_coef_cox_names50 #10
write.csv(res_coef_cox_names50, "res_coef_coxnet_names50.csv")
intersect(res_coef_cox_names50, res_coef_cox_names)

diff_genes50 <- setdiff(res_coef_cox_names50, res_coef_cox_names30)
write.csv(diff_genes50, "diff_genes50.csv")

#norint daugiau pasirinkimo: 97 genai (lambda min)
coef_xMIN <- coef(cox_fitx, s = 0.01824) # best 50
head(coef_xMIN)
coef_xMIN = coef_xMIN[coef_xMIN[,1] != 0,] 
res_coef_cox_names97 = names(coef_xMIN) # get names of the (non-zero) variables.
res_coef_cox_names97
write.csv(res_coef_cox_names97, "res_coef_coxnet_names97.csv")

diff_genes97 <- setdiff(res_coef_cox_names97, res_coef_cox_names50)
diff_genes97
write.csv(diff_genes97, "diff_genes97.csv")

#KM for other top genes
#### OTHER TOP GENES
median(surv_df$SNX21)
fit = survfit(Surv(surv_time, deceased) ~ SNX21 > 3.76, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/SNX21KM_tRAIN.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="SNX21")
dev.off()

median(surv_df$ARPC1B)
fit = survfit(Surv(surv_time, deceased) ~ ARPC1B >  1.34, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/ARPC1BKM_tRAIN.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="ARPC1B")
dev.off()

#top 50
median(surv_df$PRR11)
fit = survfit(Surv(surv_time, deceased) ~ PRR11 > 4.78, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/PRR11_tRAIN.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="PRR11")
dev.off()

median(surv_df$BAIAP2L1)
fit = survfit(Surv(surv_time, deceased) ~ BAIAP2L1 > 5.84, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/BAIAP2L1tRAIN.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="BAIAP2L1")
dev.off()

median(surv_df$POLR3GL)
fit = survfit(Surv(surv_time, deceased) ~ POLR3GL >  2.68, data=surv_df)
pval = surv_pvalue(fit, data=surv_df)$pval
print(pval)
png("figures/POLR3GLTRAIN.png", width=600, height=800, res=120)
ggsurvplot(fit, data=surv_df, pval=T, risk.table=T, title="POLR3GL")
dev.off()
#################################################################################
#SURVIVAL ROCs other genes
##33 genai
#if starting from scratch
res_coef_cox_names33 <- read.csv("res_coef_coxnet_names30.csv")
res_coef_cox_names33 <- res_coef_cox_names33$x
surv_df <- read.csv("top_glm_for_surv.csv")
rownames(surv_df) <- surv_df$X
#atsiskiriam coxnet 33 genus
coxnet.df33 <- surv_df[, (colnames(surv_df) %in% res_coef_cox_names33)]
dim(coxnet.df33)
#survroc for 33 genes
#need these for survrock
nobs <- NROW(surv_df)
cutoff <- 365

time <- surv_df$surv_time
event <- surv_df$deceased


rez_list33 <- apply(coxnet.df33, 2, survivalROC, Stime = time, status = event, predict.time = cutoff, method="KM")

for (i in seq_along(rez_list33)) {
  filename <- paste("figures/plot_survroc33", i, ".png", sep = "") 
  png(filename, width=600, height=500, res=120) # start export
  p =  plot(rez_list33[[i]]$FP, rez_list33[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=paste( "FP", "\n", "AUC = ", round(rez_list33[[i]]$AUC,3)),
            ylab="TP",
            main=paste(names(rez_list33)[i],", Method = KM, Year = 1"), )
  abline(0,1)
  print(p) 
  dev.off()
}

#########################################################################
#50 genu, bet skiria nuo 33
#atsiskiriam coxnet 33 genus
coxnet.df50 <- surv_df[, (colnames(surv_df) %in% diff_genes50)]

rez_list50 <- apply(coxnet.df50, 2, survivalROC, Stime = time, status = event, predict.time = cutoff, method="KM")
for (i in seq_along(rez_list50)) {
  filename <- paste("figures/plot_survroc_top97_diff", i, ".png", sep = "") 
  png(filename, width=600, height=500, res=120) # start export
  p =  plot(rez_list50[[i]]$FP, rez_list50[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=paste( "FP", "\n", "AUC = ", round(rez_list50[[i]]$AUC,3)),
            ylab="TP",
            main=paste(names(rez_list50)[i],", Method = KM, Year = 1"), )
  abline(0,1)
  print(p) 
  dev.off()
}

#97 genu, bet skirias nuo 50
#atsiskiriam coxnet 33 genus
coxnet.df97 <- surv_df[, (colnames(surv_df) %in% diff_genes97)]

rez_list_97 <- apply(coxnet.df97, 2, survivalROC, Stime = time, status = event, predict.time = cutoff, method="KM")
for (i in seq_along(rez_list_97)) {
  filename <- paste("figures/plot_survroc_top50_diff", i, ".png", sep = "") 
  png(filename, width=600, height=500, res=120) # start export
  p =  plot(rez_list_97[[i]]$FP, rez_list_97[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=paste( "FP", "\n", "AUC = ", round(rez_list_97[[i]]$AUC,3)),
            ylab="TP",
            main=paste(names(rez_list_97)[i],", Method = KM, Year = 1"), )
  abline(0,1)
  print(p) 
  dev.off()
}

