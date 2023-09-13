#Step10 survival:
#esmė bus ta, kad čia pasiliksiu tik TCGA žmones, nes jiems turiu survival
#ir tik 221 geną, kurį atsirenka elatic net tarp sveikų ir ne

setwd("~/rprojects/TCGA-OV-data") #wsl
library(glmnet)
library(tidyverse)
library("survival")
library(gplots)
library(survminer)
library(RegParallel)
library(survivalROC)
library(gridExtra)
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_genes <- readRDS("gtcga_elastic.RDS")
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 

gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train))
gtex_filtered_counts_train2 <- gtex_filtered_counts_train2 %>%  dplyr::select(starts_with("TCGA")) 
#336 samples #221 genes

#clinical data
pheno_train <- readRDS("tcga_pheno_train.RDS")
#336 samples, 69 variables

# we are only interested in the "Primary solid Tumor" cases for survival
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

#simple survival
fit = survfit(Surv(clin_df$overall_survival, clin_df$deceased) ~0, data=clin_df) #kol kas no model e-g- ~ gender
print(fit)
# we produce a Kaplan Meier plot, no model!
ggsurvplot(fit, data=clin_df)

# jei su stage
# remove any of the letters "a", "b" or "c", but only if they are at the end
# of the name, eg "stage iiia" would become simply "stage iii"
clin_df$tumor_stage = gsub("[ABC]$", "", clin_df$figo_stage)
clin_df[which(clin_df$tumor_stage == "Stage I"), "tumor_stage"] = NA
table(clin_df$tumor_stage, useNA = "a")

fit = survfit(Surv(overall_survival, deceased) ~ tumor_stage, data=clin_df)
pval = surv_pvalue(fit, data=clin_df)$pval
print(pval) #best ever 0.04841577
# we produce a Kaplan-Meier plot from the fitted model
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T)

###############################################################################
#okay su clinical feature modeliais done lets test genes
#for simplicity pasiimsiu 1 tiesiog pradziai ´
gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train2))
clin_df$gene_value = gtex_filtered_counts_train2[, colnames(gtex_filtered_counts_train2) %in% "ARMCX4" ]

# find the median value of the gene and print it
median_value = median(clin_df$gene_value)
print(median_value)
# divide patients in two groups, up and down regulated.
# if the patient expression is greater or equal to them median we put it among the "up-regulated", otherwise among the "down-regulated"
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

# we can fit a survival model, like we did in the previous section
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)

# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=clin_df)$pval
print(pval) #0.5695867 ns
# and finally, we produce a Kaplan-Meier plot
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T, title=paste("ARMCX4"))


#add all genes to clin_df
clin_df <- clin_df[, -11]
gtex_filtered_counts_train2$barcode <- rownames(gtex_filtered_counts_train2)
clin_df_joined <- left_join(clin_df, gtex_filtered_counts_train2, by = "barcode")
colnames(clin_df_joined) <- gsub("[.]", "_", colnames(clin_df_joined)) 
colnames(clin_df_joined) <- gsub("[-]", "_", colnames(clin_df_joined)) 
#po to kai genu vardus imsiu nepamist cia problem del RP11-398C13.8 IR PANASIU SU , IR -
res <- RegParallel(
  data = clin_df_joined,
  formula = 'Surv(overall_survival, deceased) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(clin_df_joined)[11:ncol(clin_df_joined)],
  blocksize = 221,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)
res <- res[!is.na(res$P),]
res_stat_signf <- res %>% arrange(P) %>% filter(P<0.05)
View(res_stat_signf)

write.csv(res_stat_signf$Term, "res_stat_signf1.csv")
##############################################################################
#top6 i viena modeli?
fit_COX = coxph(Surv(overall_survival, deceased) ~ GRB7+PPT2+TPM3+VPS33B+LUC7L2+PKP3, data=clin_df_joined)
fit_COX # p=4.656e-06 TIK TPM3 IR GRB7 PATIKIMI

###############################################################################
#CUTE PLOTS THAT MEAN SOMETHING?
#GRB7
median_value_GRB7 = median(clin_df_joined$GRB7)
print(median_value_GRB7)
clin_df_joined$GRB7_f= ifelse(clin_df_joined$GRB7 >= median_value_GRB7, "UP", "DOWN")
fit_GRB7 = survfit(Surv(overall_survival, deceased) ~ GRB7_f, data=clin_df_joined)
pval_GRB7 = surv_pvalue(fit_GRB7, data=clin_df_joined)$pval
print(pval_GRB7) 
GRB7 <- ggsurvplot(fit_GRB7, data=clin_df_joined, pval=T, risk.table=T, title=paste("GRB7"))
#PPT2
median_value_PPT2 = median(clin_df_joined$PPT2)
print(median_value_PPT2)
clin_df_joined$PPT2_f= ifelse(clin_df_joined$PPT2 >= median_value_PPT2, "UP", "DOWN")
fit_PPT2 = survfit(Surv(overall_survival, deceased) ~ PPT2_f, data=clin_df_joined)
pval_PPT2 = surv_pvalue(fit_PPT2, data=clin_df_joined)$pval
print(pval_PPT2) 
PPT2 <- ggsurvplot(fit_PPT2, data=clin_df_joined, pval=T, risk.table=T, title=paste("PPT2"))
#TPM3
median_value_TPM3 = median(clin_df_joined$TPM3)
print(median_value_TPM3)
clin_df_joined$TPM3_f= ifelse(clin_df_joined$TPM3 >= median_value_TPM3, "UP", "DOWN")
fit_TPM3 = survfit(Surv(overall_survival, deceased) ~ TPM3_f, data=clin_df_joined)
pval_TPM3 = surv_pvalue(fit_TPM3, data=clin_df_joined)$pval
print(pval_TPM3) 
TPM3 <- ggsurvplot(fit_TPM3, data=clin_df_joined, pval=T, risk.table=T, title=paste("TPM3"))
#VPS33B
median_value_VPS33B = median(clin_df_joined$VPS33B)
print(median_value_VPS33B)
clin_df_joined$VPS33B_f= ifelse(clin_df_joined$VPS33B >= median_value_VPS33B, "UP", "DOWN")
fit_VPS33B = survfit(Surv(overall_survival, deceased) ~ VPS33B_f, data=clin_df_joined)
pval_VPS33B = surv_pvalue(fit_VPS33B, data=clin_df_joined)$pval
print(pval_VPS33B) 
VPS33B <- ggsurvplot(fit_VPS33B, data=clin_df_joined, pval=T, risk.table=T, title=paste("VPS33B"))
#LUC7L2
median_value_LUC7L2 = median(clin_df_joined$LUC7L2)
print(median_value_LUC7L2)
clin_df_joined$LUC7L2_f= ifelse(clin_df_joined$LUC7L2 >= median_value_LUC7L2, "UP", "DOWN")
fit_LUC7L2 = survfit(Surv(overall_survival, deceased) ~ LUC7L2_f, data=clin_df_joined)
pval_LUC7L2 = surv_pvalue(fit_LUC7L2, data=clin_df_joined)$pval
print(pval_LUC7L2) 
LUC7L2 <- ggsurvplot(fit_LUC7L2, data=clin_df_joined, pval=T, risk.table=T, title=paste("LUC7L2"))
#PKP3
median_value_PKP3 = median(clin_df_joined$PKP3)
print(median_value_PKP3)
clin_df_joined$PKP3_f= ifelse(clin_df_joined$PKP3 >= median_value_PKP3, "UP", "DOWN")
fit_PKP3 = survfit(Surv(overall_survival, deceased) ~ PKP3_f, data=clin_df_joined)
pval_PKP3 = surv_pvalue(fit_PKP3, data=clin_df_joined)$pval
print(pval_PKP3) 
PKP3 <- ggsurvplot(fit_PKP3, data=clin_df_joined, pval=T, risk.table=T, title=paste("PKP3"))

list_of_plots <- list(GRB7, PPT2, TPM3, VPS33B, LUC7L2, PKP3)
jpeg(file="figures/kaplanstop6.jpeg",
width=20, height=15, units="in", res=100)
arrange_ggsurvplots(list_of_plots, ncol=3, nrow=2 )
dev.off()

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
surv_counts <- clin_df_joined2[, 11:231]


cox_fitx <- glmnet(surv_counts, y2, family="cox", maxit = 1000)
cox_fitx
plot(cox_fitx)
coef_x <- coef(cox_fitx, s = 0.05)
head(coef_x)

coef_x = coef_x[coef_x[,1] != 0,] 
res_coef_cox_names = names(coef_x) # get names of the (non-zero) variables.
res_coef_cox_names #35
#write.csv(res_coef_cox_names, "res_coef_cox_names.csv")

lasoo <- c("TTC4",    "SLC39A1", "TMEM110", "RAD50",   "ANKHD1",  "ZBTB9",
           "RPS10" , "CLDN4",   "PFDN5", "PAGR1" , "RNASEK," ,"GPS2" ,"RTEL1")

intersect(lasoo, res_coef_cox_names)

###############################################################################
#i want a df with colums: time, censor, genes exp
rownames(clin_df_joined2) <- clin_df_joined2$barcode
surv_df <- clin_df_joined2[, c(8, 232, 11:231)]
surv_df <- surv_df %>% rename(censor = decesed2, surv_time = overall_survival)
write.csv(surv_df, "top_glm_for_surv.csv")
