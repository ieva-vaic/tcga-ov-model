#Step10 survival:
#esmė bus ta, kad čia pasiliksiu tik TCGA žmones, nes jiems turiu survival
#ir tik 221 geną, kurį atsirenka elatic net tarp sveikų ir ne

setwd("~/rprojects/TCGA-OV-data") #wsl
library(glmnet)
library(tidyverse)
library(survival)
library(gplots)
library(survminer)
library(survivalROC)
library(gridExtra)
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


###############################################################################
#okay su clinical feature modeliais done lets test genes
#for simplicity pasiimsiu 1 tiesiog pradziai ´
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