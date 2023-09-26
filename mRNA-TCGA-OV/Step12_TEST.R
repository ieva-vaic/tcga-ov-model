#STEP LAST - test againts test data
#1. ELASTIC NET GTEX vs TCGA
#2. COXNET SURVIVAL!
library(glmnet)
library(yardstick)
#1st WOUPS teks perrunint modeli, biski kitaip everytime
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_counts_test <- readRDS("test_gtcga_normcounts_prot.RDS")
gtex_counts_train <- data.matrix(gtex_counts_train)
gtex_counts_test <- data.matrix(gtex_counts_test)
snames = rownames(gtex_counts_train);
group = substr(snames, 1, 4); #Sets up level information for samples.
group = as.factor(group)
#Model selection, lasso (no weak values left)
#using norm.counts
#clinical feature: gtex or TCGA data
res_gtex = cv.glmnet(
  x = gtex_counts_train,
  y = group,
  alpha = 0.5,
  family = "binomial"
)
res_gtex #atrenka 219 this time
plot(res_gtex)
# Getting genes that contribute for the prediction
res_coef_gtex = coef(res_gtex, s="lambda.min") # the "coef" function returns a sparse matrix
head(res_coef_gtex) # in a sparse matrix the "." represents the value of zero
# get coefficients with non-zero values
res_coef_gtex = res_coef_gtex[res_coef_gtex[,1] != 0,] 
res_coef_gtex = res_coef_gtex[-1]
res_coef_gtex_names = names(res_coef_gtex) # get names of the (non-zero) variables.
res_coef_gtex_names 
saveRDS(res_gtex, "ELASTIC_NET_219_modelis.RDS") #for reproducibility cia 
################################################################################
# Test/Make prediction on test dataset
y_pred = predict(res_gtex, newx=gtex_counts_test, type="class", s="lambda.min")
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
print(paste0("Precision: ",precision(confusion_matrix)))
y_pred = predict(res_gtex, type="response", newx =gtex_counts_test )
pred_class <- y_pred >0.5
test_class <- y_test == "TCGA"
v <- pred_class == test_class
sum(v) / nrow(y_pred) #loss accuracy?
###############################################################################
#coxnet redo for test
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_genes <- readRDS("gtcga_elastic.RDS")
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 

gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train))
gtex_filtered_counts_train2 <- gtex_filtered_counts_train2 %>%  dplyr::select(starts_with("TCGA")) 
dim(gtex_filtered_counts_train2)
gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train2))
#336 samples #221 genes

#clinical data
pheno_train <- readRDS("tcga_pheno_train.RDS")
dim(pheno_train)
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

#add all genes to clin_df
gtex_filtered_counts_train2$barcode <- rownames(gtex_filtered_counts_train2)
clin_df_joined <- left_join(clin_df, gtex_filtered_counts_train2, by = "barcode")
colnames(clin_df_joined) <- gsub("[.]", "_", colnames(clin_df_joined)) 
colnames(clin_df_joined) <- gsub("[-]", "_", colnames(clin_df_joined)) 
#COXNET!
clin_df_joined2 <- clin_df_joined %>% drop_na(overall_survival) %>% drop_na(deceased) 
clin_df_joined2$decesed2 <- as.integer(as.logical(clin_df_joined2$deceased))
time <- clin_df_joined2$overall_survival
status <- clin_df_joined2$decesed2
y2 <- clin_df_joined2[ ,c(8,7)]
names(y2) <- c("time", "status")
y2 <- as.matrix(y2)
head(y2)
surv_counts <- clin_df_joined2[, 9:227]
surv_counts <- as.matrix(surv_counts)
cox_fitx <- glmnet(surv_counts, y2, family="cox")
cox_fitx
plot(cox_fitx)
coef_x <- coef(cox_fitx, s = 0.08)#pasirinkau lambda 1se
head(coef_x)
coef_x = coef_x[coef_x[,1] != 0,] 
res_coef_cox_names = names(coef_x) # get names of the (non-zero) variables.
res_coef_cox_names #12
saveRDS(res_coef_cox_names, "12coxnet.1se.RDS")
saveRDS(cox_fitx, "12coxnet.cox_fitx.RDS")
set.seed(1)
cvfit <- cv.glmnet(surv_counts, y2, family = "cox", type.measure = "C")
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se

