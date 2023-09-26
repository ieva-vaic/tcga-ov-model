#random forest
library(randomForest)
library(datasets)
library(caret)
library(tidyverse)
library(venn)
setwd("~/rprojects/TCGA-OV-data") #wsl
#because the normalized hplot was bad this might not be the most appropirate?
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_counts_train <- data.frame(gtex_counts_train)
#
snames = rownames(gtex_counts_train);
group = substr(snames, 1, 4); #Sets up level information for samples.
gtex_counts_train$group = as.factor(group)

rf_gtex <- randomForest(group~., data=gtex_counts_train, proximity=TRUE) 
print(rf_gtex)

#test data
gtex_counts_test <- readRDS("test_gtcga_normcounts_prot.RDS")
snames = rownames(gtex_counts_test);
group = substr(snames, 1, 4); #Sets up level information for samples.
gtex_counts_test$group = as.factor(group)

#test the forest
dim(gtex_counts_train) #489 13682
dim(gtex_counts_test) #106 13682
setdiff(colnames(gtex_counts_train), colnames(gtex_counts_test))
setdiff(colnames(gtex_counts_test), colnames(gtex_counts_train))
colnames(gtex_counts_test) <- gsub(pattern = "-", replacement = ".",x = colnames(gtex_counts_test), fixed = TRUE)
setdiff(colnames(gtex_counts_test), colnames(gtex_counts_train))
#finally test
p1_gtex <- predict(rf_gtex, gtex_counts_test)
confusionMatrix(p1_gtex, gtex_counts_test$group)
plot(rf_gtex)
hist(treesize(rf_gtex),
     main = "No. of Nodes for the Trees",
     col = "green")

varImpPlot(rf_gtex,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
top_rf <- as.data.frame(importance(rf_gtex)) 
top_rf1 <- top_rf %>% filter(MeanDecreaseGini >1)
top1rf <- rownames(top_rf1)#52
top1rf
top_rf0<- top_rf %>% filter(MeanDecreaseGini >0)
top0rf <- rownames(top_rf0) #372
top0rf
#kaip sutampa su pvz elastic net?
elastic <- readRDS("gtcga_elastic.RDS")
venn_lits <- list(elastic_net = elastic, rf = top1rf)
venn(venn_lits)
intersect(elastic, top1rf)
intersect(elastic, top0rf)

saveRDS(rf_gtex, "rf_gtex.RDS")
saveRDS(top1rf, "top1rf.RDS")

partialPlot(rf_gtex, gtex_counts_train, TOMM5, "TCGA")


#coxnet noriu
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
                        "figo_stage")]
# create a new boolean variable that has TRUE for dead patients and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"
# are still alive
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

#get the gene expression for my rf genes (rf0 naudoju kad daugiau butu)
gtex_rf_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% top0rf] 
gtex_rf_counts_train <- as.data.frame(t(gtex_rf_counts_train))
gtex_rf_counts_train <- gtex_rf_counts_train %>%  dplyr::select(starts_with("TCGA")) 
gtex_rf_counts_train <- as.data.frame(t(gtex_rf_counts_train))

#add all genes to clin_df
gtex_rf_counts_train$barcode <- rownames(gtex_rf_counts_train)
clin_df_joined <- left_join(clin_df, gtex_rf_counts_train, by = "barcode")

#cox net?
clin_df_joined2 <- clin_df_joined %>% drop_na(overall_survival) %>% drop_na(deceased) 
clin_df_joined2$decesed2 <- as.integer(as.logical(clin_df_joined2$deceased))
time <- clin_df_joined2$overall_survival
status <- clin_df_joined2$decesed2
y2 <- clin_df_joined2[ ,c(8,7)]
names(y2) <- c("time", "status")
y2 <- as.matrix(y2)
head(y2)
surv_counts <- clin_df_joined2[, 9:380]
library(glmnet)
cox_fitx <- glmnet(surv_counts, y2, family="cox")
cox_fitx
plot(cox_fitx)
coef_x <- coef(cox_fitx, s = 0.05)
head(coef_x)
coef_x = coef_x[coef_x[,1] != 0,] 
res_coef_cox_names = names(coef_x) # get names of the (non-zero) variables.
res_coef_cox_names #51
saveRDS(cox_fitx, "cox_fitx_rf.RDS")
coxnet_elastic <- read_csv("res_coef_coxnet_names.csv")
coxnet_elastic <- coxnet_elastic$x
intersect(coxnet_elastic, res_coef_cox_names)


