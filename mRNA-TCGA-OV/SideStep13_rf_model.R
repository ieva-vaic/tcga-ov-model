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
#partition again
set.seed(18)
ids<- rbinom(nrow(gtex_counts_train), size = 1, prob = 0.8) ==1 #choose persentage
train_gtex <- gtex_counts_train[ids,]
test_gtex <- gtex_counts_train[!ids,]
rf_gtex <- randomForest(group~., data=train_gtex, proximity=TRUE) 
print(rf_gtex)
p1_gtex <- predict(rf_gtex, train_gtex)
confusionMatrix(p1_gtex, train_gtex$group)
p2_gtex <- predict(rf_gtex, test_gtex)
confusionMatrix(p2_gtex, test_gtex$group)
plot(rf_gtex)

hist(treesize(rf_gtex),
     main = "No. of Nodes for the Trees",
     col = "green")
Variable Importance
varImpPlot(rf_gtex,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
top_rf <- as.data.frame(importance(rf_gtex)) 
top_rf <- top_rf %>% filter(MeanDecreaseGini >1)

#makes no sense dar karta partitioninti.
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
top1rf <- rownames(top_rf1)#45
top_rf0<- top_rf %>% filter(MeanDecreaseGini >0)
top0rf <- rownames(top_rf0) #409
#kaip sutampa su pvz elastic net?
elastic <- readRDS("gtcga_elastic.RDS")
venn_lits <- list(elastic_net = elastic, rf = top1rf)
venn(venn_lits)
intersect(elastic, top1rf)
