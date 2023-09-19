#Step 11 side projects kinda

setwd("~/rprojects/TCGA-OV-data") #wsl
library(tidyverse)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler) 
library(survivalROC)
library(gridExtra)
library(grid)
library(ggpubr)
#try ggo on the 219 elastic net genes
elastic <- readRDS("gtcga_elastic.RDS")
elastic.entrez <- mapIds(org.Hs.eg.db, elastic, 'ENTREZID', 'SYMBOL')
elastic.entrez <- elastic.entrez[!is.na(elastic.entrez)]
ggo.elastic.entrez <- groupGO(gene = elastic.entrez,
                            OrgDb    = org.Hs.eg.db,
                            ont      = "CC",
                            level    = 3,
                            readable = TRUE)
ggo.elastic.entrez <- as.data.frame(ggo.elastic.entrez)
ggo.elastic.entrez.filt <- ggo.elastic.entrez %>% filter(Count != "0")
ggo.entez <- ggo.elastic.entrez.filt
ggo.entez <- ggo.entez %>% arrange(desc(Count))
View(ggo.entez)
ego.recorded.num <- enrichGO(gene = elastic.entrez,
                             OrgDb         = org.Hs.eg.db,
                             ont           = "CC",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = TRUE)
View(as.data.frame(ego.recorded.num))
goplot(ego.recorded.num)
#try go at the coxnet ~35 genes selected for nice pictures and stuff
res_coef_cox_names <- read.csv("res_coef_coxnet_names.csv")
res_coef_cox_names <- res_coef_cox_names$x
  
recorded.num.entrez <- mapIds(org.Hs.eg.db, res_coef_cox_names, 'ENTREZID', 'SYMBOL')
recorded.num.entrez <- recorded.num.entrez[!is.na(recorded.num.entrez)]
ggo.recorded.num <- groupGO(gene = recorded.num.entrez,
                            OrgDb    = org.Hs.eg.db,
                            ont      = "CC",
                            level    = 3,
                            readable = TRUE)
ggo.recorded_num_df <- as.data.frame(ggo.recorded.num)
ggo.filtered.recorded_num_df <- ggo.recorded_num_df %>% filter(Count != "0")
ggo <- ggo.filtered.recorded_num_df
ggo <- ggo %>% arrange(desc(Count))
write.csv(ggo, "go_result_after_survival.csv")
View(ggo)
#ego = enrichment go
ego.recorded.num <- enrichGO(gene = recorded.num.entrez,
                             OrgDb         = org.Hs.eg.db,
                             ont           = "CC",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = TRUE)
View(as.data.frame(ego.recorded.num)) #0 terms

#############################################################################
#just top genes related with survival no coxnet 
res_stat_sign_names = read_csv("res_stat_signf1.csv") 
res_stat_sign_names <- res_stat_sign_names$x

survival.num.entrez <- mapIds(org.Hs.eg.db, res_stat_sign_names, 'ENTREZID', 'SYMBOL')
ggo.survival <- groupGO(gene = survival.num.entrez,
                            OrgDb    = org.Hs.eg.db,
                            ont      = "CC",
                            level    = 3,
                            readable = TRUE)
ggo.survival <- as.data.frame(ggo.survival)
ggo.survival <- ggo.survival %>% filter(Count != "0")
ggo.survival <- ggo.survival %>% arrange(desc(Count))


##############################################################################
#rock but with survival
surv_df <- read_csv("top_glm_for_surv.csv")
surv_df <- as.data.frame(surv_df)
rownames(surv_df) <- surv_df$...1
surv_df <- surv_df[, -1]

#PAGAL PVZ
nobs <- NROW(surv_df)
cutoff <- 365
## METHOD = NNE,  Nearest Neighbor Estimation (NNE) 
SURVROC.rez= survivalROC(Stime=surv_df$surv_time,  
                     status=surv_df$censor,      
                     marker = surv_df$RAD50,     
                     predict.time = cutoff,span = 0.25*nobs^(-0.20) )
plot(SURVROC.rez$FP, SURVROC.rez$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "\n", "AUC = ",round(SURVROC.rez$AUC,3)), 
     ylab="TP",main="RAD50, Method = NNE \n  Year = 1")
abline(0,1)

## METHOD = KM, Kaplan-Meier
SURVROC.kmrez= survivalROC(Stime=surv_df$surv_time,  
                     status=surv_df$censor,      
                     marker = surv_df$RAD50,     
                     predict.time =  cutoff, method="KM")
plot(SURVROC.kmrez$FP, SURVROC.kmrez$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "\n", "AUC = ",round(SURVROC.kmrez$AUC,3)), 
     ylab="TP",main="RAD50 , Method = KM \n Year = 1")
abline(0,1)

#run for all 35
coxnet.genes <- read_csv("res_coef_coxnet_names.csv")
coxnet.genes <- coxnet.genes$x
coxnet.df <- surv_df[, (colnames(surv_df) %in%coxnet.genes)]
dim(coxnet.df)
time <- surv_df$surv_time
event <- surv_df$censor
rez_list <- apply(coxnet.df, 2, survivalROC, Stime = time, status = event, predict.time = cutoff, method="KM")
#runs seprate plots
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
