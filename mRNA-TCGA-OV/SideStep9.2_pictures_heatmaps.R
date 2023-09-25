#Heatmaps only
setwd("~/rprojects/TCGA-OV-data") #wsl
library(ggplot2)
library(viridis)
library(reshape2)
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pROC)
library(ComplexHeatmap)
library(circlize)
#gtex ir tcga sujungti counts, protein coding genes, is step5
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_counts_train <- data.frame(gtex_counts_train)
gtex_counts_train$grupe <- substr(rownames(gtex_counts_train), 1, 4)
gtex_counts_train$grupe <- as.factor(gtex_counts_train$grupe)
#is side step7 tik tcga counts ir klinikiniai:
pheno_train <- readRDS("tcga_pheno_train.RDS")
tcga_counts_train <- readRDS("tcga_counts_train.RDS")
tcga_counts_train <- as.data.frame(t(tcga_counts_train))
tcga_counts_train$barcode <- rownames(tcga_counts_train)
tcga <- right_join(pheno_train, tcga_counts_train, by = "barcode")
############################################################################
#IMPORTANT HEATMAPS
#countmatrix for TCGA and GTEX, train, protein only
gene.matrix <- as.matrix(gtex_counts_train[1:13680])
dim(gene.matrix)
#heatmap between TCGA and GTEX
my_group <- as.numeric(as.factor(gtex_counts_train$grupe))
colSide <- brewer.pal(9, "Set1")[my_group]
jpeg(file="figures/heatmap15.jpeg", height=2000, width=5000) #išsaugijimui didesniu formatu
heatmap(gene.matrix, Colv = NA, Rowv = NA, scale="col", col=rev(brewer.pal(9,"RdBu")), RowSideColors=colSide)
dev.off() #kazkodel mazokai genu uzraso
#heatmap between TCGA and GTEX, tik elastic net genai
elastic.genes <- readRDS("gtcga_elastic.RDS")
#only the elastic net genes ~200
elastic.matrix <- gene.matrix[, (colnames(gene.matrix) %in% elastic.genes)]
dim(elastic.matrix)
jpeg(file="figures/heatmap_elastic.jpeg", height=2000, width=3000) #išsaugijimui didesniu formatu
heatmap(elastic.matrix, Colv = NA, scale="col", col=rev(brewer.pal(9,"RdBu")), RowSideColors=colSide)
dev.off()
#heatmap between TCGA and GTEX, tik coxnet net genai
#only the coxnet genes ~35
coxnet.genes <- read_csv("res_coef_cox_names.csv")
coxnet.genes <- coxnet.genes$x
coxnet.genes
coxnet.matrix <- gene.matrix[, (colnames(gene.matrix) %in% coxnet.genes)]
dim(coxnet.matrix) #489  
jpeg(file="figures/heatmap_coxnet.jpeg", height=700, width=700) #išsaugijimui didesniu formatu
heatmap(coxnet.matrix, Colv = NA, scale="col", col=rev(brewer.pal(9,"RdBu")), RowSideColors=colSide)
dev.off()
################################################################################
#Heatmaps pries ir po normalizavimo
#part1: pries norm
#duomenys: nenormalizuoti gtex ir tcga, VISI meginiai, protein genai, is STEP3
counts_gtcga <- readRDS("gtcga_final_counts2.RDS") 
counts_gtcga <- data.matrix(counts_gtcga) #19193   596
#pasilinamas tas vienas meginis kuri nusalinau:
which(colnames(counts_gtcga)=="TCGA-13-1499-01A-01R-1565-13") #297
counts_gtcga <- counts_gtcga[, -297] #19193   595 liko 
snames = colnames(counts_gtcga)
group = substr(snames, 1, 4) 
group <- as.numeric(as.factor(group))
colSide <- brewer.pal(9, "Set1")[group] #dabar row zmones turetu buy
counts_gtcga <- t(counts_gtcga)
#clusters from both sides no norm
#TAKES LONG TIME, NES CLUST IS ABIEJU PUSIU!
png(file="figures/nonorm_clustbothsides.png", height=700, width=700) #išsaugijimui didesniu formatu
heatmap(counts_gtcga, scale="col", col=rev(brewer.pal(9,"RdBu")), RowSideColors=colSide)
dev.off()
#part1: po norm
#duomenys is STEP4 normalization, VISI meginiai, bet per normalizacija maziau genu
nomcounts_gtcga <- readRDS("mrna_voom_protein.RDS")
nomcounts_gtcga <- t(nomcounts_gtcga) #595 13681
png(file="figures/normcountsheatmap.png", height=700, width=700) #išsaugijimui didesniu formatu
heatmap(nomcounts_gtcga, Colv = NA, scale="col", col=rev(brewer.pal(9,"RdBu")), RowSideColors=colSide)
dev.off()
#TAKES LONG TIME, NES CLUST IS ABIEJU PUSIU!
#clusters from both sides: normcounts
png(file="figures/normcountsheatmap2.png", height=700, width=700) #išsaugijimui didesniu formatu
heatmap(nomcounts_gtcga, scale="col", col=rev(brewer.pal(9,"RdBu")), RowSideColors=colSide)
dev.off()

################################################################################
#complexheatmap
#filter for ~35 coxnet after getting the train data
res_coef_cox_names <- read.csv("res_coef_coxnet_names.csv")
res_coef_cox_names <- res_coef_cox_names$x
#complex heatmaps
top_genes38_mat <- tcga %>% select(c(res_coef_cox_names, "patient.x")) 
rownames(top_genes38_mat) <- top_genes38_mat$patient.x
top_genes38_mat <-  top_genes38_mat[, -39]
top_genes38_mat <- data.matrix(top_genes38_mat)
Heatmap(top_genes38_mat)

top_genes38_mat2 <- gtex_counts_train %>% select(res_coef_cox_names) 
top_genes38_mat2 <- data.matrix(top_genes38_mat2)
Heatmap(top_genes38_mat2)

#small heatmap with clinical
top_genes7_mat <- gtex_counts_train %>% 
  select(c("RAD50", "TTC4", "RPS10", "CLDN4", "ARPC1B", "PKP3", "GRB7")) 
top_genes7_mat$barcode <- rownames(top_genes7_mat)
top_genes7_mat_clin <- left_join(top_genes7_mat, tcga, by="barcode") 
rownames(top_genes7_mat_clin) <- top_genes7_mat_clin$barcode
top_genes7_mat_clin <- top_genes7_mat_clin %>% 
  select(c("RAD50.x", "TTC4.x", "RPS10.x", "CLDN4.x", "ARPC1B.x", "PKP3.x", "GRB7.x",
           "neoplasmhistologicgrade", 
           "figo_stage", "prior_treatment", "age_at_diagnosis", "treatment_or_therapy",
           "lymphaticinvasion", "vitalstatus", "newneoplasmeventtype", 
           "residualtumor")) 
col_fun = colorRamp2(c(-2, 0, 2), c("pink", "white", "cadetblue"))
col_fun(seq(-3, 3))
top_genes7_mat <-  top_genes7_mat[, -8]
top_genes7_mat <- data.matrix(top_genes7_mat)
rownames(top_genes7_mat) = NULL
row_ha = rowAnnotation(study = gtex_counts_train$grupe,
                       grade =top_genes7_mat_clin$neoplasmhistologicgrade, 
                       stage = top_genes7_mat_clin$figo_stage, 
                       status = top_genes7_mat_clin$vitalstatus, 
                       `residual disease` = top_genes7_mat_clin$residualtumor, 
                       gp = gpar(row_names_gp = NULL) )
Heatmap(top_genes7_mat, right_annotation = row_ha, col = col_fun ) 
