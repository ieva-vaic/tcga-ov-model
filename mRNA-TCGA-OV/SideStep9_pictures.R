#I want some plots
library(ggplot2)
library(viridis)
library(reshape2)
setwd("~/rprojects/TCGA-OV-data") #wsl
library(tidyverse)
library(RColorBrewer)
library('DESeq2')
library(pROC)
library(ComplexHeatmap)
library(circlize)
#because the normalized hplot was bad this might not be the most appropirate?
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_counts_train <- data.frame(gtex_counts_train)
gtex_counts_train$grupe <- substr(rownames(gtex_counts_train), 1, 4)
gtex_counts_train$grupe <- as.factor(gtex_counts_train$grupe)
#PICTURE FOR MY GENES
 mano_genes <- c( "PPP2R1A" ,"ARID1A", "CTNNB1", "FBXW7",  "NOTCH1" ,"NOTCH2", "NOTCH3" ,  "NOTCH4"  ,
                          "HES1"  ,"JAG2" , "DLL1",  "HOPX" , "grupe")
 gtex_counts_train %>% select(mano_genes) %>%
   pivot_longer(., cols = c(mano_genes[1:12]), names_to = "Genes", values_to = "EXPR") %>%
   ggplot(aes(x = Genes, y = EXPR, fill = grupe)) +
   geom_boxplot()+
   ylab("Normalized expression")+
   theme(
     axis.text.x = element_text(face = "italic"))+
   guides(fill=guide_legend(title="Study"))

pheno_train <- readRDS("tcga_pheno_train.RDS")
tcga_counts_train <- readRDS("tcga_counts_train.RDS")
tcga_counts_train <- as.data.frame(t(tcga_counts_train))
tcga_counts_train$barcode <- rownames(tcga_counts_train)
tcga <- right_join(pheno_train, tcga_counts_train, by = "barcode")
#PICRTURE FOR GRADE
tcga$grade <- tcga$neoplasmhistologicgrade
tcga$grade[tcga$grade %in% c("G4", "GB", "GX")] <- NA 
grade_glm_genes <- readRDS("tcga_grade_45.RDS") 
grade_genes <- c( grade_glm_genes,
                  "grade") #i need to add group so that it will be sorted
tcga %>% select(grade_genes) %>%
  na.omit() %>%
  pivot_longer(., cols = c(grade_genes[1:45]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = grade)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="Study"))
#PICTURE FOR STAGE THAT ONE GENE
table(tcga$clinicalstage2, useNA = "a")
tcga$clinicalstage2[tcga$clinicalstage2 =="Stage I"] <- NA
tcga$stage <- tcga$clinicalstage2
stage_genes <- c( "PRKRA",
                  "stage") #i need to add group so that it will be sorted
tcga %>% select(stage_genes) %>%
  na.omit()%>%
  ggplot(aes(x=stage, y=PRKRA, fill=stage))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  xlab("Stage")+
  ylab(substitute(paste(italic("PRKRA"))))
#PICTRUE FOR LYPHOVASCULAR INVASION
table(tcga$lymphaticinvasion, useNA = "a")
lumph_genes <- c( "MILR1", "VAV1",  "PLCB4", "VSIG4",
                  "lymphaticinvasion") #i need to add group so that it will be sorted
tcga %>% select(lumph_genes) %>%
  na.omit() %>%
  pivot_longer(., cols = c(lumph_genes[1:4]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = as.factor(lymphaticinvasion))) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="lypmhovascular invasion"))
#PICTRUE FOR residual disease left or not
table(tcga$tumorresidualdisease, useNA = "a")
res_dis_genes <- readRDS("tcga_resdis_13.RDS")
res_dis_genes <- c(res_dis_genes, "tumorresidualdisease")
tcga %>% select(res_dis_genes) %>%
  na.omit() %>%
  pivot_longer(., cols = c(res_dis_genes[1:13]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = as.factor(tumorresidualdisease))) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="tumorresidualdisease"))
############################################################################
#heatmap between TCGA and GTEX
gene.matrix <- as.matrix(gtex_counts_train[1:13680])
dim(gene.matrix)
# Example: grouping from the first letter:
my_group <- as.numeric(as.factor(gtex_counts_train$grupe))
colSide <- brewer.pal(9, "Set1")[my_group]

jpeg(file="figures/heatmap15.jpeg", height=2000, width=5000) #išsaugijimui didesniu formatu
heatmap(gene.matrix, Colv = NA, Rowv = NA, scale="col", col=rev(brewer.pal(9,"RdBu")), RowSideColors=colSide)
dev.off()
#man looks weird ar tbh, no difrerences
elastic.genes <- readRDS("gtcga_elastic.RDS")
#only the elastic net genes ~200
elastic.matrix <- gene.matrix[, (colnames(gene.matrix) %in% elastic.genes)]
dim(elastic.matrix)
jpeg(file="figures/heatmap_elastic.jpeg", height=2000, width=3000) #išsaugijimui didesniu formatu
heatmap(elastic.matrix, Colv = NA, scale="col", col=rev(brewer.pal(9,"RdBu")), RowSideColors=colSide)
dev.off()
#only the coxnet genes ~35
coxnet.genes <- read_csv("res_coef_cox_names.csv")
coxnet.genes <- coxnet.genes$x
coxnet.genes
coxnet.matrix <- gene.matrix[, (colnames(gene.matrix) %in% coxnet.genes)]
dim(coxnet.matrix) #489  33 -> 2 pamesti bes reikejo keist vardus vienas FAM83H_AS1 kitas PPT2_f
jpeg(file="figures/heatmap_coxnet.jpeg", height=700, width=700) #išsaugijimui didesniu formatu
heatmap(coxnet.matrix, Colv = NA, scale="col", col=rev(brewer.pal(9,"RdBu")), RowSideColors=colSide)
dev.off()
################################################################################
#side step Heatmaps pries ir po normalizavimo
#part1: pries norm
counts_gtcga <- readRDS("gtcga_final_counts2.RDS") 
counts_gtcga <- data.matrix(counts_gtcga) #19193   596
which(colnames(counts_gtcga)=="TCGA-13-1499-01A-01R-1565-13") #297
counts_gtcga <- counts_gtcga[, -297] #19193   595 liko 
snames = colnames(counts_gtcga)
group = substr(snames, 1, 4) 
group <- as.numeric(as.factor(group))
colSide <- brewer.pal(9, "Set1")[group] #dabar row zmones turetu buy
counts_gtcga <- t(counts_gtcga)
#clusters from both sides no norm
png(file="figures/nonorm_clustbothsides.png", height=700, width=700) #išsaugijimui didesniu formatu
heatmap(counts_gtcga, scale="col", col=rev(brewer.pal(9,"RdBu")), RowSideColors=colSide)
dev.off()
#part1: po norm
nomcounts_gtcga <- readRDS("mrna_voom_protein.RDS")
nomcounts_gtcga <- t(nomcounts_gtcga)
png(file="figures/normcountsheatmap.png", height=700, width=700) #išsaugijimui didesniu formatu
heatmap(nomcounts_gtcga, Colv = NA, scale="col", col=rev(brewer.pal(9,"RdBu")), RowSideColors=colSide)
dev.off()
#clusters from both sides: normcounts
png(file="figures/normcountsheatmap2.png", height=700, width=700) #išsaugijimui didesniu formatu
heatmap(nomcounts_gtcga, scale="col", col=rev(brewer.pal(9,"RdBu")), RowSideColors=colSide)
dev.off()
################################################################################
#side 2 PCA?
group_f <- factor(group, levels = c("GTEX", "TCGA"))
colType <- c("blue", "red2")[group_f]
pchType <- c(18, 16)[group_f]
pca <- prcomp(nomcounts_gtcga)
plot(
  pca$x,
  col = colType,
  cex = 1.0)
legend(
  "bottomright",
  bty = "n",
  c("GTEX", "TCGA"),
  fill = c("blue", "red2"),
  cex = 1.0)
################################################################################
#side 3 volcano plot: p-value agains log-fold change
#tam reikia log foldu ir p-values, tai suskaiciuos man DESEQ2 bet reikia ant raw
#counts_gtcga 
#i need description table sample names and gtcga or gtex
coldata <- data.frame(snames, group)
counts_gtcga <- t(counts_gtcga) #colums are patients
coldata$group <- 
  factor(coldata$group, levels = c("GTEX", "TCGA"))
dds <- DESeqDataSetFromMatrix(countData = counts_gtcga,
                              colData = coldata,
                              design= ~ 0 + group)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
# or to shrink log fold changes association with condition: (suggested by package)
res <- lfcShrink(dds, coef="groupTCGA", type="apeglm") #nezinau ar cia ok TCGA naudot
res_df <- as.data.frame(res) 
#volcano plot I wanted
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
#Adjusted P values (FDR Q values)
with(res_df, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(res_df, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(res_df$pvalue[res_df$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
###############################################################################
#redoing pictures: only the train normalized data
#First I wanna see the top gene (The ~35 coxnet after elastic net) boxplots 
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_counts_train <- data.frame(gtex_counts_train)
gtex_counts_train$grupe <- substr(rownames(gtex_counts_train), 1, 4)
gtex_counts_train$grupe <- as.factor(gtex_counts_train$grupe)
#filter for ~35 coxnet after getting the train data
res_coef_cox_names <- read.csv("res_coef_coxnet_names.csv")
res_coef_cox_names <- res_coef_cox_names$x
res_coef_cox_names <- c(res_coef_cox_names, "grupe") #cia reik priappendint kad nufiltruotu ka reik
#now for the boxplots
gtex_counts_train %>% select(res_coef_cox_names) %>%
  pivot_longer(., cols = c(res_coef_cox_names[1:38]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = grupe)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle=-90))+
  guides(fill=guide_legend(title="Study"))
suvROCgenes <- c("grupe", "TTC4", "RAD50", "RPS10", 
                 "CLDN4", "ARPC1B", "PKP3", "GRB7")
gtex_counts_train %>% select(suvROCgenes) %>%
  pivot_longer(., cols = c(suvROCgenes[2:8]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = grupe)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle=-90))+
  guides(fill=guide_legend(title="Study"))
###############################################################################
#Tiesiog ROC (ne SURV ROC) bet TCGA vs GTEX ar atspėja genų raiška grupę?
#tutorialas is statquest https://github.com/StatQuest/roc_and_auc_demo/blob/master/roc_and_auc_demo.R
cox_counts <- gtex_counts_train %>% select(res_coef_cox_names)
glm.fit=glm(grupe ~ RAD50, family=binomial, data= cox_counts)
roc(cox_counts$grupe, glm.fit$fitted.values, plot=TRUE, 
    legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage",
    ylab="True Postive Percentage", col="#377eb8", lwd=4)
roc.info <- roc(cox_counts$grupe, glm.fit$fitted.values, legacy.axes=TRUE)
str(roc.info)
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)
head(roc.df) ## head() will show us the values for the upper right-hand corner
## of the ROC graph, when the threshold is so low 
## (negative infinity) that every single sample is called "obese".
## Thus TPP = 100% and FPP = 100%

tail(roc.df) ## tail() will show us the values for the lower left-hand corner
## of the ROC graph, when the threshold is so high (infinity) 
## that every single sample is called "not obese". 
## Thus, TPP = 0% and FPP = 0%

## We can calculate the area under the curve...
roc(cox_counts$grupe, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE,
    xlab="False Positive Percentage", ylab="True Postive Percentage",
    col="#377eb8", lwd=4, print.auc=TRUE, main="RAD50")
## ...and the partial area under the curve.
roc(cox_counts$grupe, glm.fit$fitted.values, 
    plot=TRUE, legacy.axes=TRUE, percent=TRUE,
    xlab="False Positive Percentage", ylab="True Postive Percentage",
    col="#377eb8", lwd=4, print.auc=TRUE, print.auc.x=45,
    partial.auc=c(100, 90), auc.polygon = TRUE,
    auc.polygon.col = "#377eb822")
## ROC for random forest
library(randomForest) 
rf.model <- randomForest(factor(cox_counts$grupe) ~ cox_counts$RAD50)
roc(cox_counts$grupe, rf.model$votes[,1], plot=TRUE,
    legacy.axes=TRUE, percent=TRUE,
    xlab="False Positive Percentage", 
    ylab="True Postive Percentage",
    col="#4daf4a", lwd=4, print.auc=TRUE)

#############################################################################
#all gene ROCs, tutorial:https://github.com/cardiomoon/multipleROC
library(multipleROC)
a=multipleROC(group~RAD50,data=cox_counts)
b=multipleROC(group~TTC4,data=cox_counts)
c=multipleROC(group~RPS10,data=cox_counts)
d=multipleROC(group~CLDN4,data=cox_counts)
e=multipleROC(group~ARPC1B,data=cox_counts)
f=multipleROC(group~PKP3,data=cox_counts)
g=multipleROC(group~GRB7,data=cox_counts)
plot_ROC(list(a,b,c,d,e,f,g),show.eta=FALSE,show.sens=FALSE)
plot_ROC(list(a,b,c,d,e,f,g))+facet_grid(no~.)
###############################################################################
#Top ~35 genai su klinikiniais
#nuo virsaus pradedant
#grade
grade_genes38 <- c( res_coef_cox_names,
                  "grade") #i need to add group so that it will be sorted
tcga %>% select(grade_genes38) %>%
  na.omit() %>%
  pivot_longer(., cols = c(res_coef_cox_names[1:38]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = grade)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="GRADE"))
grade_genes7 <- c( "RAD50", "TTC4", "RPS10", "CLDN4", "ARPC1B", "PKP3", "GRB7",
                    "grade") #i need to add group so that it will be sorted
tcga %>% select(grade_genes7) %>%
  na.omit() %>%
  pivot_longer(., cols = c(grade_genes7[1:7]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = grade)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="GRADE"))
#stage
stage_genes38 <- c( res_coef_cox_names,
                    "stage") #i need to add group so that it will be sorted
tcga %>% select(stage_genes38) %>%
  na.omit() %>%
  pivot_longer(., cols = c(stage_genes38[1:38]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = stage)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="STAGE"))
stage_genes7 <- c( "RAD50", "TTC4", "RPS10", "CLDN4", "ARPC1B", "PKP3", "GRB7",
                   "stage") #i need to add group so that it will be sorted
tcga %>% select(stage_genes7) %>%
  na.omit() %>%
  pivot_longer(., cols = c(stage_genes7[1:7]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = stage)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="STAGE"))
###############################################################################
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
############################################################################
#