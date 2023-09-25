#Experimental plots: tcga lasso, pca, volcano, ROC (bet ne survroc)
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
#isikelti duomenys: gtex vs tcga train counts, only protein genes from step5
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_counts_train <- data.frame(gtex_counts_train)
gtex_counts_train$grupe <- substr(rownames(gtex_counts_train), 1, 4)
gtex_counts_train$grupe <- as.factor(gtex_counts_train$grupe)
#isikelti duomenys: is side step7: tik tcga train klinikiniai is XENA ir TCGA duombaziu
pheno_train <- readRDS("tcga_pheno_train.RDS")
#isikelti duomenys: only tcga train counts, be gtex is side step 7 kaip ir klinikiniai
tcga_counts_train <- readRDS("tcga_counts_train.RDS")
tcga_counts_train <- as.data.frame(t(tcga_counts_train))
tcga_counts_train$barcode <- rownames(tcga_counts_train)
tcga <- right_join(pheno_train, tcga_counts_train, by = "barcode")
#EXPERIMENTAL PLOT1: tik ansksciau GDL laboratorijoje tirti genai
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

################################################################################
#EXPERIMENTAL PLOT2: paveikslai for tik tcga lassos is step7
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

##############################################################################
#countmatrix for TCGA and GTEX
gene.matrix <- as.matrix(gtex_counts_train[1:13680])
dim(gene.matrix)
#only the elastic net genes ~200
elastic.genes <- readRDS("gtcga_elastic.RDS")
elastic.matrix <- gene.matrix[, (colnames(gene.matrix) %in% elastic.genes)]
dim(elastic.matrix)
#only the coxnet genes ~35
coxnet.genes <- read_csv("res_coef_cox_names.csv")
coxnet.genes <- coxnet.genes$x
coxnet.genes
coxnet.matrix <- gene.matrix[, (colnames(gene.matrix) %in% coxnet.genes)]
dim(coxnet.matrix) #489  

################################################################################
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
#part1: po norm
#duomenys is STEP4 normalization, VISI meginiai, bet per normalizacija maziau genu
nomcounts_gtcga <- readRDS("mrna_voom_protein.RDS")
nomcounts_gtcga <- t(nomcounts_gtcga) #595 13681
snames_n = rownames(nomcounts_gtcga)
group_n = substr(snames_n, 1, 4) 
############################################################################
#side project picture: PCA of norma gtex and tcga
group_f <- factor(group_n, levels = c("GTEX", "TCGA"))
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
#duomenys turi buti nenormalizuoti: counts_gtcga 
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
group <- cox_counts$grupe
a=multipleROC(group~RAD50,data=cox_counts)
b=multipleROC(group~TTC4,data=cox_counts)
c=multipleROC(group~RPS10,data=cox_counts)
d=multipleROC(group~CLDN4,data=cox_counts)
e=multipleROC(group~ARPC1B,data=cox_counts)
f=multipleROC(group~PKP3,data=cox_counts)
g=multipleROC(group~GRB7,data=cox_counts)
plot_ROC(list(a,b,c,d,e,f,g),show.eta=FALSE,show.sens=FALSE)
plot_ROC(list(a,b,c,d,e,f,g))+facet_grid(no~.)
