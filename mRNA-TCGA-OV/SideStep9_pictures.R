#I want some plots
library(ggplot2)
library(viridis)
library(reshape2)
setwd("~/rprojects/TCGA-OV-data") #wsl
library(tidyverse)
library(RColorBrewer)
library('DESeq2')
#because the normalized hplot was bad this might not be the most appropirate?
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_counts_train <- data.frame(gtex_counts_train)
gtex_counts_train$grupe <- substr(rownames(gtex_counts_train), 1, 4)
gtex_counts_train$grupe <- as.factor(gtex_counts_train$grupe)
gtex_counts_train %>% 
  ggplot(aes(x=grupe, y=TTC4, fill=grupe))+
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  xlab("Study")+
  ylab(substitute(paste(italic("TTC4"))))
 
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
 tcga$grade <- tcga$neoplasmhistologicgrade
 tcga$grade[tcga$grade %in% c("G4", "GB", "GX")] <- NA 
grade_genes <- c( "MTMR11" ,  "CERS2"  ,  "RNF144A",  "RBMS1"  ,  "LNPK"  , 
                  "SGCB"   ,  "ENC1"   ,  "PCDHB4" ,  "TNFRSF21" ,"PGM3"   ,  
                  "LGR4"  ,   "KNL1" ,    "CEP152" ,  "ZNF592" ,  "SHCBP1"  , 
                  "SPRED3" ,  "RASIP1" ,  "ALDH16A1" ,"KLK14"   , "TMX4"   , 
                  "MX2", 
                  "grade") #i need to add group so that it will be sorted

tcga %>% select(grade_genes) %>%
  na.omit()%>%
  pivot_longer(., cols = c(grade_genes[1:21]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = grade)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="Study"))

table(tcga$clinicalstage2, useNA = "a")
tcga$clinicalstage2[tcga$clinicalstage2 =="Stage I"] <- NA
tcga$stage <- tcga$clinicalstage2
stage_genes <- c( "GRHL3"  ,  "NRDC"  ,   "TMEM59"  , "CREG1" ,   "SLC1A4"  ,
                  "PRADC1"  , "PDK1"   ,  "PRKRA"   , "CREB1",   "RAB43" ,   
                  "TRIM59"  , "KLHL8"  ,  "HIST1H1E" ,"SBDS" ,    "TSPAN12" , 
                  "OSTF1" ,   "MAPKAP1" , "RASGEF1A","LGR4" ,    "TMEM123" ,
                  "PDGFD" ,   "BACE1" ,   "LRMP" ,    "ABHD4" ,   "ZNF821" , 
                  "SRR" ,     "PPP4R1" , "C18orf32", "ZNF236" ,  "CTDP1" ,  
                  "SLC44A2" , "CEBPG" ,   "CNFN",     "SULT2B1",  "RDH13" ,  
                  "TRIB3"  , "FAM210B" , "PCBP3" ,   "ZNF74" ,   "SLC10A3",  
                  "MT-ND5"  ,
                  "stage") #i need to add group so that it will be sorted

tcga %>% select(stage_genes) %>%
  na.omit()%>%
  pivot_longer(., cols = c(stage_genes[1:41]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = stage)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="Study"))


vital_genes <- c("CD38",  "FAM3C", "PDGFD", "FANCI", "PPL",   "IL34", 
                 "vital_status")
tcga %>% select(vital_genes) %>%
  na.omit()%>%
  pivot_longer(., cols = c(vital_genes[1:6]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = vital_status)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="Study"))

############################################################################
#heatmap between TCGA and GTEX
gene.matrix <- as.matrix(gtex_counts_train[1:13680])
dim(gene.matrix)
# Example: grouping from the first letter:
my_group <- as.numeric(as.factor(gtex_counts_train$grupe))
colSide <- brewer.pal(9, "Set1")[my_group]

jpeg(file="figures/heatmap2.jpeg", height=2000, width=3000) #išsaugijimui didesniu formatu
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
