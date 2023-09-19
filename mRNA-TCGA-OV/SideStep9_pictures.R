#I want some plots
library(ggplot2)
library(viridis)
library(reshape2)
setwd("~/rprojects/TCGA-OV-data") #wsl
library(tidyverse)
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


gtex_tcga_genes <- c( "TTC4"  ,  "SLC39A1" ,"TMEM110", "RAD50" ,  "ANKHD1"  ,"ZBTB9" ,  "RPS10" ,
   "CLDN4" ,  "PFDN5" ,  "PAGR1"  ,"RNASEK" , "GPS2"   , "RTEL1" , "grupe")
 gtex_counts_train %>% select(gtex_tcga_genes) %>%
  pivot_longer(., cols = c(gtex_tcga_genes[1:13]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = grupe)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
     axis.text.x = element_text(face = "italic"))+
   guides(fill=guide_legend(title="Study"))
 
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
library(RColorBrewer)
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