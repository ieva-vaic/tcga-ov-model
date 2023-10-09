#I want some plots: boxplots
setwd("~/rprojects/TCGA-OV-data") #wsl
library(ggplot2)
library(viridis)
library(reshape2)
library(tidyverse)
library(RColorBrewer)
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
#################################################################################
#filter for ~35 coxnet after getting the train data
res_coef_cox_names <- read.csv("res_coef_coxnet_names.csv")
res_coef_cox_names <- res_coef_cox_names$x
res_coef_cox_names <- c(res_coef_cox_names, "grupe") #cia reik priappendint kad nufiltruotu ka reik
#now for the boxplots
png("figures/train_grupes_box.png", width=600, height=500, res=120)
gtex_counts_train %>% select(res_coef_cox_names) %>%
  pivot_longer(., cols = c(res_coef_cox_names[1:10]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = grupe)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle=-90))+
  guides(fill=guide_legend(title="Study"))
dev.off()

###############################################################################
#Top ~35 genai su klinikiniais
#nuo virsaus pradedant
#grade
tcga$grade <- tcga$neoplasmhistologicgrade
tcga$grade[tcga$grade %in% c("G4", "GB", "GX")] <- NA 
res_coef_cox_names <- res_coef_cox_names[-11] #pasalinu grupe
grade_genes38 <- c(res_coef_cox_names,
                  "grade") #i need to add group so that it will be sorted
png("figures/train_grade_box.png", width=600, height=500, res=120)
tcga %>% select(grade_genes38) %>%
  na.omit() %>%
  pivot_longer(., cols = c(res_coef_cox_names[1:10]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = grade)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="GRADE"))
dev.off()

#stage
table(tcga$clinicalstage2, useNA = "a")
tcga$clinicalstage2[tcga$clinicalstage2 =="Stage I"] <- NA
tcga$stage <- tcga$clinicalstage2
stage_genes38 <- c( res_coef_cox_names,
                    "stage") #i need to add group so that it will be sorted
png("figures/train_stage_box.png", width=600, height=500, res=120)
tcga %>% select(stage_genes38) %>%
  na.omit() %>%
  pivot_longer(., cols = c(stage_genes38[1:10]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = stage)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="STAGE"))
dev.off()
###############################################################################
tcga <- tcga %>% mutate(age_diagnosis_yrs = age_at_diagnosis / 365)
hist(tcga$age_diagnosis_yrs) #nelabai galiu grupem skirstyt, nebent gal +-60?
age_genes7 <- c( res_coef_cox_names,
                   "age_diagnosis_yrs", "stage")
age_df <- tcga %>% select(age_genes7) %>%
  na.omit() %>%
  pivot_longer(., cols = c(age_genes7[1:10]), names_to = "Genes", values_to = "EXPR") %>%
  as.data.frame()
age_df$Genes <- as.factor(age_df$Genes)
png("figures/train_age_dotplot.png", width=600, height=500, res=120)
ggplot(data=age_df,  aes(x = age_diagnosis_yrs, y = EXPR))+
  geom_point(size=0.5) +
  geom_smooth(method="lm", se= FALSE, size = 1) +
  facet_wrap(~Genes, scales = "free")+
  theme(strip.text = element_text(face = "italic"))+
  xlab("Age, years") + ylab("Normalized expression")
dev.off()
#clinical data: I wanna create a table with it 
table(tcga$intermediate_dimension, useNA = "a") #nuo 0,3 iki3, 11na
table(tcga$race, useNA = "a") #283 white vs 42oher vs 11na
table(tcga$anatomicneoplasmsubdivision, useNA = "a") #  233bilateral 44left 39right, 20na
table(tcga$lymphaticinvasion, useNA = "a") # 46 no vs yes  89, 201na
table(tcga$newneoplasmeventtype, useNA = "a") #162na
table(tcga$easterncanceroncologygroup, useNA = "a") #267na
table(tcga$primarytherapyoutcomesuccess, useNA = "a") #102na
table(tcga$tumorresidualdisease, useNA = "a") #33na
#residual disease
resdis_genes38 <- c( res_coef_cox_names,
                    "tumorresidualdisease") #i need to add group so that it will be sorted
png("figures/train_res_dis_plot.png", width=800, height=500)
tcga %>% select(resdis_genes38) %>%
  na.omit() %>%
  pivot_longer(., cols = c(stage_genes38[1:10]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = tumorresidualdisease)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle = 90))+
  guides(fill=guide_legend(title="tumor residual disease"))
dev.off()

#intermediate_dimension

dim_genes7 <- c( res_coef_cox_names, "intermediate_dimension", "stage")
dim_df <- tcga %>% select(dim_genes7) %>%
  na.omit() %>%
  pivot_longer(., cols = c(dim_genes7[1:10]), names_to = "Genes", values_to = "EXPR") %>%
  as.data.frame()
dim_df$Genes <- as.factor(dim_df$Genes)
png("figures/train_dim_plot.png", width=600, height=500)
ggplot(data=dim_df,  aes(x = intermediate_dimension, y = EXPR, color= stage))+
  geom_point(size=0.5) +
  geom_smooth(method="lm", se= FALSE, size = 1) +
  facet_wrap(~Genes, scales = "free")+
  theme(strip.text = element_text(face = "italic"))+
  xlab("Intermadiate tumor dimention, mm") + ylab("Normalized expression")
dev.off()
