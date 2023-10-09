#Multivariable cox?
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)
setwd("~/rprojects/TCGA-OV-data")
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_genes <- readRDS("gtcga_elastic.RDS")
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 

gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train))
gtex_filtered_counts_train2 <- gtex_filtered_counts_train2 %>% dplyr::select(starts_with("TCGA")) 
tcga_counts <- t(gtex_filtered_counts_train2)
dim(tcga_counts)
#336 samples #214 genes
#clinical data
pheno_train <- readRDS("tcga_pheno_train.RDS")
dim(pheno_train)
#336 samples, 69 variables
pheno_train$grade <- pheno_train$neoplasmhistologicgrade
pheno_train$grade[pheno_train$grade %in% c("G4", "GB", "GX")] <- NA 
# survival df
clin_df = pheno_train[,
                      c("barcode",
                        "vital_status",
                        "days_to_death",
                        "days_to_last_follow_up",
                        "neoplasmhistologicgrade",
                        "figo_stage",
                        "age_at_diagnosis", 
                        "tumorresidualdisease", 
                        "intermediate_dimension", "grade")]

# create a new boolean variable that has TRUE for dead patients and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)
# of the name, eg "stage iiia" would become simply "stage iii"
clin_df$tumor_stage = gsub("[ABC]$", "", clin_df$figo_stage)
clin_df[which(clin_df$tumor_stage == "Stage I"), "tumor_stage"] = NA
table(clin_df$tumor_stage, useNA = "a")

tcga_counts <- as.data.frame(tcga_counts)
tcga_counts$barcode <- rownames(tcga_counts)
clin_df_joined <- left_join(clin_df, tcga_counts, by = "barcode")
clin_df_joined$censor <- as.integer(as.logical(clin_df_joined$deceased))
clin_df_joined <- as.data.frame(clin_df_joined)



covariate_names <- c(age_at_diagnosis="Age at Dx",
                     tumor_stage="Figo stage",
                     grade = "Grade",
                     tumorresidualdisease = "Residual disease",
                     intermediate_dimension = "Tumor intermediate dimention", 
                     EXO1  = "EXO1",
                      RAD50 = "RAD50",
                      PPT2 = "PPT2",
                      LUC7L2 = "LUC7L2",
                      PKP3 = "PKP3",
                      CDCA5 = "CDCA5",
                      ZFPL1 = "ZFPL1",
                      VPS33B = "VPS33B",
                      GRB7 = "GRB7",
                      TCEAL4 = "TCEAL4")

clin_df_joined %>%
  mutate(tumor_stage = as.factor(tumor_stage), 
         tumorresidualdisease = as.factor(tumorresidualdisease), 
         grade = as.factor(grade)) %>%
  analyse_multivariate(vars(overall_survival, censor),
                       covariates = vars(age_at_diagnosis,tumor_stage, 
                                         tumorresidualdisease,
                                         grade, intermediate_dimension),
                       covariate_name_dict = covariate_names,
                       reference_level_dict=c(tumorresidualdisease="No Macroscopic disease")) ->
  result2
print(result2)
png("figures/covaites_clinical.png", width=600, height=500)
forest_plot(result2)
dev.off()

clin_df_joined %>%
  mutate(tumor_stage = as.factor(tumor_stage), 
         tumorresidualdisease = as.factor(tumorresidualdisease), 
         grade = as.factor(grade),
         EXO1 = as.numeric(EXO1),
          RAD50 = as.numeric(RAD50),
         PPT2 = as.numeric(PPT2),
         LUC7L2 = as.numeric(LUC7L2),
         PKP3 = as.numeric(PKP3),
         CDCA5 = as.numeric(CDCA5),
         VPS33B = as.numeric(VPS33B),
         TCEAL4 = as.numeric(TCEAL4),
         ZFPL1 = as.numeric(ZFPL1),
         GRB7 = as.numeric(GRB7)
) %>%
  analyse_multivariate(vars(overall_survival, censor),
                       covariates = vars(age_at_diagnosis,tumor_stage, 
                                         tumorresidualdisease,
                                         grade, intermediate_dimension, 
                                         EXO1,
                                          RAD50 ,
                                          PPT2 ,
                                          LUC7L2 ,
                                          PKP3 ,
                                          CDCA5 ,
                                          ZFPL1 ,
                                          VPS33B,
                                          GRB7 ,
                                          TCEAL4 ),
                       covariate_name_dict = covariate_names,
                       reference_level_dict=c(tumorresidualdisease="No Macroscopic disease")) ->
  result9
print(result9)
png("figures/covaites_clinical_genes.png", width=800, height=500,)
forest_plot(result9,  labels_displayed = c("endpoint", "factor", "n"),)
dev.off()

#Age HR=1?
ggplot(data=clin_df_joined,  aes(x = age_at_diagnosis, y = overall_survival))+
  geom_point() +
  geom_smooth()

#age ne dienom, o metais?
clin_df_joined <-clin_df_joined %>%
  mutate(age_at_diagnosis_ys = age_at_diagnosis/365)

ggplot(data=clin_df_joined,  aes(x = age_at_diagnosis_ys, y = overall_survival))+
  geom_point() +
  geom_smooth()

age.cox <- coxph(Surv(overall_survival, censor) ~ age_at_diagnosis_ys, data = clin_df_joined)
summary(age.cox)
