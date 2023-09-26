#Multivariable cox?
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot.RDS")
gtex_genes <- readRDS("gtcga_elastic.RDS")
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 

gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train))
gtex_filtered_counts_train2 <- gtex_filtered_counts_train2 %>% dplyr::select(starts_with("TCGA")) 
tcga_counts <- t(gtex_filtered_counts_train2)

#336 samples #221 genes
#clinical data
pheno_train <- readRDS("tcga_pheno_train.RDS")
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
                     RAD50="RAD50",
                     TTC4="TTC4",
                     CLDN4="CLDN4", 
                     RPS10= "RPS10", 
                     ARPC1B = "ARPC1B", 
                     PKP3="PKP3",
                     GRB7="GRB7")

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
forest_plot(result2)

clin_df_joined %>%
  mutate(tumor_stage = as.factor(tumor_stage), 
         tumorresidualdisease = as.factor(tumorresidualdisease), 
         grade = as.factor(grade),
         RAD50 = as.numeric(RAD50),
         TTC4 = as.numeric(TTC4),
         CLDN4 = as.numeric(CLDN4),
         RPS10 = as.numeric(RPS10),
         ARPC1B = as.numeric(ARPC1B),
         PKP3 = as.numeric(PKP3), 
         GRB7 = as.numeric(GRB7)) %>%
  analyse_multivariate(vars(overall_survival, censor),
                       covariates = vars(age_at_diagnosis,tumor_stage, 
                                         tumorresidualdisease,
                                         grade, intermediate_dimension, 
                                         RAD50,
                                         TTC4,
                                         CLDN4, 
                                         RPS10, 
                                         ARPC1B, 
                                         PKP3,
                                         GRB7),
                       covariate_name_dict = covariate_names,
                       reference_level_dict=c(tumorresidualdisease="No Macroscopic disease")) ->
  result
print(result)
forest_plot(result)

clin_df_joined %>%
  mutate(RAD50 = as.numeric(RAD50),
         TTC4 = as.numeric(TTC4),
         CLDN4 = as.numeric(CLDN4),
         RPS10 = as.numeric(RPS10),
         ARPC1B = as.numeric(ARPC1B),
         PKP3 = as.numeric(PKP3), 
         GRB7 = as.numeric(GRB7)) %>%
  analyse_multivariate(vars(overall_survival, censor),
                       covariates = vars( RAD50,
                                          TTC4,
                                          CLDN4, 
                                          RPS10, 
                                          ARPC1B, 
                                          PKP3,
                                          GRB7),
                       covariate_name_dict = covariate_names) ->
  result_genes
print(result_genes)
forest_plot(result_genes)

#choosing diffrerent genes? I wanna choose ne top survROC bet tiesiog top 9 cox univarite genes
# GRB7
# PPT2
# TPM3
# VPS33B
# LUC7L2
# PKP3
# FANCG
# EXO1
# CDCA5


covariate_names <- c(age_at_diagnosis="Age at Dx",
                     tumor_stage="Figo stage",
                     grade = "Grade",
                     tumorresidualdisease = "Residual disease",
                     intermediate_dimension = "Tumor intermediate dimention", 
                     GRB70="GRB70",
                     PPT2="PPT2",
                     TPM3="TPM3",
                     VPS33B="VPS33B",
                     LUC7L2="LUC7L2",
                     PKP3="PKP3",
                     FANCG="FANCG",
                     EXO1="EXO1",
                     CDCA5="CDCA5")
clin_df_joined %>%
  mutate(tumor_stage = as.factor(tumor_stage), 
         tumorresidualdisease = as.factor(tumorresidualdisease), 
         grade = as.factor(grade),
         GRB7 = as.numeric(GRB7),
         PPT2 = as.numeric(PPT2),
         TPM3 = as.numeric(TPM3),
         VPS33B = as.numeric(VPS33B),
         LUC7L2 = as.numeric(LUC7L2),
         PKP3 = as.numeric(PKP3),
         FANCG = as.numeric(FANCG),
         EXO1 = as.numeric(EXO1),
         CDCA5 = as.numeric(CDCA5),
) %>%
  analyse_multivariate(vars(overall_survival, censor),
                       covariates = vars(age_at_diagnosis,tumor_stage, 
                                         tumorresidualdisease,
                                         grade, intermediate_dimension, 
                                         GRB7,
                                         PPT2,
                                         TPM3,
                                         VPS33B,
                                         LUC7L2,
                                         PKP3,
                                         FANCG,
                                         EXO1,
                                         CDCA5),
                       covariate_name_dict = covariate_names,
                       reference_level_dict=c(tumorresidualdisease="No Macroscopic disease")) ->
  result9
print(result9)
forest_plot(result9,  labels_displayed = c("endpoint", "factor", "n"),)

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
#top12, 1se
clin_df_joined %>%
  mutate(tumor_stage = as.factor(tumor_stage), 
         tumorresidualdisease = as.factor(tumorresidualdisease), 
         grade = as.factor(grade),
         KDF1 = as.numeric(KDF1),
         PPT2 = as.numeric(PPT2),
         TPM3 = as.numeric(TPM3),
         VPS33B = as.numeric(VPS33B),
         LUC7L2 = as.numeric(LUC7L2),
         PKP3 = as.numeric(PKP3),
         RAD50 = as.numeric(RAD50),
         EXO1 = as.numeric(EXO1),
         CDCA5 = as.numeric(CDCA5),
         ZFPL1 = as.numeric(ZFPL1),
         GRB7 = as.numeric(GRB7),
         SNX21 = as.numeric(SNX21)
  ) %>%
  analyse_multivariate(vars(overall_survival, censor),
                       covariates = vars(age_at_diagnosis,tumor_stage, 
                                         tumorresidualdisease,
                                         grade, intermediate_dimension, 
                                         KDF1,
                                         PPT2,
                                         TPM3,
                                         VPS33B,
                                         LUC7L2,
                                         PKP3,
                                         RAD50,
                                         EXO1,
                                         CDCA5,
                                         ZFPL1,
                                         GRB7,
                                         SNX21),
                       covariate_name_dict = covariate_names,
                       reference_level_dict=c(tumorresidualdisease="No Macroscopic disease")) ->
  result12
print(result12)
forest_plot(result12,  labels_displayed = c("endpoint", "factor", "n"),)





