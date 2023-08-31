#Step 1.2: clean pheno data
library(tidyverse)
setwd("~/rprojects/TCGA-OV-data")
pheno <- readRDS("pheno.RDS")
#Step 2.1 parse and clean pheno data
non_informative <- c("tumor_descriptor", "state", "is_ffpe", "tissue_type",
                     "days_to_diagnosis", "last_known_disease_status", "tissue_or_organ_of_origin", "prior_malignancy", 
                     "classification_of_tumor", "tumor_grade", "progression_or_recurrence", "alcohol_history", "gender",
                     "ethnicity", "disease_type", "primary_site", "project_id", "name", "releasable",
                     "released", "days_to_sample_procurement", "treatments", "sample_type_id", "composition")
new_pheno <- pheno[ ,!names(pheno) %in% non_informative] #new pheno data with no annoying colums
treatments_dataframe <- bind_rows(pheno$treatment, .id = "id") #reikia kad treatmentu dataframe butu vardai kurie po to sutaps su pheno
treatments_dataframe$patient_names<- treatments_dataframe$submitter_id  
treatments_dataframe$patient_names <-gsub("_treatment.*", "", treatments_dataframe$patient_names)
drop_from_treamtment <- c("days_to_treatment_end", "days_to_treatment_start", "treatment_id",
                          "regimen_or_line_of_therapy", "treatment_effect", "therapeutic_agents", "initial_disease_status", 
                          "treatment_intent_type", "treatment_anatomic_site", "treatment_outcome", "state")
drop_from_treamtment %in% colnames(treatments_dataframe) 
treatments_dataframe <- treatments_dataframe[ ,!names(treatments_dataframe) %in% drop_from_treamtment]
#every person has 2 rows for treatments except for those who had no tretment
#here I drop the rows where the reatment was not given (indicated by no/na in treatment_or_therapy)
table(treatments_dataframe$treatment_or_therapy, useNA='a') #429 yes kiti no ar na, todel bandau nufiltruot tik yes
had_treatment <- subset(treatments_dataframe, treatment_or_therapy == 'yes') #new treatments dataframe
table(had_treatment$patient_names, useNA='a') #this showes that some are double
had_treatment$patient <- had_treatment$patient_names #add patient names for easy merging
#paste treatment types to have 1 row per person
had_treatment <- had_treatment %>% group_by(patient_names) %>% mutate(undergone_treatments=paste(sort(treatment_type), collapse="_"))   #galutinej eilutej treatment eiga
had_treatment_collaped <- had_treatment[!duplicated(had_treatment$patient_names), ] # remove duplications, resulting smaller dataframe
#merge treatments and normal pheno data frame (final pheno dataframe)
pheno_final <- merge(x = new_pheno, y = had_treatment_collaped, by = "patient", all = TRUE) #429 pac ir 47 clinikiniai
rm(new_pheno)
rm(treatments_dataframe) 
rm(had_treatment)
rm(had_treatment_collaped)
rm(drop_from_treamtment)
rm(non_informative)

saveRDS(pheno_final, "pheno_no_empty_data.RDS")

