#Step 2.1: clean pheno data
library(tidyverse)
setwd("~/rprojects/TCGA-OV-data")
#load pheno data from step1
pheno <- readRDS("pheno.RDS")
#Step 2.1 parse and clean pheno data
#remove non informative colums:
non_informative <- c("shortLetterCode",  "tumor_descriptor", "state", "is_ffpe", "tissue_type",
                     "days_to_diagnosis", "last_known_disease_status", "tissue_or_organ_of_origin", "prior_malignancy", 
                     "classification_of_tumor", "tumor_grade", "progression_or_recurrence", "alcohol_history", "gender",
                     "ethnicity", "disease_type", "primary_site", "project_id", "name", "releasable",
                     "released", "days_to_sample_procurement", "treatments", "sample_type_id", "composition",
                     "pathology_report_uuid", "year_of_diagnosis" , "diagnosis_id", "exposure_id",
                     "demographic_id", "oct_embedded", "preservation_method" )
new_pheno <- pheno[ ,!names(pheno) %in% non_informative] #new pheno data with no annoying colums
#parse treatment column:
treatments_dataframe <- bind_rows(pheno$treatment, .id = "id") #reikia kad treatmentu dataframe butu vardai kurie po to sutaps su pheno
treatments_dataframe$patient_names<- treatments_dataframe$submitter_id  
treatments_dataframe$patient_names <-gsub("_treatment.*", "", treatments_dataframe$patient_names)
drop_from_treamtment <- c("days_to_treatment_end", "days_to_treatment_start", "treatment_id",
                          "regimen_or_line_of_therapy", "treatment_effect", "therapeutic_agents", "initial_disease_status", 
                          "treatment_intent_type", "treatment_anatomic_site", "treatment_outcome", "state",
                          "created_datetime", "updated_datetime", "id", "submitter_id")
drop_from_treamtment %in% colnames(treatments_dataframe) 
treatments_dataframe <- treatments_dataframe[ ,!names(treatments_dataframe) %in% drop_from_treamtment]
#every person has 2 rows for treatments except for those who had no treatment
#here I drop the rows where the treatment was not given (indicated by no/na in treatment_or_therapy)
table(treatments_dataframe$treatment_or_therapy, useNA='a') #429 yes kiti no ar na, todel bandau nufiltruot tik yes
had_treatment <- subset(treatments_dataframe, treatment_or_therapy == 'yes') #new treatments dataframe
table(had_treatment$patient_names, useNA='a') #this showed that some are double
had_treatment$patient <- had_treatment$patient_names #add patient names for easy merging
#paste treatment types to have 1 row per person
had_treatment <- had_treatment %>% group_by(patient_names) %>% mutate(undergone_treatments=paste(sort(treatment_type), collapse="_"))   #galutinej eilutej treatment eiga
had_treatment_collaped <- had_treatment[!duplicated(had_treatment$patient_names), ] # remove duplications, resulting smaller dataframe
#merge treatments and normal pheno data frame (final pheno dataframe)
pheno_final <- merge(x = new_pheno, y = had_treatment_collaped, by = "patient", all = TRUE) #429 pac ir 35 clinikiniai
dim(pheno_final)
#remove extra dfs
rm(new_pheno)
rm(treatments_dataframe) 
rm(had_treatment)
rm(had_treatment_collaped)
rm(drop_from_treamtment)
rm(non_informative)
#rearange pheno_final into pheno_final2
pheno_final2 <- pheno_final %>%
  relocate(definition, .after = sample_type)  %>%
  relocate(primary_diagnosis, .after = icd_10_code) %>%
  relocate(treatment_type, .before = undergone_treatments) %>%
  select(-(submitter_id))%>%
  relocate(days_to_last_follow_up, .before = days_to_death) %>%
  relocate(age_at_diagnosis, .before = age_at_index) %>%
  relocate(bcr_patient_barcode, .before = sample) %>%
  relocate(sample.aux, .before = sample) %>%
  select(-(patient_names)) %>%
  select(-(sample_submitter_id)) 
dim(pheno_final2)

#CLEAN figo_stage (remove letters) 
table(pheno_final2$figo_stage, useNA = "a");
pheno_final2$STAGE = gsub("[ABC]$", "", pheno_final2$figo_stage)
table(pheno_final2$STAGE, useNA = "a") 
#Recode the variable "STAGE" from Stage IA.. to.. -2 -1 0 +1.
pheno_final2$STAGE_num <- pheno_final2$STAGE
pheno_final2$STAGE_num  <- recode(pheno_final2$STAGE_num , "Stage I" = -2, "Stage II" = -1, "Stage III" = 0, "Stage IV" = 1)
table(pheno_final2$STAGE_num , useNA = "a") 

#save TCGA PHENO BE NEREIKALINGŲ DUOMENŲ -> "pheno_no_empty_data.RDS"
dim(pheno_final2)
saveRDS(pheno_final2, "pheno_no_empty_data.RDS")

########################################################################
#also add clinical from XENA
XENAclin <- read_csv("~/rprojects/TCGA-OV-data/00_ClinTraits.csv")
XENAclin <- XENAclin[, -1]
XENAclin$sample.aux <- XENAclin$sampleID
dim(XENAclin) #584 samples 36 vars
XENAclin <- XENAclin[(XENAclin$sample.aux %in% pheno_final2$sample.aux), ] 
dim(XENAclin)#lieka 420
joined_clin <- left_join(pheno_final2, XENAclin, by = "sample.aux" ) 
dim(joined_clin)#nepamirst kad 9 be data
#save APJUNGTAS TCGA IR XENA PHENO -> "joinedTCGA_XENA_clinical.RDS"
saveRDS(joined_clin, "joinedTCGA_XENA_clinical.RDS")
table(joined_clin$figo_stage, useNA="a")
table(joined_clin$clinicalstage, useNA="a")
#chek ids
justIDS <- joined_clin[, c(1:6,35)]
