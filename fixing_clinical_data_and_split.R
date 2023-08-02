#scriptas kuris sutvarko phenodata i reikalingus stulpelius ir nusplitina mano trai. test ir weird setus
#uzduotys:
#secondary navikus ir prior treatment atsiskirt i weird
#atsiskiri jei 3 data train, test ir weird, kaip minimum vektoriai atskiri su vardais
#parsinti ta keista  treatmentstulpeli 
setwd("~/tcga-ov-data")
# Load packages
library("TCGAbiolinks")
library("caret")
library("dplyr")
library("tdyr")

tcga_data = readRDS(file = "tcga_data.RDS") #pakeiciau pavadinima
pheno <- readRDS(file = "/home/ieva/tcga-ov-data//pheno.RDS")
#####################################################################
#apsiziuriu data
table(pheno$figo_stage, pheno$vital_status, useNA='a')
table(pheno$days_to_death, pheno$vital_status, useNA='a')
ggplot(pheno, aes(x=figo_stage  , y=days_to_death)) +
  geom_boxplot()

#neinformatyvus stulpeliai (visi na, visi vienodi) tumor_descriptor, state, is_ffpe, tissue_type,
#days_to_diagnosis, last_known_disease_status, tissue_or_organ_of_origin, prior_malignancy, 
# classification_of_tumor, tumor_grade, progression_or_recurrence, alcohol_history, gender
#ethnicity, disease_type, primary_site, project_id, name, releasable, released, days_to_sample_procurement, sample_type_id

#informative stulpeliai (skaiciai) intermediate_dimention, shortest_dimention, longest_dimention,
#days_to_last_follow_up, age_at_diagnosis, year_of_diagnosis, age_at_index, days_to_birth, 
# year_of_birth, days_to_death, year_of_death, 
#informative stulpeliai (factors) figo_stage, race, vital_status

#susipiciaus visi like zemiau
table(pheno$synchronous_malignancy, useNA='a')#pas mane irgi 2 buvo

table(pheno$primary_diagnosis, useNA='a') #1 cystadenocarcinoma ir 4 papillary visos HGSOC yra todel geriau neatsikirt
table(pheno$morphology, useNA='a') #same as primary diagnosis 
table(pheno$icd_10_code, useNA='a') #same as primary diagnosis

table(pheno$prior_treatment, useNA='a') #tik vienas prior treatment #noriu ismest ta 1 zmogu

table(pheno$site_of_resection_or_biopsy, useNA='a') #tik 6 "Specified parts of peritoneum" 
table(pheno$days_to_collection, useNA='a')  #Sample can be collected can be prospectively or retrospectively. This can be a negative value for samples taken retrospectively.
table(pheno$initial_weight, useNA='a') #sample, specimen weight in grams
table(pheno$oct_embedded, useNA='a') #A boolean value indicating whether the Optimal Cutting Temperature compound (OCT) is used to embed tissue samples prior to frozen sectioning on a microtome-cryostat.
table(pheno$preservation_method, useNA='a') #one is ffpe after all
table(pheno$is_ffpe, useNA='a') #nors cia raso kad ner ffpe, bet paliksiu does not matter

#noriu pasalinti visus stulpelius kurie mane erzina nes nera informatyvus
non_informative <- c("tumor_descriptor", "state", "is_ffpe", "tissue_type",
                     "days_to_diagnosis", "last_known_disease_status", "tissue_or_organ_of_origin", "prior_malignancy", 
                     "classification_of_tumor", "tumor_grade", "progression_or_recurrence", "alcohol_history", "gender",
                     "ethnicity", "disease_type", "primary_site", "project_id", "name", "releasable", "released", "days_to_sample_procurement", "treatments", "sample_type_id", "composition")
non_informative %in% colnames(pheno) # true if no spelling mistakes were made
new_pheno <- pheno[ ,!names(pheno) %in% non_informative]
##############################################################################
#reikia isparsinti treatments 
treatments_dataframe <- bind_rows(pheno$treatment, .id = "id")
#reikia kad treatmentu dataframe butu vardai kurie po to sutaps su pheno
treatments_dataframe$patient_names<- treatments_dataframe$submitter_id  
treatments_dataframe$patient_names <-gsub("_treatment.*", "", treatments_dataframe$patient_names)
#dabar tikrinu kurie stulpeliai useless
table(treatments_dataframe$days_to_treatment_end, useNA='a') # visi na, drop
table(treatments_dataframe$days_to_treatment_start, useNA='a') # visi na, drop
table(treatments_dataframe$treatment_id, useNA='a') # all unique, maybe drop
table(treatments_dataframe$submitter_id, useNA='a') #perkeliau i patient names ir nusplitinau gala
table(treatments_dataframe$treatment_type, useNA='a')  #kiekvienam zmogui yra po radiation ir parmasiutical eilute
table(treatments_dataframe$regimen_or_line_of_therapy, useNA='a') # visi na, drop
table(treatments_dataframe$treatment_effect, useNA='a')  # visi na, drop
table(treatments_dataframe$therapeutic_agents, useNA='a')  # visi na, drop
table(treatments_dataframe$initial_disease_status, useNA='a')  # visi na, drop
table(treatments_dataframe$treatment_intent_type, useNA='a')  # visi na, drop
table(treatments_dataframe$treatment_anatomic_site, useNA='a')  # visi na, drop
table(treatments_dataframe$treatment_outcome, useNA='a')  # visi na, drop
table(treatments_dataframe$state, useNA='a')  # visi resleased, drop
drop_from_treamtment <- c("days_to_treatment_end", "days_to_treatment_start", "treatment_id",
                          "regimen_or_line_of_therapy", "treatment_effect", "therapeutic_agents", "initial_disease_status", 
                          "treatment_intent_type", "treatment_anatomic_site", "treatment_outcome", "state")
drop_from_treamtment %in% colnames(treatments_dataframe) # true if no spelling mistakes were made
treatments_dataframe <- treatments_dataframe[ ,!names(treatments_dataframe) %in% drop_from_treamtment]

###############################################################################
#pasalinu be treatment
table(treatments_dataframe$treatment_or_therapy, useNA='a') #429 yes kiti no ar na, todel bandau nufiltruot tik yes
had_treatment <- subset(treatments_dataframe, treatment_or_therapy == 'yes') #nudropinau
table(had_treatment$patient_names, useNA='a') #rodo kad yra double kai kur
had_treatment$patient <- had_treatment$patient_names #reikes po to sumerginimui

#group_by(treatments_dataframe, patient_names) %>% summarize(ppp=paste(sort(treatment_or_therapy), collapse="_"), n=n()) %>% group_by(ppp) %>% summarize(n=n())# jul

had_treatment <- had_treatment %>% group_by(patient_names) %>% mutate(undergone_treatments=paste(sort(treatment_type), collapse="_"))   #galutinej eilutej treatment eiga
had_treatment_collaped <- had_treatment[!duplicated(had_treatment$patient_names), ] # removinu duplications, susimazins dataframas

################################################################################
#merge treatments ir normal pheno data frame (sumazinta)
pheno_final <- merge(x = new_pheno, y = had_treatment_collaped, by = "patient", all = TRUE) 

#nusplitint weird: definition ir treatments
split_by_definition <- split(pheno_final, f = pheno_final$definition, drop = T)
weird_cases <- split_by_definition$`Recurrent Solid Tumor` # weird grupe 7 zmones
primary_tumor <- split_by_definition$`Primary solid Tumor`

split_prior_treatment <- split(primary_tumor, f = primary_tumor$prior_treatment, drop = T)
primary_tumor <- split_prior_treatment$No
weird_cases <- rbind(weird_cases, split_prior_treatment$Yes) 

###################################################################################
# Set (random-number-generator) seed so that results are consistent between runs
set.seed(18)
train_ids <-rbinom(nrow(primary_tumor), size = 1, prob = 0.8) ==1

x_train = primary_tumor[train_ids, ]
x_test  = primary_tumor[!train_ids, ]

#save 3 failus i csv
write.csv(weird_cases, "weird_cases.csv" )
write.csv(x_train, "train_data.csv" )
write.csv(x_test, "test_data.csv" )