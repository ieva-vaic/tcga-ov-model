#This script is for: 1) cheking and filtering phenodata (clinical data) 
#2) spliting data into rain, test nad weird sets
#tasks in this script:
#secondary navikus ir prior treatment atsiskirt i weird
#atsiskiri jei 3 data train, test ir weird, kaip minimum vektoriai atskiri su vardais
#parsinti ta keista  treatmentstulpeli 
setwd("~/tcga-ov-data")
# Load packages
library("TCGAbiolinks")
library("caret")
library("dplyr")
library("tidyr")
library("tidyverse")
tcga_data = readRDS(file = "tcga_data.RDS") 
pheno <- readRDS(file = "/home/ieva/tcga-ov-data//pheno.RDS")
#####################################################################
#Chek the data 
table(pheno$figo_stage, pheno$vital_status, useNA='a')
table(pheno$days_to_death, pheno$vital_status, useNA='a')
ggplot(pheno, aes(x=figo_stage  , y=days_to_death)) +
  geom_boxplot()
############################################################
#noninformative  columns (all na, all the same value):
#tumor_descriptor, state, is_ffpe, tissue_type, days_to_diagnosis, 
#last_known_disease_status, tissue_or_organ_of_origin, prior_malignancy,
#classification_of_tumor, tumor_grade, progression_or_recurrence, 
#alcohol_history, gender, ethnicity, disease_type, primary_site, project_id, 
#name, releasable, released, days_to_sample_procurement, sample_type_id
#informative colums (numerical);
#intermediate_dimention, shortest_dimention, longest_dimention, days_to_last_follow_up,
#age_at_diagnosis, year_of_diagnosis, age_at_index, days_to_birth, 
# year_of_birth, days_to_death, year_of_death, initial_weight
#informative colums (factors) figo_stage, race, vital_status

#Inspect susipiciaus colums:
#synchonous maliganacy
table(pheno$synchronous_malignancy, useNA='a')#pas mane irgi 2 buvo
#morphology/diagnosis all three are the same
table(pheno$primary_diagnosis, useNA='a') #1 cystadenocarcinoma ir 4 papillary visos HGSOC yra todel geriau neatsikirt
table(pheno$morphology, useNA='a') #same as primary diagnosis 
table(pheno$icd_10_code, useNA='a') #same as primary diagnosis
#prior treatment
table(pheno$prior_treatment, useNA='a') #tik vienas prior treatment #noriu ismest ta 1 zmogu
# suspicous colums but not worth deletion:
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
                     "ethnicity", "disease_type", "primary_site", "project_id", "name", "releasable",
                     "released", "days_to_sample_procurement", "treatments", "sample_type_id", "composition")
non_informative %in% colnames(pheno) # true if no spelling mistakes were made
new_pheno <- pheno[ ,!names(pheno) %in% non_informative] #new pheno data with no annoying colums
##############################################################################
#Parsing treatments colum
#first convert the treatments colum to dataframe
treatments_dataframe <- bind_rows(pheno$treatment, .id = "id")
#reikia kad treatmentu dataframe butu vardai kurie po to sutaps su pheno
treatments_dataframe$patient_names<- treatments_dataframe$submitter_id  
treatments_dataframe$patient_names <-gsub("_treatment.*", "", treatments_dataframe$patient_names)
#dabar tikrinu kurie stulpeliai useless
table(treatments_dataframe$days_to_treatment_end, useNA='a') # visi na, drop
table(treatments_dataframe$days_to_treatment_start, useNA='a') # visi na, drop
table(treatments_dataframe$treatment_id, useNA='a') # all unique, maybe drop
table(treatments_dataframe$submitter_id, useNA='a') #perkeliau i patient names ir nusplitinau gala
table(treatments_dataframe$treatment_type, useNA='a')  #kiekvienam zmogui yra po radiation ir Pharmaceutical eilute
table(treatments_dataframe$regimen_or_line_of_therapy, useNA='a') # visi na, drop
table(treatments_dataframe$treatment_effect, useNA='a')  # visi na, drop
table(treatments_dataframe$therapeutic_agents, useNA='a')  # visi na, drop
table(treatments_dataframe$initial_disease_status, useNA='a')  # visi na, drop
table(treatments_dataframe$treatment_intent_type, useNA='a')  # visi na, drop
table(treatments_dataframe$treatment_anatomic_site, useNA='a')  # visi na, drop
table(treatments_dataframe$treatment_outcome, useNA='a')  # visi na, drop
table(treatments_dataframe$state, useNA='a')  # visi released, drop
drop_from_treamtment <- c("days_to_treatment_end", "days_to_treatment_start", "treatment_id",
                          "regimen_or_line_of_therapy", "treatment_effect", "therapeutic_agents", "initial_disease_status", 
                          "treatment_intent_type", "treatment_anatomic_site", "treatment_outcome", "state")
drop_from_treamtment %in% colnames(treatments_dataframe) 
treatments_dataframe <- treatments_dataframe[ ,!names(treatments_dataframe) %in% drop_from_treamtment]

###############################################################################
#every person has 2 rows for treatments except for those who had no tretment
#here I drop the rows where the reatment was not given (indicated by no/na in treatment_or_therapy)
table(treatments_dataframe$treatment_or_therapy, useNA='a') #429 yes kiti no ar na, todel bandau nufiltruot tik yes
had_treatment <- subset(treatments_dataframe, treatment_or_therapy == 'yes') #new treatments dataframe
table(had_treatment$patient_names, useNA='a') #this showes that some are double
had_treatment$patient <- had_treatment$patient_names #add patient names for easy merging
#paste treatment types to have 1 row per person
had_treatment <- had_treatment %>% group_by(patient_names) %>% mutate(undergone_treatments=paste(sort(treatment_type), collapse="_"))   #galutinej eilutej treatment eiga
had_treatment_collaped <- had_treatment[!duplicated(had_treatment$patient_names), ] # remove duplications, resulting smaller dataframe

################################################################################
#merge treatments and normal pheno data frame (final pheno dataframe)
pheno_final <- merge(x = new_pheno, y = had_treatment_collaped, by = "patient", all = TRUE) 
#############################################################################
#split weird cases from the rest of the data by definition and treatments
split_by_definition <- split(pheno_final, f = pheno_final$definition, drop = T)
weird_cases <- split_by_definition$`Recurrent Solid Tumor` # weird group 7 people
primary_tumor <- split_by_definition$`Primary solid Tumor` #new new phenodata

split_prior_treatment <- split(primary_tumor, f = primary_tumor$prior_treatment, drop = T)
primary_tumor <- split_prior_treatment$No
weird_cases <- rbind(weird_cases, split_prior_treatment$Yes) 

###################################################################################
# Set (random-number-generator) seed so that results are consistent between runs
set.seed(18) #choose favorite number
train_ids <-rbinom(nrow(primary_tumor), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
x_train = primary_tumor[train_ids, ] #new train data
x_test  = primary_tumor[!train_ids, ] #new test data 
#save 3 files in csv
write.csv(weird_cases, "weird_cases.csv", quote = T, row.names = F) #normal quote F uzdedam bet kadangi kai kur gali but kableliu neveiks
write.csv(x_train, "train_data.csv", quote = T, row.names = F)
write.csv(x_test, "test_data.csv", quote = T, row.names = F )

#split the rest of the data according to split clinical dataframes
# for most normal analyses I want a matrix as it is given by assay function
#option 1: go back and filter original matrix
counts_ov <- assay(tcga_data)

train_counts2 <- counts_ov[, colnames(counts_ov) %in% train_barcodes]
saveRDS(train_counts2, "train_count_matrix.RDS")

#same with other count matrixes
test_counts2 <- counts_ov[, colnames(counts_ov) %in% test_barcodes]
saveRDS(test_counts2, "test_count_matrix.RDS")

weird_counts2 <- counts_ov[, colnames(counts_ov) %in% weird_barcodes]
saveRDS(weird_counts2, "weird_count_matrix.RDS") #this is weird because 1 variable too much

#####################################################################################
#spliting other kinds of count data just in case
#assay data converted to dataframe
counts_ov <- assay(tcga_data)
counts_ov <- as.data.frame(t(counts_ov))
counts_ov$barcode <- rownames(counts_ov)

#counts files in dataframes 
test_barcodes <- x_test[['barcode']]
test_counts <- counts_ov[rownames(counts_ov) %in% test_barcodes, ]

train_barcodes <- x_train[['barcode']]
train_counts <- counts_ov[rownames(counts_ov) %in% train_barcodes, ]

weird_barcodes <- weird_cases[['barcode']]
weird_counts <- counts_ov[rownames(counts_ov) %in% weird_barcodes, ]

#save everthing to cvs
write.csv(test_counts, "test_counts.csv", quote = F, row.names = F)
write.csv(train_counts, "train_counts.csv", quote = F, row.names = F)
write.csv(weird_counts, "weird_counts.csv", quote = F, row.names = F )

####################################################################################
#useless but left for reference if I need semijoin that filters while joining
# #jsemioin counts and filtered clinical
# counts_ov <- tibble::rownames_to_column(counts_ov, "row_names")
# counts_ov$barcode <- counts_ov$row_names
# counts_ov <- counts_ov[!names(counts_ov) %in% 'row_names']
# #"test data"
# test_data_df = list(x_test, counts_ov)
# test_data_df <- test_data_df %>% reduce(semi_join, by='barcode')
# View(test_data_df)
# 
# train_data_df = list(x_train, counts_ov)
# train_data_df <- train_data_df %>% reduce(semi_join, by='barcode')
# View(train_data_df)
# 
# weird_data_df = list(weird_cases, counts_ov)
# weird_data_df <- weird_data_df %>% reduce(semi_join, by='barcode')
# View(weird_data_df)

# write.csv(test_data_df, "test_data_df.csv", quote = F, row.names = F)
# write.csv(train_data_df, "train_data_df.csv", quote = F, row.names = F)
# write.csv(weird_data_df, "weird_data_df.csv", quote = F, row.names = F )


