#Step 2.2.0 remove weird people from tcga
library(tidyverse)
library(SummarizedExperiment)
setwd("~/rprojects/TCGA-OV-data")
#read into tcga and pheno data
pheno_final <- readRDS("joinedTCGA_XENA_clinical.RDS") #or use "pheno_no_empty_data.RDS"
tcga_data <- readRDS("tcga_data.RDS")
tcga_data <- assay(tcga_data)
#fix tcga_data
tcga_counts_t <- t(tcga_data) #praradom colnames
tcga_counts_t <- as.data.frame(tcga_counts_t) #non-numeric
tcga_counts_t$barcode <- rownames(tcga_counts_t)
#splice tcga_data with joinedTCGA_XENA_clinical
tgca_pheno <- right_join(pheno_final, tcga_counts_t, by = "barcode")
dim(tgca_pheno) #429 zmones ir 60729 (60660 genai + 69 clinical)
rownames(tgca_pheno) <- tgca_pheno$barcode

################################################################################
#split off weird cases
#split weird cases from the rest of the data by definition and treatments
mRNA_full1 <- tgca_pheno #temporarily change name
#split off 7 Reacurrent tumors
split_by_definition <- split(mRNA_full1, f = mRNA_full1$definition, drop = T)
weird_cases <- split_by_definition$`Recurrent Solid Tumor` # weird group 7 people
mRNA_full <- split_by_definition$`Primary solid Tumor` #new new phenodata
dim(weird_cases)
#split off 1 case with prior treatment
split_prior_treatment <- split(mRNA_full, f = mRNA_full$prior_treatment, drop = T)
mRNA_full <- split_prior_treatment$No #now the mRNA_full are good cases
weird_cases <- rbind(weird_cases, split_prior_treatment$Yes) 
dim(weird_cases)
#split off 5 non C56.9 cases
table(mRNA_full$primary_diagnosis, useNA="a")
table(mRNA_full$icd_10_code, useNA="a")
split_non_serous <- split(mRNA_full, f = mRNA_full$icd_10_code, drop = T)
mRNA_full <- split_non_serous$C56.9 #now the mRNA_full are good cases
weird_cases <- rbind(weird_cases, split_non_serous$C48.1)
weird_cases <- rbind(weird_cases, split_non_serous$C48.2) 
#chek if splits are correct
dim(weird_cases) # 13 zmones nusišalina
dim(mRNA_full) # 416 lieka iš 429, so cheks out
#saveRDS to pass on to add names
saveRDS(mRNA_full, "tcga_no_weird_full.RDS")

#split back to pheno and tcga_data
pheno <- mRNA_full[, 1:69]
tcga_data <- mRNA_full[, 70:60729]

saveRDS(tcga_data, "tcga_no_weird_counts.RDS")
saveRDS(pheno, "tcga_no_weird_pheno_XENA_TCGA.RDS")
