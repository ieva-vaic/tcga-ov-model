#Step 2: splice gtex and TCGA-OV mRNA counts together
library(tidyverse)
library("TCGAbiolinks")
library("SummarizedExperiment")
setwd("~/rprojects/TCGA-OV-data")
gtex_counts <- readRDS("GTEX/gtex_counts.RDS")

#the gtex is now a matrix that looks like this:
#rows: ENSEMBLIDS (transcipt level), all 56200 of them 
head(colnames(gtex_counts))
#name: ENSEMBL IDS
#Description: real names 
#gtex_counts[3:181] #is case names 

#make a numeric dataframe from gtex
gtex_df <- as.data.frame(gtex_counts)
rownames(gtex_df) <- gtex_df$Name
gtex_df <- gtex_df[, -c(1,2)]
gtex_df[,1:180] <- lapply(gtex_df[,1:180], as.numeric)
dim(gtex_df) #56200 genes

#read in tcga counts data for mRNA 
tcga_data <- readRDS("tcga_data.RDS")

##assay data converted to dataframe
mRNA_counts <- assay(tcga_data)
mRNA_counts <- as.data.frame(mRNA_counts) #60660 transcripts (rows), 429 zmones (cols) 
mRNA_counts[,1:429] <- lapply(mRNA_counts[,1:429], as.numeric)
str(mRNA_counts) #60660 genes

#filter rows to match the tcga transcript ids
gtex_names <- rownames(gtex_df)
filtered_mRNA <- filter(mRNA_counts, rownames(mRNA_counts) %in% gtex_names) #tidyverse command
dim(filtered_mRNA) #liko 35117 genai

###############################################################################
#splice with pheno
pheno_final <- readRDS("pheno_no_empty_data.RDS")
filtered_mRNA_t <- as.data.frame(t(filtered_mRNA))
filtered_mRNA_t$barcode <- rownames(filtered_mRNA_t)
tail(colnames(filtered_mRNA_t))
mRNA_full <- left_join(pheno_final, filtered_mRNA_t, by = "barcode")
dim(mRNA_full) #429 zmones ir 35164 (35118 genai + 47 clinical)
rownames(mRNA_full) <- mRNA_full$barcode

################################################################################
#split off weird cases
#split weird cases from the rest of the data by definition and treatments
mRNA_full1 <- mRNA_full #temporarily change name
split_by_definition <- split(mRNA_full1, f = mRNA_full1$definition, drop = T)
weird_cases <- split_by_definition$`Recurrent Solid Tumor` # weird group 7 people
mRNA_full <- split_by_definition$`Primary solid Tumor` #new new phenodata

split_prior_treatment <- split(mRNA_full, f = mRNA_full$prior_treatment, drop = T)
mRNA_full <- split_prior_treatment$No #now the mRNA_full are good cases
weird_cases <- rbind(weird_cases, split_prior_treatment$Yes) 

dim(weird_cases) # 8 zmones nusišalina
dim(mRNA_full) #421 lieka iš 429, so cheks out
##############################################################################
#finally splice all together

#tik counts sujungt
dim(filtered_mRNA) #primenu kad liko 35117 genai

mRNA_filtered_names <- rownames(filtered_mRNA)
filtered_gtex <- filter(gtex_df, rownames(gtex_df) %in% mRNA_filtered_names) #liko 35117 genai
dim(filtered_gtex) #35117   180

#transpose GTEX
gtex_t <- as.data.frame(t(filtered_gtex)) # 180 zmoniu 56200genu transcriptu
dim(mRNA_full) #421 (nes ner weirds) 35164 mrnr counts + pheno 47
dim(gtex_t) #180 35117 #gtex
dim(filtered_mRNA_t) #mRNA counts +"barcode"  su dar likusiais wierd

#create gtex column
mRNA_full$gtex <- "ovarian cancer"
gtex_t$gtex <- "control"

#dabar reikia prie kiekvieno geno stulpelio prideti naujas eilutes, zmones :rbind
gtcga <- bind_rows(mRNA_full, gtex_t)
dim(gtcga) #601(421 + 180) 35165(sutampa su mRNA_full+ gtex stulpelis)
saveRDS(gtcga, "joined_gtex_tcga.RDS")