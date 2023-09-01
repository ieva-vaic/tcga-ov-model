#Step 1.0 remove weird people from tcga
library(tidyverse)
setwd("~/rprojects/TCGA-OV-data")
#po to i prieki nukelt
#splice with pheno
pheno_final <- readRDS("pheno_no_empty_data.RDS")
tcga_data <- readRDS("tcga_with_names.RDS")
#kol kas man nereikia names filtravimui, tiesiog toks failas yra, skip if using no name data
rownames(tcga_data) <- tcga_data$ensembl
tcga_data <- tcga_data[, -c(430, 431, 432, 433)]
tcga_counts_t <- t(tcga_data) #praradom colnames
tcga_counts_t <- as.data.frame(tcga_counts_t) #non-numeric
tcga_counts_t$barcode <- rownames(tcga_counts_t)

tgca_pheno <- right_join(pheno_final, tcga_counts_t, by = "barcode")
dim(tgca_pheno) #429 zmones ir 60510 (60464 genai + 47 clinical)
rownames(tgca_pheno) <- tgca_pheno$barcode

################################################################################
#split off weird cases
#split weird cases from the rest of the data by definition and treatments
mRNA_full1 <- tgca_pheno #temporarily change name
split_by_definition <- split(mRNA_full1, f = mRNA_full1$definition, drop = T)
weird_cases <- split_by_definition$`Recurrent Solid Tumor` # weird group 7 people
mRNA_full <- split_by_definition$`Primary solid Tumor` #new new phenodata

split_prior_treatment <- split(mRNA_full, f = mRNA_full$prior_treatment, drop = T)
mRNA_full <- split_prior_treatment$No #now the mRNA_full are good cases
weird_cases <- rbind(weird_cases, split_prior_treatment$Yes) 

dim(weird_cases) # 8 zmones nusišalina
dim(mRNA_full) #421 lieka iš 429, so cheks out
#saveRDS to pass on to add names
saveRDS(mRNA_full, "tcga_no_weird.RDS")

#get the names of the 8 removed, so I could filter the gtex later and not repeat everything
weird_cases$barcode
#i leave it like this, will be easier to just copy
c("TCGA-13-0913-02A-01R-1564-13", "TCGA-13-1489-02A-01R-1565-13", "TCGA-29-1705-02A-01R-1567-13",
 "TCGA-29-1707-02A-01R-1567-13", "TCGA-29-1770-02A-01R-1567-13", "TCGA-29-2414-02A-01R-1569-13",
 "TCGA-61-2008-02A-01R-1568-13", "TCGA-61-1721-01A-01R-1569-13")
