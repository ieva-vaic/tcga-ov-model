#Step 5 model stage
library(tidyverse)
library(glmnet)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler) #go
tcga_counts_train <- readRDS( "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/train_mRNA_tcga_normcounts.RDS")
tcga_counts_train <- t(tcga_counts_train)
colData <- readRDS("C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/train_tcga_coldata.RDS")

#new coding for the Stage
table(colData$figo_stage, useNA = "a") #2na, stae 1 nera?
colData$figo_stage_f = gsub("[ABC]$", "", colData$figo_stage)
table(colData$figo_stage_f, useNA = "a") #stage na =2

#remove na
colData <- filter(colData, !is.na(colData$figo_stage_f))#be tidyverse blogai bus

#new coding: stage II is -1, stage III is 0, stage 3 is +1
colData$figo_recode <- colData$figo_stage_f
str(colData$figo_recode)
colData$figo_recode <- ordered(colData$figo_recode, levels = c("Stage II", "Stage III", "Stage IV"))
table(colData$figo_recode)

#new coding with numbers
colData$figo_recode_num <- colData$figo_stage_f
colData$figo_recode_num <- recode(colData$figo_recode_num, "Stage II" = -1, "Stage III" = 0, "Stage IV" = 1)
colData$figo_recode_num

#now that the 2 NA and 1 stage I is deleted, make sure the counts are equal as well
all(rownames(colData) %in% colnames(tcga_counts_train))
nrow(colData) == ncol(tcga_counts_train)
stage_samples <- rownames(colData) #make a list of new samples
stage_data <- tcga_counts_train[, stage_samples] #subset by sample name
dim(stage_data) #60660   324
nrow(colData) == ncol(stage_data)
all(rownames(colData) == colnames(stage_data)) #all good


###############################################################################
#model #1
## GLMET No.1
# clinical feature: figo_recode (II>III>IV)
# model family gaussian -> erroras, todel multinomial

norm.data <- t(stage_data) #need apversto
train_response_recode <- colData$figo_recode
res_factor = cv.glmnet(
  x = norm.data,
  y = train_response_recode,
  alpha = 0.5, 
  family = "multinomial"
)
res_factor # nieko neselektina

#model #2
## GLMET No.2
# clinical feature: figo_recode_num (-1>0>1)
# model family gaussian

norm.data <- t(stage_data) #need apversto
train_response_recode_num <- colData$figo_recode_num
res_num = cv.glmnet(
  x = norm.data,
  y = train_response_recode_num,
  alpha = 0.5, 
  family = "gaussian"
)
res_num # 49 min, 0 1se

# Getting genes that contribute for the prediction
res_coef_recoded_num = coef(res_num, s="lambda.min") # the "coef" function returns a sparse matrix
head(res_coef_recoded_num) # in a sparse matrix the "." represents the value of zero
# get coefficients with non-zero values
res_coef_recoded_num = res_coef_recoded_num[res_coef_recoded_num[,1] != 0,] 
# note how performing this operation changed the type of the variable
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef_recoded_num = res_coef_recoded_num[-1]
relevant_genes_recoded_num = names(res_coef_recoded_num) # get names of the (non-zero) variables.
length(relevant_genes_recoded_num) # 49

#convert relavant gene names
relevant_genes_recoded_num <- gsub("\\..*", "",relevant_genes_recoded_num)
relevant_genes_recoded_num
#biomart
# cia uztenka karta padaryt
listEnsembl() #shows available databases
ensembl <- useEnsembl(biomart = "genes") #unsuportive gali but
datasets <- listDatasets(ensembl) #issirinkti human genes
ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl') #name of data base ir data set
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)
##
recorded.num.genes <- getBM(attributes = c('ensembl_gene_id','external_gene_name'), #values to retreve
                            filters = "ensembl_gene_id", #input on quary
                            values = relevant_genes_recoded_num,
                            mart = ensembl.con) #sukurtas conection object
recorded.num.genes
write(recorded.num.symbols1, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/Figures mRNA/stage_genes.txt")

#############################################################################
#GO for 49 genes
#from org.Hs.eg.db
recorded.num.symbols <- recorded.num.genes$external_gene_name
recorded.num.symbols1 <- recorded.num.symbols[-c(26, 28, 31:35, 38, 39, 42:49)] #no na's allowed
recorded.num.entrez <- mapIds(org.Hs.eg.db, recorded.num.symbols1, 'ENTREZID', 'SYMBOL')
#no na's allowed - this time na's are in go terms
#recorded.num.entrez <- recorded.num.entrez[!is.na(recorded.num.entrez)] #jei butu na
#names chr vector of entrez ids

ggo.recorded.num <- groupGO(gene = recorded.num.entrez,
                            OrgDb    = org.Hs.eg.db,
                            ont      = "CC",
                            level    = 3,
                            readable = TRUE)

head(ggo.recorded.num) # count is 0 because i provided no counts!
ggo.recorded_num_df <- as.data.frame(ggo.recorded.num)
ggo.filtered.recorded_num_df <- ggo.recorded_num_df %>% filter(Count != "0")
ggo_stages <- ggo.filtered.recorded_num_df
View(ggo_stages)
#ego = enrichment go
ego.recorded.num <- enrichGO(gene = recorded.num.entrez,
                             OrgDb         = org.Hs.eg.db,
                             ont           = "CC",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = TRUE)
as.data.frame(ego.recorded.num) #0 terms


# ensembl_gene_id external_gene_name
# 1  ENSG00000068615              REEP1 eceptor expression-enhancing protein 1 (REEP1), neural cell gene
# 2  ENSG00000078618               NRDC nardilysin convertase, pseudogene
# 3  ENSG00000101255              TRIB3 TRIB3 is a pseudokinase protein that acts as a negative feedback regulator of the ATF4-dependent transcription during the integrated stress response (ISR
# 4  ENSG00000105251                SHD pseudogene
# 5  ENSG00000105341              DMAC2 Distal Membrane Arm Assembly Component 2, is a protein-coding gene that plays a role in the assembly of the mitochondrial respiratory chain complex I
# 6  ENSG00000105427               CNFN cornifelin, no protein product
# 7  ENSG00000108688               CCL7 cytokine
# 8  ENSG00000119421             NDUFA8 NADH dehydrogenase [ubiquinone] 1 alpha subcomplex subunit 8, is a protein-coding gene that plays a role in mitochondrial function and energy production
# 9  ENSG00000129353            SLC44A2 SLC44A2 (Solute Carrier Family 44 Member 2) is a protein-coding gene that plays a role in choline transport
# 10 ENSG00000130962              PRRG1 PRRG1 (Proline Rich And Gla Domain 1) no known function
# 11 ENSG00000137673               MMP7 Matrix metalloproteinase 7 (MMP7), also known as matrilysin
# 12 ENSG00000143162              CREG1 cellular repressor of E1A stimulated genes 1
# 13 ENSG00000149573              MPZL2 transmembrane glycoprotein and a member of the immunoglobulin superfamily that mediates homophilic cell-cell adhesion
# 14 ENSG00000152558            TMEM123 transmembrane protein 123, unclear function
# 15 ENSG00000162998               FRZB SFRP3, is a protein-coding gene that encodes for frizzled-related protein, a Wnt-binding protein
# 16 ENSG00000174562              KLK15 (Kallikrein-Related Peptidase 15)
# 17 ENSG00000175600              SUGCT for succinyl-CoA:glutarate-CoA transferase, an enzyme that is involved in the metabolism of glutarate
# 18 ENSG00000175701               MTLN mitoregulin, is a gene that encodes for a long non-coding RNA (lncRNA)
# 19 ENSG00000180228              PRKRA protein kinase, interferon-inducible 
# 20 ENSG00000181374              CCL13 chemokine
# 21 ENSG00000185294             SPPL2C signal peptide peptidase like 2C, testis gene
# 22 ENSG00000186318              BACE1 beta-secretase 1, is a protein-coding gene that encodes for a transmembrane aspartyl protease
# 23 ENSG00000189320            FAM180A no function known
# 24 ENSG00000211979           IGHV7-81 Immunoglobulin Heavy Variable 7-81 (Non-Functional)
# 25 ENSG00000228041            EIF3JP1 pseudogene
# 26 ENSG00000232855                   
# 27 ENSG00000234030           TMEM97P1 pseudogene
# 28 ENSG00000234519                   
# 29 ENSG00000236430            KRT8P29 pseudogene
# 30 ENSG00000237708            SPIN2P1 pseudogene
# 31 ENSG00000237870                   
# 32 ENSG00000238906                   
# 33 ENSG00000240440                   
# 34 ENSG00000243129                   
# 35 ENSG00000250243                   
# 36 ENSG00000251301          LINC02384 lncRNA
# 37 ENSG00000256356            HSPA8P5 hsp70 pseudogene
# 38 ENSG00000256734                   
# 39 ENSG00000256984                   
# 40 ENSG00000258972           NDUFB8P1 pseudogene 
# 41 ENSG00000263407            MIR4660 microRNA
# 42 ENSG00000267922                   
# 43 ENSG00000268864                   
# 44 ENSG00000269807                   
# 45 ENSG00000275115                   
# 46 ENSG00000280167                   
# 47 ENSG00000282021                   
# 48 ENSG00000286085                   
# 49 ENSG00000286089           

################################################################################
#sutapimu gaudymas su recorded.num.symbols1
gtcga_names2 <- read.delim("C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/Figures mRNA/gtex_vs_tcga_genes.txt", header = T)
gtcga_names2 <- gtcga_names2$gtcga_genes_names
gtcga_names2

intersect(recorded.num.symbols1, gtcga_names2) #no interect
