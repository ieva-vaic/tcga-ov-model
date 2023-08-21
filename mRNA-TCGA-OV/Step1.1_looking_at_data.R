#in this script I will do some descriptive statistics on mRNA
#so that I could do the same on the mir data as well
setwd("C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data")
#raw counts?
library(biomaRt)
library(tidyverse)
library(gridExtra)
library("TCGAbiolinks")
library("GDCRNAtools") #BUVO INSTALINTA DARBO KOMPE TIK 

library(devtools)
devtools::install_github(repo='Jialab-UCR/GDCRNATools')
library(GDCRNATools)

tcga_data <- readRDS("tcga_data.RDS")
counts <- as.data.frame(assay(tcga_data))
#Can i change ensemble ids to gene names?
ensembl_names <- rownames(counts)
ensembl_names2 <- gsub("\\..*", "", ensembl_names)
ensembl_names2 <- as.data.frame(ensembl_names2)
# cia uztenka karta padaryt
listEnsembl() #shows available databases
ensembl <- useEnsembl(biomart = "genes") #unsuportive gali but
datasets <- listDatasets(ensembl) #issirinkti human genes
ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl') #name of data base ir data set
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)
###
gene_names <- getBM(attributes = c('ensembl_gene_id','external_gene_name'), #values to retreve
                    filters = "ensembl_gene_id", #input on quary
                    values = ensembl_names2, #idet dataframa
                    mart = ensembl.con) #sukurtas conection object
gene_names
###
counts$ensembl_names <- rownames(counts)
counts$ensembl_gene_id <-  gsub("\\..*", "", counts$ensembl_names)
counts <- left_join(counts, gene_names, by = "ensembl_gene_id")


#dabar in counts turim vardus uždėtus genų
#ARID1A
arid1a <- counts %>% 
  filter(external_gene_name == "ARID1A") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("ARID1A" = "V1") 
arid1a$ARID1A <- as.numeric(arid1a$ARID1A)

arid1A <- ggplot(arid1a, aes(x=ARID1A)) + 
  geom_histogram()
#HPRT
#dabar in counts turim vardus uždėtus genų
hprt <- counts %>% 
  filter(external_gene_name == "HPRT1") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("HPRT1" = "V1") 
hprt$HPRT1 <- as.numeric(hprt$HPRT1)

hPRT1 <- ggplot(hprt, aes(x=HPRT1)) + 
  geom_histogram()

#GAPDH
#dabar in counts turim vardus uždėtus genų
gapdh <- counts %>% 
  filter(external_gene_name == "GAPDH") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("GAPDH" = "V1") 
gapdh$GAPDH <- as.numeric(gapdh$GAPDH)

Gapdh <- ggplot(gapdh, aes(x=GAPDH)) + 
  geom_histogram()

#CTNNB1
#dabar in counts turim vardus uždėtus genų
betacat <- counts %>% 
  filter(external_gene_name == "CTNNB1") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("CTNNB1" = "V1") 
betacat$CTNNB1 <- as.numeric(betacat$CTNNB1)

betaC <- ggplot(betacat, aes(x=CTNNB1)) + 
  geom_histogram()

#NOTCH1
#dabar in counts turim vardus uždėtus genų
notch1 <- counts %>% 
  filter(external_gene_name == "NOTCH1") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("NOTCH1" = "V1") 
notch1$NOTCH1 <- as.numeric(notch1$NOTCH1)

notch <- ggplot(notch1, aes(x=NOTCH1)) + 
  geom_histogram()

#MZB1
#dabar in counts turim vardus uždėtus genų
MZB <- counts %>% 
  filter(external_gene_name == "MZB1") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("MZB1" = "V1") 
MZB$MZB1 <- as.numeric(MZB$MZB1)

Mzb1 <- ggplot(MZB, aes(x=MZB1)) + 
  geom_histogram()

grid.arrange(arid1A, notch, betaC, Gapdh, hPRT1, Mzb1)

###############################################################################
#what if I normalize the data?
count_matrix <- assay(tcga_data) #started with 60660 genes
#1 type 1: from TCGAbiolinks TCGAnormalize
train_dataNorm <- TCGAanalyze_Normalization(
  tabDF = count_matrix, 
  geneInfo =  geneInfoHT
) #(46462)
#clean low and high counts
# quantile filter of genes #34766 genes
train_dataNorm_filtered <- TCGAanalyze_Filtering(
  tabDF = train_dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
) #(34834)

train_dataNorm_filtered <- as.data.frame(train_dataNorm_filtered)
train_dataNorm_filtered$ensembl_gene_id  <- rownames(train_dataNorm_filtered)

ensembl_names3 <- rownames(train_dataNorm_filtered)
ensembl_names3 <- as.data.frame(ensembl_names3)
ensembl_names3 <- getBM(attributes = c('ensembl_gene_id','external_gene_name'), #values to retreve
                    filters = "ensembl_gene_id", #input on quary
                    values = ensembl_names3, #idet dataframa
                    mart = ensembl.con) #sukurtas conection object
gene_names
###
counts_biolinks <- left_join(train_dataNorm_filtered, ensembl_names3, by = "ensembl_gene_id") #34834   

#dabar in counts turim vardus uždėtus genų, tik iklijuosiu biolinks_counts
#ARID1A
arid1a <- counts_biolinks %>% 
  filter(external_gene_name == "ARID1A") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("ARID1A" = "V1") 
arid1a$ARID1A <- as.numeric(arid1a$ARID1A)

arid1A <- ggplot(arid1a, aes(x=ARID1A)) + 
  geom_histogram()
#HPRT
#dabar in counts turim vardus uždėtus genų
hprt <- counts_biolinks %>% 
  filter(external_gene_name == "HPRT1") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("HPRT1" = "V1") 
hprt$HPRT1 <- as.numeric(hprt$HPRT1)

hPRT1 <- ggplot(hprt, aes(x=HPRT1)) + 
  geom_histogram()

#GAPDH
#dabar in counts turim vardus uždėtus genų
gapdh <- counts_biolinks %>% 
  filter(external_gene_name == "GAPDH") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("GAPDH" = "V1") 
gapdh$GAPDH <- as.numeric(gapdh$GAPDH)

Gapdh <- ggplot(gapdh, aes(x=GAPDH)) + 
  geom_histogram()

#CTNNB1
#dabar in counts turim vardus uždėtus genų
betacat <- counts_biolinks %>% 
  filter(external_gene_name == "CTNNB1") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("CTNNB1" = "V1") 
betacat$CTNNB1 <- as.numeric(betacat$CTNNB1)

betaC <- ggplot(betacat, aes(x=CTNNB1)) + 
  geom_histogram()

#NOTCH1
#dabar in counts turim vardus uždėtus genų
notch1 <- counts_biolinks %>% 
  filter(external_gene_name == "NOTCH1") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("NOTCH1" = "V1") 
notch1$NOTCH1 <- as.numeric(notch1$NOTCH1)

notch <- ggplot(notch1, aes(x=NOTCH1)) + 
  geom_histogram()

#MZB1
#dabar in counts turim vardus uždėtus genų
MZB <- counts_biolinks %>% 
  filter(external_gene_name == "MZB1") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("MZB1" = "V1") 
MZB$MZB1 <- as.numeric(MZB$MZB1)

Mzb1 <- ggplot(MZB, aes(x=MZB1)) + 
  geom_histogram()

grid.arrange(arid1A, notch, betaC, Gapdh, hPRT1, Mzb1)

###############################################################################
#what if I normalize the data?
count_matrix <- assay(tcga_data) #started with 60660 genes
#1 type 2: from GDCRNAtools: gdcVoomNormalization: does TNB ir voom transformatio
trainrnaExpr <- gdcVoomNormalization(counts = count_matrix, filter = FALSE) #60660 genes
#okay I can't do it ant sito kompo daba