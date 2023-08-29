#in this script I will do some descriptive statistics on mRNA
#so that I could do the same on the mir data as well
setwd("C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/nesutvarkyti failai TCGA-OV analysis/RDS files/")
#raw counts?
library(biomaRt)
library(tidyverse)
library(gridExtra)
library("TCGAbiolinks")
library(SummarizedExperiment)
library("GDCRNATools") #BUVO INSTALINTA DARBO KOMPE TIK 

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

be_nieko <- grid.arrange(arid1A, notch, betaC, Gapdh, hPRT1, Mzb1)
be_nieko
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
ensembl_names3
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

tcga_normalizavimas <- grid.arrange(arid1A, notch, betaC, Gapdh, hPRT1, Mzb1)
tcga_normalizavimas
###############################################################################
#what if I normalize the data another way?
count_matrix <- assay(tcga_data) #started with 60660 genes
#1 type 2: from GDCRNAtools: gdcVoomNormalization: does TNB ir voom transformatio
trainrnaExpr <- gdcVoomNormalization(counts = count_matrix, filter = FALSE) #60660 genes
#pridedam normal vardus
trainrnaExpr1 <- trainrnaExpr
trainrnaExpr1 <- as.data.frame(trainrnaExpr1)
trainrnaExpr1$ensembl_gene_id  <- rownames(trainrnaExpr1)
trainrnaExpr1$ensembl_gene_id <- gsub("\\..*", "", trainrnaExpr1$ensembl_gene_id )

ensembl_namesLIMMA <- trainrnaExpr1$ensembl_gene_id
ensembl_namesLIMMA <- as.data.frame(ensembl_namesLIMMA)

ensembl_namesLIMMA <- getBM(attributes = c('ensembl_gene_id','external_gene_name'), #values to retreve
                        filters = "ensembl_gene_id", #input on quary
                        values = ensembl_namesLIMMA, #idet dataframa
                        mart = ensembl.con) #sukurtas conection object
ensembl_namesLIMMA
###
counts_limma <- left_join(trainrnaExpr1, ensembl_namesLIMMA, by = "ensembl_gene_id") #34834   

#dabar in counts turim vardus uždėtus genų
#ARID1A
arid1a <- counts_limma %>% 
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
hprt <- counts_limma %>% 
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
gapdh <- counts_limma %>% 
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
betacat <- counts_limma %>% 
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
notch1 <- counts_limma %>% 
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
MZB <- counts_limma %>% 
  filter(external_gene_name == "MZB1") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("MZB1" = "V1") 
MZB$MZB1 <- as.numeric(MZB$MZB1)

Mzb1 <- ggplot(MZB, aes(x=MZB1)) + 
  geom_histogram()

limma_plot <- grid.arrange(arid1A, notch, betaC, Gapdh, hPRT1, Mzb1)
limma_plot 
##############################################################################
#what if I normalize the data another way?
count_matrix <- assay(tcga_data) #started with 60660 genes
coldata <-  readRDS("pheno.RDS")
#1 type 3: from DESEQ2

# create dds
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ 1) # not spcifying model

## remove all genes with counts < 15 in more than 75% of samples (429*0.75=321)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 321,] #321
nrow(dds75) # 17824 genes


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) 
#%>%t() #šiam kartui nereik

#pridedam normal vardus
norm.counts <- norm.counts
norm.counts <- as.data.frame(norm.counts)
norm.counts$ensembl_gene_id  <- rownames(norm.counts)
norm.counts$ensembl_gene_id <- gsub("\\..*", "", norm.counts$ensembl_gene_id )

ensembl_namesDESEQ <- norm.counts$ensembl_gene_id
ensembl_namesDESEQ <- as.data.frame(ensembl_namesDESEQ)

ensembl_namesDESEQ <- getBM(attributes = c('ensembl_gene_id','external_gene_name'), #values to retreve
                            filters = "ensembl_gene_id", #input on quary
                            values = ensembl_namesDESEQ, #idet dataframa
                            mart = ensembl.con) #sukurtas conection object
ensembl_namesDESEQ
###
counts_deseq <- left_join(norm.counts, ensembl_namesDESEQ, by = "ensembl_gene_id") #17824      

#dabar in counts turim vardus uždėtus genų
#ARID1A
arid1a <- counts_deseq %>% 
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
hprt <- counts_deseq %>% 
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
gapdh <- counts_deseq %>% 
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
betacat <- counts_deseq %>% 
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
notch1 <- counts_deseq %>% 
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
MZB <- counts_deseq %>% 
  filter(external_gene_name == "MZB1") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("MZB1" = "V1") 
MZB$MZB1 <- as.numeric(MZB$MZB1)

Mzb1 <- ggplot(MZB, aes(x=MZB1)) + 
  geom_histogram()

deseq_plot <- grid.arrange(arid1A, notch, betaC, Gapdh, hPRT1, Mzb1)
deseq_plot 

##########################
#DESEQ no WST
dds_simple <- assay(dds75)
#pridedam normal vardus
dds_simple <- as.data.frame(dds_simple)
dds_simple$ensembl_gene_id  <- rownames(dds_simple)
dds_simple$ensembl_gene_id <- gsub("\\..*", "", dds_simple$ensembl_gene_id )

ensembl_namesDESEQ_simple <- dds_simple$ensembl_gene_id
ensembl_namesDESEQ_simple <- as.data.frame(ensembl_namesDESEQ_simple)

ensembl_namesDESEQ_simple <- getBM(attributes = c('ensembl_gene_id','external_gene_name'), #values to retreve
                            filters = "ensembl_gene_id", #input on quary
                            values = ensembl_namesDESEQ_simple, #idet dataframa
                            mart = ensembl.con) #sukurtas conection object
ensembl_namesDESEQ_simple
###
counts_deseq_simple <- left_join(dds_simple, ensembl_namesDESEQ_simple, by = "ensembl_gene_id") #17824      

#dabar in counts turim vardus uždėtus genų
#ARID1A
arid1a <- counts_deseq_simple %>% 
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
hprt <- counts_deseq_simple %>% 
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
gapdh <- counts_deseq_simple %>% 
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
betacat <- counts_deseq_simple %>% 
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
notch1 <- counts_deseq_simple %>% 
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
MZB <- counts_deseq_simple %>% 
  filter(external_gene_name == "MZB1") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("MZB1" = "V1") 
MZB$MZB1 <- as.numeric(MZB$MZB1)

Mzb1 <- ggplot(MZB, aes(x=MZB1)) + 
  geom_histogram()

deseq_simple_plot <- grid.arrange(arid1A, notch, betaC, Gapdh, hPRT1, Mzb1)
deseq_simple_plot 

###############################################################################
#TCGA norm vs deseq without any model (with vst) - my most used normalizations
# I choose ARID1A for the fist plot
arid1a_deseq <- counts_deseq %>% 
  filter(external_gene_name == "ARID1A") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("ARID1A" = "V1") 
arid1a_deseq_ARID1A <- as.numeric(arid1a_deseq$ARID1A)

arid1a_biolinks <- counts_biolinks %>% 
  filter(external_gene_name == "ARID1A") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("ARID1A" = "V1") 
arid1a_biolinks_ARID1A <- as.numeric(arid1a_biolinks$ARID1A)


arid1as <- data.frame(arid1a_deseq_ARID1A, arid1a_biolinks_ARID1A)
#thankfully borth have the same amount of samples


# Basic scatter plot 
ggplot(arid1as, aes(x=arid1a_deseq_ARID1A, y=arid1a_biolinks_ARID1A)) + geom_point()


#what wst does?
arid1a_deseq_no_wst <- counts_deseq_simple %>% 
  filter(external_gene_name == "ARID1A") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("ARID1A" = "V1") 
arid1a_deseq_no_wst <- as.numeric(arid1a_deseq_no_wst$ARID1A)
arid1as2 <- data.frame(arid1a_deseq_ARID1A, arid1a_deseq_no_wst)
ggplot(arid1as2, aes(x=arid1a_deseq_ARID1A, y=arid1a_deseq_no_wst)) + geom_point()

#limma vs wst_deseq
#what wst does?
arid1a_limma <- counts_limma %>% 
  filter(external_gene_name == "ARID1A") %>%
  t() %>%
  as.data.frame() %>%
  filter(!row_number() %in% c(431, 432, 430)) %>%
  rename("ARID1A" = "V1") 
arid1a_limma <- as.numeric(arid1a_limma$ARID1A)
arid1as3 <- data.frame(arid1a_deseq_ARID1A, arid1a_limma)
ggplot(arid1as3, aes(x=arid1a_deseq_ARID1A, y=arid1a_limma)) + geom_point()

#biolinks_vs_limma
arid1as4 <- data.frame(arid1a_biolinks_ARID1A, arid1a_limma)
ggplot(arid1as4, aes(x=arid1a_biolinks_ARID1A, y=arid1a_limma)) + geom_point()


