#Step8 WGCNA to confirm
setwd("~/rprojects/TCGA-OV-data") #wsl
library(WGCNA)
library(gridExtra)
library(glmnet)
library(ggplot2)
library(CorLevelPlot)
library(venn)
library(tidyverse)
allowWGCNAThreads() 

tcga_counts <- readRDS("tcga_counts_train.RDS")
pheno_train <- readRDS("tcga_pheno_train.RDS")


tcga_counts <- as.data.frame(t(tcga_counts))
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
#stf runs quite slow!
sft <- pickSoftThreshold(tcga_counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)
#ziurim max r square ir minimum mean conectivity

sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2) #library gridExtra
#choose above red line (above 0.8) and low conectivity (low point on the 2nd plot)

# convert matrix to numeric
tcga_counts[] <- sapply(tcga_counts, as.numeric) #reik šitu dalyku po to esančiai eilutei 

soft_power <- 6
temp_cor <- cor #cor function pasirinkti kitaip namespace error bus
cor <- WGCNA::cor #naudosim sita


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(tcga_counts,
                          maxBlockSize = 14000, #for fast clustering: how many genes in a block, visi genai viename bloke
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

saveRDS(bwnet, "bwnet.mRNA.TCGA-OV.RDS") 
#a save is necessary because bwnet can run for hours and then quit because of memory issues
cor <- temp_cor #sugražinam

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs #MEs yra module egenges


# Print out a preview
head(module_eigengenes) #clusters are colors


# get number of genes for each module
table(bwnet$colors) 

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


#traits 
# create traits file - binarize categorical variables
traits <- pheno_train %>% 
  mutate(vital_status_bin = ifelse(grepl('Dead', vital_status), 1, 0)) %>% 
  dplyr::select(70) #select the new colum, after the last one

pheno_train$figo_stage_f = gsub("[ABC]$", "", pheno_train$figo_stage)
table(pheno_train$figo_stage_f, useNA = "a") #stage na =2
pheno_train$figo_stage_f[pheno_train$figo_stage_f =="Stage I"] <- NA
#for figo stage, padaryti palyginimą tarp visų galimų stadijų
severity.out <- binarizeCategoricalColumns(pheno_train$figo_stage_f,
                                           includePairwise = TRUE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)
traits <- cbind(traits, severity.out)

pheno_train$neoplasmhistologicgrade_F = pheno_train$neoplasmhistologicgrade
table(pheno_train$neoplasmhistologicgrade_F, useNA = "a") #i want g2 vs G3
pheno_train$neoplasmhistologicgrade_F[pheno_train$neoplasmhistologicgrade_F %in% c("G4", "GB", "GX")] <- NA
severity.out2 <- binarizeCategoricalColumns(pheno_train$neoplasmhistologicgrade_F,
                                           includePairwise = TRUE,
                                           includeLevelVsAll = F,
                                           minCount = 1)
traits <- cbind(traits, severity.out2)

# Define numbers of genes and samples
nSamples <- nrow(tcga_counts)
nGenes <- ncol(tcga_counts)
#find corelation between eigengenes
module.trait.corr <- cor(module_eigengenes, traits, use = 'p') #pearson cor between eigengenes and traits
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples) #cor p values


heatmap.data <- merge(module_eigengenes, traits, by = 'row.names') #combine  to data frame

colnames(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names') #first colum negali buti numeric todel convertuojam i rownames


CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[113:119], #CLINICAL
             y = names(heatmap.data)[1:112], #"COLORS"
             col = c("blue1", "skyblue", "white", "pink", "red"))
#maybe open separately 20X45inch portrait as there are a 101 eigengenes :)
#significance denoted by *

module.gene.mapping <- as.data.frame(bwnet$colors) #modules stored as colors

mgmstage <- module.gene.mapping %>% 
  filter(`bwnet$colors` == "blue") %>% 
  rownames() #filtruosim genus kurie yra musu norimos "spalvos" 

mgm_grade <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'greenyellow') %>% 
  rownames() #filtruosim genus kurie yra musu norimos "spalvos" 

#Identifying driver genes
module.membership.measure <- cor(module_eigengenes, tcga_counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

#grade
gene.signf.corr.grade <- cor(tcga_counts, traits$data.G3.vs.G2, use = 'p')
gene.signf.corr.pvals.grade <- corPvalueStudent(gene.signf.corr.grade, nSamples)
ensembl.ids.grade <- gene.signf.corr.pvals.grade %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25) %>% 
  rownames()
ensembl.ids.grade


grade_glm <- relevant_genes_grade #is step7
intersect(ensembl.ids.grade, relevant_genes_grade) #12 intersect
#"MX2" "MTMR11" "SPRED3" "LNPK" "SGCB" "PCDHB4" "KNL1" "LGR4" "CEP152" "RBMS1" "KLK14" "TMX4" "SUV39H1"
#okay why do all the fancy clustering jei po to tiesiog top genes nusifiltruoji neatsizvelgiant?
intersect(ensembl.ids.grade, mgmGRADE) #"MX2"  "OAS3" "OASL"

#stage
gene.signf.corr.stage <- cor(tcga_counts, traits$`data.Stage IV.vs.all`, use = 'p')
gene.signf.corr.pvals.stage <- corPvalueStudent(gene.signf.corr.stage, nSamples)
ensembl.ids.stage <- gene.signf.corr.pvals.stage %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25) %>% 
  rownames()
ensembl.ids.stage

intersect(ensembl.ids.stage, mgmSTAGE4) #none

########################
mgmstage
stage_counts <- tcga_counts[colnames(tcga_counts) %in% mgmstage]

pheno_train$figo_stage
pheno_train$stage_num <- gsub("[ABC]$", "", pheno_train$figo_stage)
table(pheno_train$stage_num, useNA = "a") #stage na =2
pheno_train$stage_num[pheno_train$stage_num =="Stage I"] <- NA

stage <- pheno_train[pheno_train$stage_num %in% c("Stage II", "Stage III", "Stage IV"), ]
dim(stage) #lieka 334 zmones
#na trukdys
stage_counts <- stage_counts[ (rownames(stage_counts) %in% rownames(stage)), ]
stage_counts <- data.matrix(stage_counts)
dim(stage_counts)
grade_factor <- recode(stage$stage_num, "Stage II" = -1, "Stage III" = 0, "Stage IV" = 1)
stage_glm = cv.glmnet(
  x = stage_counts,
  y = grade_factor,
  alpha = 1, 
  family = "gaussian"
)
stage_glm #21 
stage_coef= coef(stage_glm, s="lambda.min") # the "coef" function returns a sparse matrix
stage_coef = stage_coef[stage_coef[,1] != 0,] 
stage_coef = stage_coef[-1]
relevant_genes_stage= names(stage_coef) # get names of the (non-zero) variables.
relevant_genes_stage  #

stage_genes_no_cluster <- c("GRHL3"  ,  "NRDC"  ,   "TMEM59"  , "CREG1" ,  
                            "SLC1A4"  , "PRADC1"  , "PDK1"   ,  "PRKRA"   ,
                            "CREB1",   "RAB43" ,   "TRIM59"  , "KLHL8"  ,
                            "HIST1H1E" ,"SBDS" ,    "TSPAN12" , "OSTF1" , 
                            "MAPKAP1" , "RASGEF1A","LGR4" ,    "TMEM123" , 
                            "PDGFD" ,   "BACE1" ,   "LRMP" ,    "ABHD4" , 
                            "ZNF821" ,  "SRR" ,     "PPP4R1" , "C18orf32",
                            "ZNF236" ,  "CTDP1" ,   "SLC44A2" , "CEBPG" , 
                            "CNFN",     "SULT2B1",  "RDH13" ,   "TRIB3"  ,
                            "FAM210B" , "PCBP3" ,   "ZNF74" ,   "SLC10A3",
                            "MT-ND5" )
intersect(stage_genes_no_cluster, relevant_genes_stage)
stage_intersect_list <- list(cluster_blue =relevant_genes_stage, no_clustering= stage_genes_no_cluster )
venn(stage_intersect_list)

