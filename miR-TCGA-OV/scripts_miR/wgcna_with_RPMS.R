#this script is trying to do deseq2 normalization 
#and outlier selection from the WGCNA workfolw
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(CorLevelPlot) #this is from github, not in conda
#install.packages("remotes")
#remotes::install_github("kevinblighe/CorLevelPlot") #if running for the first time
library(gridExtra)

full_data_mirs <- readRDS("MERGED.rpms.RDS") #if not using snakemake
dim(full_data_mirs)
colnames(full_data_mirs[1338:1405]) 
pheno <- full_data_mirs[1339:1405] #atsiskiriu sutvarkyta pheno
rownames(pheno) <- pheno$bcr_patient_barcode
counts <- full_data_mirs[1:1339] #atsiskiriu counts, bet reik grazint t
counts <- as.data.frame(t(counts))
tail(rownames(counts)) #pameciau colnames bet tuoj atstatysiu
colnames(counts) <- counts[1339,]
counts <- counts[-1339,] #paliekam tik skaièius aplink
str(counts) #it is chr
counts[,1:489] <- sapply(counts[,1:489], as.numeric)
str(counts) #now it is numeric


#èia ðitas basically padaro mano counts matrica tai useless 
#counts$mir_names <- rownames(counts)
# data <- counts %>% 
#   gather(key = 'samples', value = 'counts', -mir_names) %>% 
#  # mutate(samples = gsub('\\.', '-', samples)) %>% 
#   inner_join(., pheno, by = c('samples' = 'full_barcode')) %>% 
#   select(1,2,3) %>% #4 pakeiciau i 3
#   spread(key = 'samples', value = 'counts') %>% 
#   column_to_rownames(var = 'mir_names') 

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(LOGcounts)) #finds outliers rows samples, colums genes, t for transposing
summary(gsg) #weather there are outliers
gsg$allOK #fals tells there are outliers

table(gsg$goodGenes) #false is the number of outliers
table(gsg$goodSamples)

# remove genes that are detectd as outliers
counts <- LOGcounts[gsg$goodGenes == TRUE,]
str(counts) #liko 1561 is 1881 genu

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(counts)), method = "average")
plot(htree) #distant samples should be excluded 

# pca - method 2

pca <- prcomp(t(counts))
pca.dat <- pca$x #x is information abouth principal components

pca.var <- pca$sdev^2 #variance explained by principal  components
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2) #percentage of the variance explained

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
#I don't wanna loose samples right now so no filtering by name
# exclude outlier samples
samples.to.be.excluded <- c('TCGA-59-A5PD', 'TCGA-61-1724', 'TCGA-61-2088',
"TCGA-59-2363", "TCGA-25-1329", "TCGA-24-1419", "TCGA-59-2348",
"TCGA-09-2044", "TCGA-61-1918", "TCGA-09-2054", 
"TCGA-24-2290", "TCGA-04-1652", 
"TCGA-30-1857", "TCGA-61-1721")
data.subset <- counts[,!(colnames(counts) %in% samples.to.be.excluded)]

#data.subset <- counts #i do this because I exluded no samples
# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset
# exclude outlier samples
colData <- pheno %>% filter(!row.names(.) %in% samples.to.be.excluded)

# fixing column names in colData
#colData <- pheno #i do this because I exlucded no samples
names(colData) #i have no problems thus no fixing
#names(colData) <- gsub(':ch1', '', names(colData)) #remove "ch1
#names(colData) <- gsub('\\s', '_', names(colData)) #remove spaces!

#add colData names because it has none
rownames(colData) <- colData$bcr_patient_barcode

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))  

# not creating dds as log counts are already normalized
# dds <- DESeqDataSetFromMatrix(countData = data.subset,
#                               colData = colData,
#                               design = ~ 1) # not spcifying model
# ## remove all genes with counts < 15 in more than 75% of samples 
# #èia reikia paskaièiuot 2 skaièiø!!!!!!! (489samples*0.75=XX)
# nrow(colData) * 0.75 #todel 2 sk yra 367
# ## suggested by WGCNA on RNAseq FAQ
# 
# dds75 <- dds[rowSums(counts(dds) >= 15) >= 366,]
# nrow(dds75) # 229 genes
# 
# # perform variance stabilization
# #dds_norm <- vst(dds75) #neveikia, NES MAÞAI GENØ todël:
# dds_norm <- varianceStabilizingTransformation(dds75, fitType = "local") #chek fittype
# # get normalized counts
# norm.counts <- assay(dds_norm) %>% 
#   t()

# 4. Network Construction 
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2)) #create a list of powers

# Call the network topology analysis function
sft <- pickSoftThreshold(counts,
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
#siuo atveju labai keistai pasiskirto topology modelis,
#pasirinksiu todel didziausia skaièiø arèiausiai kreives
#choose above red line (above 0.8) and low conectivity (low point on the 2nd plot) 
#I'll choose 22 on this

# convert matrix to numeric
counts[] <- sapply(counts, as.numeric) #reik ðitu dalyku po to esanèiai eilutei 

soft_power <- 14
temp_cor <- cor #cor function pasirinkti kitaip namespace error bus
cor <- WGCNA::cor #naudosim sita

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(counts,
                          maxBlockSize = 300, #for fast clustering: how many genes in a block, visi genai viename bloke
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,#buvo 0.25
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3) #buvo 3


cor <- temp_cor #sugraþinam

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

#everything is turqoise?
