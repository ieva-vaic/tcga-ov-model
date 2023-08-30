# STEP 0.0.1. looking at microRNA data
library("tidyverse")
library("SummarizedExperiment")
library("TCGAbiolinks")
library(gridExtra)
library(WGCNA)
#1.counts
#2.rpms
#3.after low count filtering, log-transformed rpms 
## (this: as.matrix(log2(tcga_mir_data.keep+0.5)))
setwd("C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-MiR data/RPMS_miR/")
miRcounts <- readRDS("tcga_mir_data.counts.RDS")
mirRMPs <- readRDS("tcga_mir_data.rpms.RDS")
miR_logs <- readRDS("tcga_mir_data.log.RDS")

miRcounts <- as.data.frame(miRcounts) # 1881, sveiki skaičiai 
mirRMPs <- as.data.frame(mirRMPs) # 1881 nuo labai mažu iki dideliu 
miR_logs <- as.data.frame(miR_logs) #Šitas pafiltruotas 1338 mirs log transformed tai su kableliais, dalis minusiniai

#išsirinksiu miR-200c, nes ji yra laikoma svarbia KV
"hsa-mir-200c" %in% rownames(miRcounts)

mir_200c_counts <- miRcounts %>% 
  filter(rownames(miRcounts) == "hsa-mir-200c") %>%
  t() %>%
  as.data.frame() 
mir_200c_counts$`hsa-mir-200c-counts` <- as.numeric(mir_200c_counts$`hsa-mir-200c`)

mir_200c_rpms <- mirRMPs %>% 
  filter(rownames(mirRMPs) == "hsa-mir-200c") %>%
  t() %>%
  as.data.frame() 
mir_200c_rpms$`hsa-mir-200c-rpms` <- as.numeric(mir_200c_rpms$`hsa-mir-200c`)

mir_200c_logrmps <- miR_logs %>% 
  filter(rownames(miR_logs) == "hsa-mir-200c") %>%
  t() %>%
  as.data.frame() 
mir_200c_logrmps$`hsa-mir-200c-logrpms` <- as.numeric(mir_200c_logrmps$`hsa-mir-200c`)

#histograms
logrpms <- ggplot(mir_200c_logrmps, aes(x=`hsa-mir-200c-logrpms` )) + 
  geom_histogram()
rpms <- ggplot(mir_200c_rpms, aes(x=`hsa-mir-200c-rpms` )) + 
  geom_histogram()
counts <- ggplot(mir_200c_counts, aes(x=`hsa-mir-200c-counts` )) + 
  geom_histogram()
grid.arrange(logrpms, rpms, counts)

##merged yra padarytas su log-transformed data
full_data_mirs <- readRDS("MERGED.rpms.RDS")

rownames(full_data_mirs) <- full_data_mirs$bcr_patient_barcode

dim(full_data_mirs) #489 1405
colnames(full_data_mirs[1338:1405]) 
pheno <- full_data_mirs[1339:1405] #atsiskiriu sutvarkyta pheno
rownames(pheno) <- pheno$bcr_patient_barcode
counts <- full_data_mirs[1:1338] #atsiskiriu counts, bet reik grazint t
counts <- as.data.frame(t(counts))
tail(rownames(counts)) #pameciau colnames bet tuoj atstatysiu
str(counts) #it is chr
LOGcounts <- counts


gsg <- goodSamplesGenes(t(LOGcounts)) #finds outliers rows samples, colums genes, t for transposing
summary(gsg) #weather there are outliers
gsg$allOK #fals tells there are outliers

table(gsg$goodGenes) #false is the number of outliers
table(gsg$goodSamples)

# remove genes that are detectd as outliers
counts <- LOGcounts[gsg$goodGenes == TRUE,]
str(counts) #liko 1561 is 1881 genu
htree <- hclust(dist(t(counts)), method = "average")
plot(htree) #distant samples should be excluded 


samples.to.be.excluded <- c("TCGA-23-1116", "TCGA-13-0730")
data.subset <- counts[,!(colnames(counts) %in% samples.to.be.excluded)]

#data.subset <- counts #i do this because I exluded no samples
# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset
# exclude outlier samples
colData <- pheno %>% filter(!row.names(.) %in% samples.to.be.excluded)
names(colData) 
making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))  



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
#biski wrid bet not as bad kaip buvo, renkuosi 14

counts[] <- sapply(counts, as.numeric) #reik ?itu dalyku po to esan?iai eilutei 

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


cor <- temp_cor #sugra?inam
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

#############################################################################

