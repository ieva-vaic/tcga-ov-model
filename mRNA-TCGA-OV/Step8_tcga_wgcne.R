#Step8 WGCNA to confirm
setwd("~/rprojects/TCGA-OV-data") #wsl
library(WGCNA)
library(gridExtra)
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