#Step 3: get tcga count data -> outlier chek with hplots -> lognormalize
#steps 1 (tcga data download) and 2 (pheno parse) is done in the gtex/tcga folder
library(tidyverse)
library(dendextend)
library(WGCNA)
library(GDCRNATools)
library(SummarizedExperiment)
#get count matrix, rows as genes
#working with full data! previously I was doing this step after test/train split,
#now doing this so not to repeat the step and to stay consistent with gtex analysis
coldata_tcga <- readRDS("C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/pheno_no_empty_data.RDS")
rownames(coldata_tcga) <- coldata_tcga$barcode
tcga_data <- readRDS("C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/tcga_data.RDS")
#make a dataframe 
mRNA_counts <- assay(tcga_data)
mRNA_counts <- as.data.frame(mRNA_counts) #60660 transcripts (rows), 429 zmones (cols) 
mRNA_counts[,1:429] <- lapply(mRNA_counts[,1:429], as.numeric)
str(mRNA_counts) #60660 genes
###################################################################################
#OUTLIERS

# detect outlier genes with gsg
gsg <- goodSamplesGenes(mRNA_counts) #no transpose, wgcna package
summary(gsg) #weather there are outliers
gsg$allOK #fals tells there are outliers

table(gsg$goodGenes) #false is the number of outliers
table(gsg$goodSamples) #no outliers


# detect outlier samples - hierarchical clustering 

htree_mrna <- hclust(dist(t(mRNA_counts)), method = "average") #can take some time 

pdf(file="C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/Figures mRNA/htree_tcga_only.pdf", height=30, width=60) #i?saugijimui didesniu formatu
plot(htree_mrna) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off()

#take out the outliers
# exclude outlier samples
samples.to.be.excluded <- c("TCGA−13−0919−01A−01R−1564−13") #maybe this one is unussually high
mRNA_counts <- mRNA_counts[,!(colnames(mRNA_counts) %in% samples.to.be.excluded)] #lieka 429 zmones bendrai

# exclude outlier samples
coldata_tcga <- coldata_tcga %>% 
  filter(!row.names(.) %in% samples.to.be.excluded) #lieka 429 zmones

all(rownames(coldata_tcga) %in% colnames(mRNA_counts))
nrow(coldata_tcga) == ncol(mRNA_counts)
all(rownames(coldata_tcga) == colnames(mRNA_counts)) #eile nesutampa?

target <- colnames(mRNA_counts)
coldata_tcga <- coldata_tcga[match(target, rownames(coldata_tcga)),]
all(rownames(coldata_tcga) == colnames(mRNA_counts)) #dabar sutampa
rm(target)

saveRDS(coldata_tcga, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/coldata_tcga_mrna.RDS")
###############################################################################
# Normalize, I will use GDC RNA tools
#1 type 2: from GDCRNAtools: gdcVoomNormalization: does TNB ir voom transformatio
tcga_mRNA_voom <- gdcVoomNormalization(counts = mRNA_counts, filter = FALSE) #
#pridedam normal vardus
tcga_mRNA_voomt <- t(tcga_mRNA_voom)

saveRDS(tcga_mRNA_voom, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/tcga_mRNA_voom.RDS")
###############################################################################

