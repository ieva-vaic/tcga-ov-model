#Step 4: get joined tcga/gtex count data -> outlier chek with hplots -> lognormalize
library(tidyverse)
library(dendextend)
library(WGCNA)
library(GDCRNATools)
#get count matrix, rows as genes
#working with full data!
setwd("~/rprojects/TCGA-OV-data") #wsl
gtcga <- readRDS("joined_gtex_tcga.RDS")
counts_gtcga <- readRDS("filtered_gtgca_counts.RDS") #large numeric matrix with rows as genes
coldata_gtcga <- readRDS("coldata_gtcga_full.RDS")

##############################################################################
#OUTLIERS

# detect outlier genes with gsg
gsg <- goodSamplesGenes(counts_gtcga) #no transpose, wgcna package
summary(gsg) #weather there are outliers
gsg$allOK #fals tells there are outliers

table(gsg$goodGenes) #false is the number of outliers
table(gsg$goodSamples) #no outliers


# detect outlier samples - hierarchical clustering 

htree <- hclust(dist(t(counts_gtcga)), method = "average") #can take some time 

pdf(file="figures/htree_no_pseudo_gtex_tcga.pdf", height=30, width=60) #išsaugijimui didesniu formatu
plot(htree) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off()

#add colors to hclust, (library dendextend)
dend <- as.dendrogram(htree)
col_aa_red <- ifelse(grepl("GTEX", labels(dend)), "red", "blue")
dend2 <- assign_values_to_leaves_edgePar(dend=dend, value = col_aa_red, edgePar = "col") #

pdf(file="figures/htree_no_pseudo_gtex_tcga_red.pdf", height=80, width=100) #išsaugijimui didesniu formatu
plot(dend2) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off()

#I think TCGA−13−0760−01A−01R−1564−13 should be removed but later


###############################################################################
# Normalize, I will use GDC RNA tools
#1 type 2: from GDCRNAtools: gdcVoomNormalization: does TNB ir voom transformatio
mRNA_voom <- gdcVoomNormalization(counts = counts_gtcga, filter = FALSE) #
#pridedam normal vardus
mRNA_voomt <- t(mRNA_voom)

saveRDS(mRNA_voom, "mrna_voom_no_pseudo.RDS")
###############################################################################
#see normalization 
htree_norm <- hclust(dist(mRNA_voomt), method = "average") #can take some time
dend3 <- as.dendrogram(htree_norm)
col_aa_red <- ifelse(grepl("GTEX", labels(dend3)), "red", "blue")
dend4 <- assign_values_to_leaves_edgePar(dend=dend3, value = col_aa_red, edgePar = "col") 
pdf(file="figures/htree_no_pseudo_gtex_tcga_norm_red.pdf", height=80, width=100) #išsaugijimui didesniu formatu
plot(dend4) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off()
#it is very bad, viskas idealiai atsiskyrę, su deseq irgi idealiai, o tcga nieko nepadaro lieka same kaip neclusterintas


