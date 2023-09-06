#Step 4: get joined tcga/gtex count data -> outlier chek with hplots -> lognormalize
library(tidyverse)
library(dendextend)
library(WGCNA)
library(GDCRNATools)
library(edgeR)
library(limma)
library(DESeq2)
#get count matrix, rows as genes
#working with full data!
setwd("~/rprojects/TCGA-OV-data") #wsl

counts_gtcga <- readRDS("gtcga_final_counts.RDS") #large numeric df with rows as genes
##############################################################################
#OUTLIERS
counts_gtcga <- data.matrix(counts_gtcga) #19193   596
# detect outlier genes with gsg
gsg <- goodSamplesGenes(counts_gtcga) #no transpose, wgcna package
summary(gsg) #weather there are outliers
gsg$allOK #fals tells there are outliers

table(gsg$goodGenes) #false is the number of outliers
table(gsg$goodSamples) #no outliers


# detect outlier samples - hierarchical clustering 

htree <- hclust(dist(t(counts_gtcga)), method = "average") #can take some time 

pdf(file="figures/htree_proteins_gtex_tcga.pdf", height=30, width=60) #išsaugijimui didesniu formatu
plot(htree) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off()

##probably ta tolimiausia zmogu nusalint reiktu. 
#XENA tutotiale plius buvo nusalintas TCGA.25.1870.01
#TCGA−13−1499−01A−01R−1565−13 nusalinsiu
which(colnames(counts_gtcga)=="TCGA-13-1499-01A-01R-1565-13") #297
counts_gtcga <- counts_gtcga[, -297] #19193   595 liko 

#nusibraizysiu dendograma pagal gtex vs tcga
dend <- as.dendrogram(htree)
col_aa_red <- ifelse(grepl("GTEX", labels(dend)), "red", "blue")
dend2 <- assign_values_to_leaves_edgePar(dend=dend, value = col_aa_red, edgePar = "col") #

pdf(file="figures/htree_protein_gtex_tcga_red.pdf", height=80, width=100) #išsaugijimui didesniu formatu
plot(dend2) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off()
###############################################################################
# Normalize, I will use GDC RNA tools
#GDC RNA tools, with filter, model added (issitraukus GDC RNA tools funkcija)
counts <- counts_gtcga
expr = DGEList(counts = counts)
expr = calcNormFactors(expr)

snames = colnames(counts_gtcga);
group = substr(snames, 1, 4); #Sets up level information for samples.

# #unfiltered
# v_unfiltered <- voom(expr, design= model.matrix(~0 + group), plot = TRUE)$E

#filtered
keepALL <- rowSums(cpm(expr) > 1) >= 0.5*ncol(counts)
nGenes <- as.numeric(summary(keepALL)[2]) + 
  as.numeric(summary(keepALL)[3])
nKeep <- summary(keepALL)[3]

cat (paste('Total Number of genes: ', nGenes, '\n', sep=''))
cat (paste('Number of genes for downstream analysis: ', nKeep, 
           '\n', sep=''))

exprALL <- expr[keepALL,,keep.lib.sizes = TRUE]
v_filtered <- voom(exprALL, design= model.matrix(~0 + group), plot = TRUE)$E

saveRDS(v_filtered, "mrna_voom_protein.RDS")
###############################################################################
#see normalization 
mRNA_voomt <- t(v_filtered)
htree_norm <- hclust(dist(mRNA_voomt), method = "average") #can take some time
dend3 <- as.dendrogram(htree_norm)
col_aa_red <- ifelse(grepl("GTEX", labels(dend3)), "red", "blue")
dend4 <- assign_values_to_leaves_edgePar(dend=dend3, value = col_aa_red, edgePar = "col") 
pdf(file="figures/htree_protein_gtex_tcga_norm_red.pdf", height=80, width=100) #išsaugijimui didesniu formatu
plot(dend4) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off()
#atsiskyre bet oh well