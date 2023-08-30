#Step 4: get joined tcga/gtex count data -> outlier chek with hplots -> lognormalize
library(tidyverse)
library(dendextend)
library(WGCNA)
library(GDCRNATools)
#get count matrix, rows as genes
#working with full data!
gtcga <- readRDS("C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/joined_gtex_tcga.RDS")
# regain count data: genes in rows
counts_gtcga <- gtcga[, 48:35165] #nusimažinu tik countus, bet dar liko gtex
which( colnames(counts_gtcga)=="gtex" ) #35118 stulpelis
counts_gtcga <- counts_gtcga[, -35118]
#transpose!
counts_gtcga <- t(counts_gtcga) #large numeric matrix with rows as genes
#regain pheno/coldata
which(colnames(gtcga)=="gtex" ) #35165 stulpelis
coldata_gtcga <- gtcga[, c(1:47, 35165)]
colnames(coldata_gtcga) #48
saveRDS(coldata_gtcga, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/coldata_gtcga_full.RDS")
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

pdf(file="C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/Figures mRNA/htree_full_gtex_tcga.pdf", height=30, width=60) #išsaugijimui didesniu formatu
plot(htree) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off()

#add colors to hclust, (library dendextend)
dend <- as.dendrogram(htree)
col_aa_red <- ifelse(grepl("GTEX", labels(dend)), "red", "blue")
dend2 <- assign_values_to_leaves_edgePar(dend=dend, value = col_aa_red, edgePar = "col") #

pdf(file="C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/Figures mRNA/htree_full_gtex_tcga_red.pdf", height=80, width=100) #išsaugijimui didesniu formatu
plot(dend2) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off()

#
#h <- heatmap(counts_gtcga, scale="column", col = cm.colors(256) )


###############################################################################
# Normalize, I will use GDC RNA tools
#1 type 2: from GDCRNAtools: gdcVoomNormalization: does TNB ir voom transformatio
mRNA_voom <- gdcVoomNormalization(counts = counts_gtcga, filter = FALSE) #
#pridedam normal vardus
mRNA_voomt <- t(mRNA_voom)

saveRDS(mRNA_voom, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/mrna_voom.RDS")
###############################################################################
#see normalization 
htree_norm <- hclust(dist(mRNA_voomt), method = "average") #can take some time
dend3 <- as.dendrogram(htree_norm)
col_aa_red <- ifelse(grepl("GTEX", labels(dend3)), "red", "blue")
dend4 <- assign_values_to_leaves_edgePar(dend=dend3, value = col_aa_red, edgePar = "col") 
pdf(file="C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/Figures mRNA/htree_full_gtex_tcga_norm_red.pdf", height=80, width=100) #išsaugijimui didesniu formatu
plot(dend4) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off()
#it is very bad, viskas idealiai atsiskyrę, su deseq irgi idealiai, o tcga nieko nepadaro lieka same kaip neclusterintas


