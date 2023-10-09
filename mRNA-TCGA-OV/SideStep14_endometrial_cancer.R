#download other cancer data
# Load packages
library("TCGAbiolinks")
library("SummarizedExperiment")
setwd("/home/ieva/rprojects/TCGA-OV-data/") #wsl
library(tidyverse)
library(biomaRt)
library(tidyverse)
library(dendextend)
library(WGCNA)
library(GDCRNATools)
library(edgeR)
library(limma)
library(DESeq2)
## Build your query
query_TCGA = GDCquery(
  project = "TCGA-UCEC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts", 
  access = "open" )

##Download your data (do not forget to setwd to data folder, creates a big folder)
GDCdownload(query = query_TCGA, files.per.chunk = 200) 
##Prepare your data
tcga_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
ucec <- assay(tcga_data)
dim(ucec)
saveRDS(ucec, "ucec.RDS")
pheno_ucec <-  as.data.frame(colData(tcga_data)) 
dim(pheno_ucec)
#################################################################################################
ucec <- readRDS("ucec.RDS")
ucec_df <- as.data.frame(ucec) 
dim(ucec_df)#60660   589 
ucec_df[,1:589] <- lapply(ucec_df[,1:589], as.numeric)
str(ucec_df) #60660 genes
#get gene names in the count df
ucec_df$ensembl <- rownames(ucec_df)
ucec_df$ensembl_gene_id <- gsub("\\..*", "",ucec_df$ensembl)
length(ucec_df$ensembl_gene_id) # 60660 
#biomart
listEnsembl() #shows available databases
ensembl <- useEnsembl(biomart = "genes") #unsuportive gali but
datasets <- listDatasets(ensembl) #issirinkti human genes
ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl') #name of data base ir data set
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)
ucec_genes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', "gene_biotype"), #values to retreve
                     filters = "ensembl_gene_id", #input on quary
                     values = ucec_df$ensembl_gene_id,
                     mart = ensembl.con) #sukurtas conection object
dim(ucec_genes) #60419 = 241 pamesta
counts_ucec_with_gene_names <- left_join(ucec_df, ucec_genes, by= "ensembl_gene_id")
dim(counts_ucec_with_gene_names) #60660   593
tail(colnames(counts_ucec_with_gene_names))
saveRDS(counts_ucec_with_gene_names, "counts_ucec_with_gene_names.RDS")
#####################################################################################################
counts_ucec_with_gene_names <- readRDS("counts_ucec_with_gene_names.RDS")
gtcga <- readRDS("gtcga_final_protein_w_biomart_names.RDS")
tail(colnames(gtcga))
tail(colnames(counts_ucec_with_gene_names))
names(counts_ucec_with_gene_names)[1:589] <- paste('UCEC',
 names(counts_ucec_with_gene_names)[1:589], sep = '_')
colnames(counts_ucec_with_gene_names) 
gtcga_ucec <- left_join(counts_ucec_with_gene_names, gtcga, by = "ensembl_gene_id")
dim(gtcga_ucec) #zodziu visi genai liko, dabar as noresiu kad liktu tik tie, kurie yra gtcga
final_genes <- rownames(gtcga)
final_genes
gtcga_ucec_filtered <- gtcga_ucec[gtcga_ucec$external_gene_name.y %in% final_genes, ]
colnames(gtcga_ucec_filtered)
rownames(gtcga_ucec_filtered) 

ENTG_YG <- gtcga_ucec_filtered[grepl('_PAR_Y', gtcga_ucec_filtered$ensembl.x), ]
ENTG_YG
ENTG_YG_names <- ENTG_YG$ensembl.x #turbut tie patys bet anyways
gtcga_ucec_filtered <- gtcga_ucec_filtered[!(gtcga_ucec_filtered$ensembl.x %in% ENTG_YG_names), ]
rownames(gtcga_ucec_filtered) <- gtcga_ucec_filtered$external_gene_name.y
colnames(gtcga_ucec_filtered)
gtcga_ucec_filtered <- gtcga_ucec_filtered[,-c(590:595, 776, 1193:1195)]
gtcga_ucec <- as.data.frame(t(gtcga_ucec_filtered))
gtcga_ucec$study <- substr(rownames(gtcga_ucec), 1, 4)
gtcga_ucec$study
saveRDS(gtcga_ucec, "gtcga_ucec_counts.RDS")
#normalize
gtcga_ucec_counts <- gtcga_ucec[, -19198]
gtcga_ucec_counts <- data.matrix(t(gtcga_ucec_counts))
group <- gtcga_ucec$study

expr = DGEList(counts = gtcga_ucec_counts)
expr = calcNormFactors(expr)
keepALL <- rowSums(cpm(expr) > 1) >= 0.5*ncol(gtcga_ucec_counts)
nGenes <- as.numeric(summary(keepALL)[2]) + 
  as.numeric(summary(keepALL)[3])
nKeep <- summary(keepALL)[3]

cat (paste('Total Number of genes: ', nGenes, '\n', sep=''))
cat (paste('Number of genes for downstream analysis: ', nKeep, 
           '\n', sep=''))

exprALL <- expr[keepALL,,keep.lib.sizes = TRUE]
v_filtered <- voom(exprALL, design= model.matrix(~0 + group), plot = TRUE)$E
saveRDS(v_filtered, "ucec_gtex_count_voom.RDS")
v_filtered_df <- data.frame(t(v_filtered))
################################################################################
#VEL BEDA KAD NER GRUPIU
v_filtered_df$group <- substr(rownames(gtcga_ucec), 1, 4)
################################################################################
res_coef_cox_names <- read.csv("res_coef_coxnet_names.csv")
res_coef_cox_names <- res_coef_cox_names$x
res_coef_cox_names <- c(res_coef_cox_names, "group") #cia reik priappendint kad nufiltruotu ka reik
#now for the boxplots
png("figures/ucec_boxplot.png", width=600, height=500, res=120)
v_filtered_df %>% dplyr::select(res_coef_cox_names) %>%
  pivot_longer(., cols = c(res_coef_cox_names[1:10]), names_to = "Genes", values_to = "EXPR") %>%
  ggplot(aes(x = Genes, y = EXPR, fill = group)) +
  geom_boxplot()+
  ylab("Normalized expression")+
  theme(
    axis.text.x = element_text(face = "italic", angle=-90))+
  guides(fill=guide_legend(title="Study"))
dev.off()