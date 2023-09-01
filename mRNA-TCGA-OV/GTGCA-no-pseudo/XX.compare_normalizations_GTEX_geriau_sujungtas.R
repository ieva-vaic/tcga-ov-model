#cia prie mano jau geriau apjungto Gtex seto darau kitus normalizavimus ir paleiziu lasso
library(GDCRNATools)
#get count matrix, rows as genes
#working with full data!
setwd("~/rprojects/TCGA-OV-data") #wsl
counts_gtcga <- readRDS("gtcga_final_counts.RDS") #large numeric matrix with rows as genes
##############################################################################
counts_gtcga <- data.matrix(counts_gtcga)

# Normalize, GDC RNA tools, with filter
#1 type 2: from GDCRNAtools: gdcVoomNormalization: does TNB ir voom transformation
mRNA_voom2 <- gdcVoomNormalization(counts = counts_gtcga, filter = TRUE) #
#pridedam normal vardus
mRNA_voom2t <- t(mRNA_voom2)

saveRDS(mRNA_voom2, "mrna_voom_no_pseudo_filter.RDS")

#split
gtcga_counts2 <- t(mRNA_voom2)
gtcga_counts2 <- as.data.frame(gtcga_counts2)
# Set (random-number-generator) seed so that results are consistent between runs
set.seed(18) #choose favorite number
train_ids <- rbinom(nrow(gtcga_counts2), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
gtex_counts_train2 = gtcga_counts2[train_ids, ] #new train data 107 35117 ->   494 14776
gtex_counts_test2  = gtcga_counts2[!train_ids, ] #new test data 107 35117 -> 107 14776

#
snames = rownames(gtex_counts_train2);
group = substr(snames, 1, 4); #Sets up level information for samples.
group = as.factor(group)
voom2_counts_train <- data.matrix(gtex_counts_train2)
res_num = cv.glmnet(
  x = voom2_counts_train,
  y = group,
  alpha = 1, 
  family = "binomial"
)
res_num

res_coef_recoded_num = coef(res_num, s="lambda.min") # the "coef" function returns a sparse matrix
# get coefficients with non-zero values
res_coef_recoded_num = res_coef_recoded_num[res_coef_recoded_num[,1] != 0,] 
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef_recoded_num = res_coef_recoded_num[-1]
relevant_genes_recoded_num = names(res_coef_recoded_num) # get names of the (non-zero) variables.
relevant_genes_recoded_num # 12

#Ar kazkas sutampa su mano "gerai" sujungtu gtex?

nefiltruotas <- c("EVA1B" ,  "COQ8A" ,  "NICN1" ,  "ACAD11",  "CCDC39"  ,"ZBTB9"  , 
                  "TOMM6" ,  "ILK"  ,   "EEF1G"   ,"TAX1BP3",
                  "GABARAP", "WDR83OS" ,"MT-TP")

intersect(relevant_genes_recoded_num, nefiltruotas)

################################################################################
#add model, pvz is XENA tutorialo
library(edgeR)
library(limma)
library(DESeq2)
counts_gtcga <- data.matrix(counts_gtcga)

x = DGEList(counts_gtcga)
snames = colnames(x);
group = substr(snames, 1, 4); #Sets up level information for samples.
x$samples$group = group #Assigns samples to appropriate group.
dim(counts_gtcga) #10925   601

cpm = cpm(x);
lcpm = cpm(x, log = TRUE);
L = mean(x$samples$lib.size) * 1e-6;
L #Displays average library size in millions.

M = median(x$samples$lib.size) * 1e-6;
M #Displays median library size in millions.

table(x$samples$group) #Returns number of samples per group.

keep.exprs = filterByExpr(x, group = group);
x = x[keep.exprs, , keep.lib.sizes = FALSE];
dim(x) #Returns number of genes and samples retained.

x = calcNormFactors(x, method = "upperquartile")
design = model.matrix(~0 + group)
colnames(design) = gsub("group", "", colnames(design));
contr.matrix = makeContrasts(TCGAvsGTEX = TCGA - GTEX,
                             levels = colnames(design))
v = voom(x, design, plot = TRUE)
voomExpr = v$E

#split
gtcga_counts3 <- t(voomExpr)
gtcga_counts3 <- as.data.frame(gtcga_counts3)
# Set (random-number-generator) seed so that results are consistent between runs
set.seed(18) #choose favorite number
train_ids <- rbinom(nrow(gtcga_counts2), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
gtex_counts_train3 = gtcga_counts3[train_ids, ] #new train data ->   494 21894
gtex_counts_test3  = gtcga_counts3[!train_ids, ] #new test data ->  107 21894

####
gtex_counts_train3 <- t(gtex_counts_train3)
#model
snames1 = rownames(gtex_counts_train3);
group1 = substr(snames1, 1, 4); #Sets up level information for samples.
group1 = as.factor(group1) #494
res_num1 = cv.glmnet(
  x = gtex_counts_train3,
  y = group1,
  alpha = 1, 
  family = "binomial"
)
res_num1

res_coef_recoded_num1 = coef(res_num1, s="lambda.min") # the "coef" function returns a sparse matrix
# get coefficients with non-zero values
res_coef_recoded_num1 = res_coef_recoded_num1[res_coef_recoded_num1[,1] != 0,] 
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef_recoded_num1 = res_coef_recoded_num1[-1]
relevant_genes_recoded_num1 = names(res_coef_recoded_num1) # get names of the (non-zero) variables.
relevant_genes_recoded_num1 # 10

intersect(relevant_genes_recoded_num, relevant_genes_recoded_num1)
intersect(relevant_genes_recoded_num1, nefiltruotas)
