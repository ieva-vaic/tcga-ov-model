#Step 4 split
coldata_tcga <- readRDS( "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/coldata_tcga_mrna.RDS")
tcga_mRNA_voom <-  readRDS("C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/tcga_mRNA_voom.RDS")
tcga_mRNA_voom <-t(tcga_mRNA_voom)
# Set (random-number-generator) seed so that results are consistent between runs
set.seed(18) #choose favorite number
train_ids <- rbinom(nrow(tcga_mRNA_voom), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
tcga_counts_train = tcga_mRNA_voom[train_ids, ] #new train data  326 60660
tcga_counts_test  = tcga_mRNA_voom[!train_ids, ] #new test data 103 60660

#add clinical data 

tcga_pheno_train = coldata_tcga[train_ids, ] #new train data 326  47
tcga_pheno_test  = coldata_tcga[!train_ids, ] #new test data 103  47

saveRDS(tcga_counts_train, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/train_mRNA_tcga_normcounts.RDS")
saveRDS(tcga_counts_test, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/test_mRNA_tcga_normcounts.RDS")

saveRDS(tcga_pheno_train, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/train_tcga_coldata.RDS")
saveRDS(tcga_pheno_test, "C:/Users/Ieva/Desktop/NVI GDL/R projetcs/tcga-ov-data/TCGA-OV-mRNA data/mRNA and clinical data TCGA-OV/test_tcga_coldata.RDS")