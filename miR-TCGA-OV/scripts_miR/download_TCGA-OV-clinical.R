#this is a STEP1.2: a script for downloading TCGA-OV clinical data
setwd("/home/ieva/rprojects/TCGA-OV-data/") #home wsl - data folder 
#load libraries
library("TCGAbiolinks")
#library(argparse)
# parse data
#parser <- ArgumentParser(description= 'this program downloads TCGA-OV clinical data')
#parser$add_argument('--output', '-o', help= 'Output file')
#xargs<- parser$parse_args()
#I need to download the clinial stuff seperatly. 
query_TCGA_clinical = GDCquery_clinic(
  project = "TCGA-OV" )
colnames(query_TCGA_clinical)
dim(query_TCGA_clinical) #too many samples
#add rownames for later
rownames(query_TCGA_clinical) <- query_TCGA_clinical$bcr_patient_barcode
#I will narrow down 608 samples later when combining data set with microRNA counts
#let's save!
saveRDS(object = query_TCGA_clinical,
        file = snakemake@output[[1]],
        compress = FALSE)