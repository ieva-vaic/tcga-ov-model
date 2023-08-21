#this is a STEP1.2: a script for downloading TCGA-OV clinical data
#setwd("/home/ieva/rprojects/TCGA-OV-data/") #home wsl - data folder 
#load libraries
library(tidyverse)

tcga_mir_LOGdata <- readRDS(snakemake@input[[1]])
full_clinical <- readRDS(snakemake@input[[2]])

#tcga_mir_LOGdata = readRDS(file = "tcga_mir_LOGdata.RDS") 
#full_clinical <- readRDS(file = "full_clinical.RDS")
#i will first transpose the count data 
tcga_mir_LOGdata_t <- as.data.frame(t(tcga_mir_LOGdata))
dim(tcga_mir_LOGdata_t) #499 cases ir 1881 mirs
#as there are no full barcodes in clinical data need to fix
tcga_mir_LOGdata_t$bcr_patient_barcode <- rownames(tcga_mir_LOGdata_t)
tcga_mir_LOGdata_t$bcr_patient_barcode <- gsub("-01.*","", tcga_mir_LOGdata_t$bcr_patient_barcode)
#these 2A samples should be removed as they are "Recurrent Solid Tumor"
tcga_mir_LOGdata_t <- tcga_mir_LOGdata_t[!grepl("-02A", tcga_mir_LOGdata_t$bcr_patient_barcode),] #491 cases
str(tcga_mir_LOGdata_t$bcr_patient_barcode) #491 names liko
#still 2 non unique zmones yra:
x <- data.frame(table(tcga_mir_LOGdata_t$bcr_patient_barcode))
x <- x[x$Freq == 2, ]
xnames <- c("TCGA-09-0366",
            "TCGA-23-1023")
x <- filter(tcga_mir_LOGdata_t, tcga_mir_LOGdata_t$bcr_patient_barcode %in% xnames)
rownames(x)
#pašalinu nes su 1R žymėm pavadinimas "TCGA-23-1023-01R-01R-1564-13"
#"TCGA-09-0366-01A-01R-1986-13" is in the "removed samples list in firehose:
#https://gdac.broadinstitute.org/runs/tmp/sample_report__2018_01_17/Replicate_Samples.html
x <- c("TCGA-09-0366-01A-01R-1986-13", "TCGA-23-1023-01R-01R-1564-13" )
tcga_mir_LOGdata_t <- tcga_mir_LOGdata_t[!(row.names(tcga_mir_LOGdata_t) %in% x),]
dim(tcga_mir_LOGdata_t) #489 liko
#filter clinical data set (so that the extra clinical data cases are not present)
sum(full_clinical$bcr_patient_barcode %in% tcga_mir_LOGdata_t$bcr_patient_barcode) #bus 489 zmones
mir_data_LOGdf <- semi_join(full_clinical, tcga_mir_LOGdata_t, by='bcr_patient_barcode')
dim(mir_data_LOGdf) #489 cases
#first I wanna add the full barcodes to the phenodata
#best way is to join clinical data (mir_data_df) with count data (tcga_mir_data_t)
#first add back the full barcodes to 
tcga_mir_LOGdata_t$full_barcode <- rownames(tcga_mir_LOGdata_t)
full_LOGdata_mirs <- full_join(tcga_mir_LOGdata_t, mir_data_LOGdf, "bcr_patient_barcode" )
dim(full_LOGdata_mirs) # 1340 + 66 = 1405 nes liko tik vienas stulpelis per kuri junge!
colnames(full_LOGdata_mirs[1337:1405]) #see the join at the end
#saveRDS(full_LOGdata_mirs, xargs$output)
saveRDS(full_LOGdata_mirs, snakemake@output[[1]])