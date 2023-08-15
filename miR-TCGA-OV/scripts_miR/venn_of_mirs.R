#this script produces the interconnections between datasets
#setwd("/home/ieva/rprojects/TCGA-OV-data") #home wsl - data folder 
#load libraries
library(tidyverse)
library(ComplexHeatmap)
library(gplots)

relevant_mirs_log_status <- read.csv(snakemake@input[[1]], header = T)
relevant_mirs_log_stage <- read.csv(snakemake@input[[2]], header = T)
relevant_mirs_3figo_log <- read.csv(snakemake@input[[3]], header = T)

#no snakemake
#relevant_mirs_log_status <- read.csv("mirs_status.csv", header = T)
#relevant_mirs_log_stage <- read.csv("mirs_stage2.csv", header = T)
#relevant_mirs_3figo_log <- read.csv("mirs_stage3.csv", header = T)

#turn back into vectors
relevant_mirs_log_status <- relevant_mirs_log_status$x
relevant_mirs_log_stage <- relevant_mirs_log_stage$x
relevant_mirs_3figo_log <- relevant_mirs_3figo_log$x

#from this script
relevant_mirs_log_status #43
relevant_mirs_log_stage #24
relevant_mirs_3figo_log #75

#upset is a venn plot
list_of_relevant_mirs <- list(
  log_stage2 = relevant_mirs_log_stage,
  log_stage3 = relevant_mirs_3figo_log, 
  log_status = relevant_mirs_log_status
)
matrix_of_relevant_mirs <- list_to_matrix(list_of_relevant_mirs)
m = make_comb_mat(matrix_of_relevant_mirs)
upsetM <- UpSet(m)
png(filename=snakemake@output[[1]])
plot(upsetM)
dev.off()

#find same mirs between 3 lists
venn_diagram <- venn(list_of_relevant_mirs)
png(filename=snakemake@output[[2]])
plot(venn_diagram)
dev.off()

mirs_status_ <- intersect(relevant_mirs_log_stage, relevant_mirs_3figo_log )
mirs_status_
write.csv(mirs_status_, snakemake@output[[3]], quote = F, row.names = T)


