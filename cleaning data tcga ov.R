#valymas, quality control
setwd("~/tcga-ov-data")
train_counts <- read.csv( "train_counts.csv", header = T) #ant readu 
write.csv(train_counts, "train_counts_cleaned.csv", quote = T, row.names = F)

train_clinical <- read.csv( "train_data.csv", header = T) #ant readu 
write.csv(train_clinical, "train_data_cleaned.csv", quote = T, row.names = F)


