#l
setwd("~/tcga-ov-data")
train_data_df <- read.csv( "train_data_df_cleaned.csv", header = T) # train readus ir clinical atskirai
#sujoinint clinal su readais-countais
#paleist modeli, outputas bus predictions aba jo matas

setwd("~/tcga-ov-data")
train_counts_cleanded <- read.csv( "train_counts_cleaned.csv", header = T) #ant readu 
train_clinical_cleanded <- read.csv( "train_data_cleaned.csv", header = T) #ant clinical 

train_cleanded_joined <- semi_join(train_clinical_cleanded, train_counts_cleanded, by = c("barcode", "barcode"))


