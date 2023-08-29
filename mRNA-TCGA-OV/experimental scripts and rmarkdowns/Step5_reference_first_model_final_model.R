#l
library(dplyr)
library(ggplot2)
setwd("~/tcga-ov-data")
#sujoinint clinal su readais-countais
#paleist modeli, outputas bus predictions aba jo matas

train_counts_cleaned <- readRDS("train_counts_cleaned.RDS")
train_clinical_cleaned <- read.csv( "train_data_cleaned.csv", header = T) #ant clinical 

train_cleaned_joined <- full_join(train_clinical_cleaned, train_counts_cleaned, by = "barcode")

nrow(train_cleaned_joined)
ncol(train_cleaned_joined)

#siaip clean etape geriau sita perkelti
#rasti outcome dead alive stulpeli, padaryt table patikrinimui
table(train_cleaned_joined$vital_status, useNA='a')
#glm(outcome ~ x)
model <- glm(vital_status == "Dead" ~ 1, data= train_cleaned_joined, family = "binomial" ) #vidutine tikimybe, nuo efektines grupes patikirint art tikrai
summary(model) #intercepta rodo, rodo p ar daugiau uz 0

pred_prob <- predict(model, type = "response")
# dead/nrow vidutine tikimybe yra 219/347, todel base probability yra 0.6311239

pred_classes <- pred_prob > 0.5 #padaro true false, siuo atveju visi true
table(pred_classes, train_cleaned_joined$vital_status) #paglyginimas su real reiksmem 
a <- train_cleaned_joined$vital_status == "Dead" #cia reikalinga palyginimui tikru reiksmiu
v <- pred_classes == a  #cia lyginam predictions su reiksmem
sum(v) / nrow(train_cleaned_joined) #loss accurary siuo atveju paskaiciuoja same skaiciu 0.6311239

#2nd model 
model_clinical <- glm(vital_status == "Dead" ~ figo_stage, data= train_cleaned_joined,
             na.action = na.exclude, family = "binomial" ) #vidutine tikimybe, nuo efektines grupes bet patikirint ar tikrai
summary(model_clinical) #intercepta rodo

ggplot(train_cleaned_joined, aes(x=vital_status  , y=figo_stage)) +
  geom_boxplot() #optional see your stuff mapped

pred_prob_clin <- predict(model_clinical, type = "response")

pred_classes_clin <- pred_prob_clin > 0.5
table(pred_classes_clin, train_cleaned_joined$vital_status)
a <- train_cleaned_joined$vital_status == "Dead"
v <- pred_classes_clin == a  
sum(v, na.rm = T) / nrow(train_cleaned_joined) #loss accurary for figo not cleaned 0.6455331


#pagalvisim ar glm modelis is the best
#pagalvisim base clinikiniu duomenu 
#pagalvot outcome gal cotinuous gal yr iki progreso 
#looss iverti paglvot