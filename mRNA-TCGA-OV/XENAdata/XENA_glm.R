##my offf the rails analysis
library(glmnet)
library(dendextend)
#set dir to XENA
setwd("~/rprojects/XENA")

v <- readRDS("v.RDS") #iki glm
v_counts <- v$E
v_counts <- as.data.frame(v_counts)

#see normalization 
htree_norm <- hclust(dist(t(v_counts)), method = "average") #can take some time
dend3 <- as.dendrogram(htree_norm)
col_aa_red <- ifelse(grepl("GTEX", labels(dend3)), "red", "blue")
dend4 <- assign_values_to_leaves_edgePar(dend=dend3, value = col_aa_red, edgePar = "col") 
pdf(file="hclust_mano.pdf", height=80, width=100) #išsaugijimui didesniu formatu
plot(dend4) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off() #it is the same shiiiit

#get test data
xena_counts <- t(v_counts)
xena_counts <-as.data.frame(xena_counts)
set.seed(18) #choose favorite number
train_ids <- rbinom(nrow(xena_counts), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
xena_counts_train = xena_counts[train_ids, ] #new train data 407 3943
xena_counts_test  = xena_counts[!train_ids, ] #new test data  99 3943

snames = rownames(xena_counts_train);
group = substr(snames, 1, 4); #Sets up level information for samples.
group = as.factor(group)
xena_counts_train <- data.matrix(xena_counts_train)
res_num = cv.glmnet(
  x = xena_counts_train,
  y = group,
  alpha = 1, 
  family = "binomial"
)
res_num #

res_coef_recoded_num = coef(res_num, s="lambda.min") # the "coef" function returns a sparse matrix
# get coefficients with non-zero values
res_coef_recoded_num = res_coef_recoded_num[res_coef_recoded_num[,1] != 0,] 
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef_recoded_num = res_coef_recoded_num[-1]
relevant_genes_recoded_num = names(res_coef_recoded_num) # get names of the (non-zero) variables.
relevant_genes_recoded_num # 18


