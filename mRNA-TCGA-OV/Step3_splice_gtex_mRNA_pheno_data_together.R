#Step 3: splice gtex and TCGA-OV mRNA counts together
library(tidyverse)
library(RColorBrewer)
setwd("~/rprojects/TCGA-OV-data")
gtex_counts <- readRDS("GTEX/gtex_counts.RDS")

#the gtex is now a matrix that looks like this:
#rows: ENSEMBLIDS (transcipt level), all 56200 of them 
head(colnames(gtex_counts))
dim(gtex_counts)
#name: ENSEMBL IDS
#Description: real names 
#gtex_counts[3:181] #is case names 

#make a numeric dataframe from gtex
gtex_df <- as.data.frame(gtex_counts)
#read in tcga counts with names data for mRNA 
#tcga_counts <- readRDS("tcga_with_names.RDS") #nusalinti bevardziai pagal biomart
tcga_counts <- readRDS("tcga_with_names_all.RDS") ####

#get how many gene names match between gtex and tcga
tcga_names <- tcga_counts$external_gene_name
gtex_names <- gtex_df$Description
sum(gtex_names %in% tcga_names) #36056
sum(tcga_names %in% gtex_names) #35815
intersect_genes <- intersect(gtex_names, tcga_names)
##what are the genes that do not match?
not_matching <- tcga_counts[!(tcga_counts$external_gene_name %in% intersect_genes), ]
table(not_matching$gene_biotype, useNA="a") #1583 tiek protein coding prarandama!
##
rownames(gtex_df) <- gtex_df$Name
gtex_df[,3:182] <- lapply(gtex_df[,3:182], as.numeric)
dim(gtex_df) #56200 genes
gtex_df$external_gene_name <- gtex_df$Description
gtex_df$ensembl_gene_id <- gsub("\\..*", "",gtex_df$Name)
sum(gtex_df$ensembl_gene_id %in% tcga_counts$ensembl_gene_id) #55661 sutampa tiek, dingsta 539

##############################################################################
#the problem when joinging on gene names is that not all gene names is unique,
#this might be because some genes might have haplotypes
#ENSG00000117713.20 tcga
#ENSG00000117713.18 gtex
#Thus I need single line per gene on both count matrixes!
#choosing to splice on the ensemble gene names
#some genes (PAR_Y) has the same ensg name, thus need to delete them
dup_tcga <- duplicated(tcga_counts$ensembl_gene_id)
dup_tcga <- tcga_counts[dup_tcga, ] #visi parY
dup_tcga$ensembl
#all of these Y chromosome genes has no counts, I don´t want them anayways
dup_gtex <- duplicated(gtex_df$ensembl_gene_id) ##44
dup_gtex <- gtex_df[dup_gtex, ] 
dup_gtex$Name
#again, Y chromosome
ENTG_Y <- tcga_counts[grepl('_PAR_Y', tcga_counts$ensembl), ]
ENTG_Y_names <- ENTG_Y$ensembl
tcga_counts_filt <- tcga_counts[!(tcga_counts$ensembl %in% ENTG_Y_names), ]
rownames(tcga_counts_filt) <- tcga_counts_filt$ensembl_gene_id #
dim(tcga_counts_filt) #left:  60419   420
#now for gtex
ENTG_YG <- gtex_df[grepl('_PAR_Y', gtex_df$Name), ]
ENTG_YG_names <- ENTG_YG$Name 
gtex_counts_filt <- gtex_df[!(gtex_df$Name %in% ENTG_YG_names), ]
rownames(gtex_counts_filt) <- gtex_counts_filt$ensembl_gene_id #veik rownames 
dim(gtex_counts_filt) #left 56156   184
#join!
gtcga <- left_join(gtex_counts_filt, tcga_counts_filt, by = "ensembl_gene_id")
dim(gtcga) #56156   603
#some genes are only in GTEX!
sum(is.na(gtcga$external_gene_name.y)) #tiek genu ner tcga: 687, sakyciau ner daug
gtcga_na <- gtcga[is.na(gtcga$external_gene_name.y), ] 
#loose genes that are not in both gtex and tcga:
gtcga_final <- gtcga[!(is.na(gtcga$external_gene_name.y)), ] 
dim(gtcga_final)#final 55469 genes
saveRDS(gtcga_final, "gtcga_final.RDS")
##############################################################################
#remove non-protein-coding
biotypes <- table(gtcga_final$gene_biotype, useNA = "a")
biotypes <- as.data.frame(biotypes)
biotypes
png(file="figures/biotypes.png", height=1000, width=1000)
  ggplot(biotypes,
       aes(x = biotypes$Var1, y = biotypes$Freq, color = biotypes$Var1, fill=biotypes$Var1)) +
  geom_col() +
  labs(x = NULL, y = "Number of genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()

gtgca_protein <- gtcga_final[(gtcga_final$gene_biotype == "protein_coding" ), ]
dim(gtgca_protein) #liko 19197 genu is 55469     
table(gtgca_protein$gene_biotype, useNA = "a") 

################################################################################
#find repeating names
#I prefer biomart names, as some names in GTEX are synonyms in HUGO nomenclature
#however, some biomart names are missing, thus are going to be filled by gtex provided names
gtgca_protein["external_gene_name.y"][gtgca_protein["external_gene_name.y"] == ''] <- NA
any(is.na(gtgca_protein$external_gene_name.y))
gtgca_protein <- gtgca_protein %>% 
  mutate(external_gene_name.y = coalesce(external_gene_name.y, Description))  
any(is.na(gtgca_protein$external_gene_name.y))
#fix NPIPA9
which(gtgca_protein$external_gene_name.y == 'NPIPA9')
gtgca_protein[13803, 602] <- NA
gtgca_protein[13807, 602] <- NA
gtgca_protein <- gtgca_protein %>% 
  mutate(external_gene_name.y = coalesce(external_gene_name.y, Description))  
any(is.na(gtgca_protein$external_gene_name.y))
#add gene names as rownames
rownames(gtgca_protein) <- gtgca_protein$external_gene_name.y

#SAVE!
saveRDS(gtgca_protein, "gtcga_final_protein_w_biomart_names.RDS")

#loose the gene descriptions now
gtgca_final_no_names <- gtgca_protein[, -c(1:2, 183:184, 601:603)] 
dim(gtgca_final_no_names) #19197   596
saveRDS(gtgca_final_no_names, "gtcga_final_counts.RDS")



