#Step 2: splice gtex and TCGA-OV mRNA counts together
library(tidyverse)
library("TCGAbiolinks")
library("SummarizedExperiment")
setwd("~/rprojects/TCGA-OV-data")
gtex_counts <- readRDS("GTEX/gtex_counts.RDS")

#the gtex is now a matrix that looks like this:
#rows: ENSEMBLIDS (transcipt level), all 56200 of them 
head(colnames(gtex_counts))
#name: ENSEMBL IDS
#Description: real names 
#gtex_counts[3:181] #is case names 

#make a numeric dataframe from gtex
gtex_df <- as.data.frame(gtex_counts)
#read in tcga counts with names data for mRNA 
tcga_counts <- readRDS("tcga_with_names.RDS") ####

#get how many gene names match between gtex and tcga
tcga_names <- tcga_counts$external_gene_name
gtex_names <- gtex_df$Description
sum(gtex_names %in% tcga_names) #36056
sum(tcga_names %in% gtex_names) #35815
sutampa <- intersect(gtex_names, tcga_names)
"ARID1A" %in% sutampa #yra zinomi genai
"NOTCH1" %in% sutampa #yra zinomi genai

##what are the genes that do not match?
not_matching <- tcga_counts[!(tcga_counts$external_gene_name %in% sutampa), ]
table(not_matching$gene_biotype, useNA="a") #1583 tiek protein coding prarandama!
##
rownames(gtex_df) <- gtex_df$Name
gtex_df[,3:180] <- lapply(gtex_df[,3:180], as.numeric)
dim(gtex_df) #56200 genes
gtex_df$external_gene_name <- gtex_df$Description
gtex_df$ensembl_gene_id <- gsub("\\..*", "",gtex_df$Name)
sum(gtex_df$ensembl_gene_id %in% tcga_counts$ensembl_gene_id) #55513 sutampa tiek

#also multiple matches in left_join on ensembl_gene_id, so not an option
##############################################################################
#the problem when joinging on gene names is that not all gene names is unique,
#this is because some gene might have haplotypes
#jei jungt per transcriptus irgi not good
gtex_ARID <- gtex_df %>% 
  filter(external_gene_name == "ARID1A")
tcg_arid <- tcga_counts %>% 
  filter(external_gene_name == "ARID1A")
#skiriasi transciprtai
#ENSG00000117713.20 tcga
#ENSG00000117713.18 gtex
################################################################################
#I need single line per gene on both count matrixes!!!!!!!!
#one option is for the conflicted genes leave the higest expression transcript
#rownames(tcga_counts) <- tcga_counts$external_gene_name
dup_tcga <- duplicated(tcga_counts$ensembl_gene_id)
dup_tcga <- tcga_counts[dup_tcga, ] #visi parY
dup_tcga_counts <- dup_tcga[,1:429]
rowSums(dup_tcga_counts) 
#all of these genes has no counts, I don´t want them anayways
dup_gtex <- duplicated(gtex_df$ensembl_gene_id) ##44
dup_gtex <- gtex_df[dup_gtex, ] #visi turi "PAR_Y", t.y yra Y chormosomos, neturetu but raiskos moteryje
dup_gtex_counts <- dup_gtex[,3:180]
rowSums(dup_gtex_counts)  #visi 0, kaip ir turi but
## delete dupplicated - gal uzteks visus parY nutrinti?
ENTG_Y <- tcga_counts[grepl('_PAR_Y', tcga_counts$ensembl), ]
ENTG_Y_names <- ENTG_Y$ensembl
tcga_counts_filt <- tcga_counts[!(tcga_counts$ensembl %in% ENTG_Y_names), ]
rownames(tcga_counts_filt) <- tcga_counts_filt$ensembl_gene_id #
dim(tcga_counts_filt) #lieka  60419   433
#now for gtex
ENTG_YG <- gtex_df[grepl('_PAR_Y', gtex_df$Name), ]
ENTG_YG_names <- ENTG_YG$Name #turbut tie patys bet anyways
gtex_counts_filt <- gtex_df[!(gtex_df$Name %in% ENTG_YG_names), ]
rownames(gtex_counts_filt) <- gtex_counts_filt$ensembl_gene_id #veik rownames 
#join!
gtcga <- left_join(gtex_counts_filt, tcga_counts_filt, by = "ensembl_gene_id")
dim(gtcga) #56156   616
#some genes are only in GTEX!
sum(is.na(gtcga$external_gene_name.y)) #tiek genu ner tcga: 687, sakyciau ner daug
gtcga_na <- gtcga[is.na(gtcga$external_gene_name.y), ] 
#paziurejus visokie random genai, vistiek per juos skaiciuot nieko negalesiu todel pasalinu
gtcga_final <- gtcga[!(is.na(gtcga$external_gene_name.y)), ] #final 55469 genes
saveRDS(gtcga_final, "gtcga_final.RDS")
###############################################################################
#remove weird cases, the list is from other script
weird_cases <- c("TCGA-13-0913-02A-01R-1564-13", "TCGA-13-1489-02A-01R-1565-13", "TCGA-29-1705-02A-01R-1567-13",
                 "TCGA-29-1707-02A-01R-1567-13", "TCGA-29-1770-02A-01R-1567-13", "TCGA-29-2414-02A-01R-1569-13",
                 "TCGA-61-2008-02A-01R-1568-13", "TCGA-61-1721-01A-01R-1569-13")
gtcga_final_8less <- gtcga_final[, !(colnames(gtcga_final) %in% weird_cases)] 
saveRDS(gtcga_final_8less, "gtcga_final_without_weird.RDS")
##############################################################################
#remove pseudogenes
table(gtcga_final_8less$gene_biotype, useNA = "a")
#tec = to be experimentaly confirmed
unvanted <- c("unprocessed_pseudogene", "unitary_pseudogene", 
              "transcribed_unprocessed_pseudogene", "transcribed_unitary_pseudogene",
              "transcribed_processed_pseudogene",  "TR_V_pseudogene", "TR_J_pseudogene",
              "pseudogene", "TEC", "processed_pseudogene", 'IG_V_pseudogene', 'IG_J_pseudogene',
              "IG_pseudogene", "IG_C_pseudogene", "artifact", "rRNA_pseudogene",
              "unitary_pseudogene", "translated_processed_pseudogene")
gtgca_pseudo <- gtcga_final_8less[!(gtcga_final_8less$gene_biotype %in% unvanted), ]
dim(gtgca_pseudo) #liko 39630 genu is 55469     
table(gtgca_pseudo$gene_biotype, useNA = "a") 

################################################################################
#find repeating names
n_occur <- data.frame(table(gtgca_pseudo$Description))
xx <- n_occur[n_occur$Freq > 1,]
reaccuring_gtcga <- gtgca_pseudo[gtgca_pseudo$Description %in% n_occur$Var1[n_occur$Freq > 1],] 
#1700 transcriptu priklauso genams su 37168 priklausantiems 147
table(reaccuring_gtcga$gene_biotype) #10 tik protein coiding
protein_miss <- reaccuring_gtcga %>% filter(gene_biotype == "protein_coding")
#[1] "RGS5"     "C2orf61"  "C2orf61"  "CYB561D2" "MAL2"     "LYNX1"    "LYNX1"    "SPATA13"  "GOLGA8M" "ELFN2" 

#
gtgca_final <- gtgca_pseudo[!(gtgca_pseudo$Description %in% n_occur$Var1[n_occur$Freq > 1]),] 
rownames(gtgca_final) <- gtgca_final$Description
saveRDS(gtgca_final, "gtcga_final_without_weird_genes_and_samples.RDS")
#loose the gene descriptions now
gtgca_final_no_names <- gtgca_final[, -c(1:2, 183:184, 606:608)] #37930   601
saveRDS(gtgca_final_no_names, "gtcga_final_counts.RDS")



################################################################################
#only protein-coding names 
gtcga_protein <- filter(gtcga_final_8less,  gtcga_final_8less$gene_biotype == "protein_coding")
dim(gtcga_protein) #lieka 19197
#idomu, kad xena tutoriale taip pat yra protein coding ju pateikiamoje lenteleje panasiai genu, t.y. 19116
#find repeating names
n_occur <- data.frame(table(gtcga_protein$Description))
xx <- n_occur[n_occur$Freq > 1,]
reaccuring_gtcga <- gtcga_protein[gtcga_protein$Description %in% n_occur$Var1[n_occur$Freq > 1],] 
reaccuring_gtcga$Description #"C2orf61" "LYNX1" kartojas 2 kartus, mesiu lauk
repeating_genes <- c("C2orf61" ,"LYNX1")
gtcga_protein <- gtcga_protein[!(gtcga_protein$Description %in% n_occur$Var1[n_occur$Freq > 1]),] 
rownames(gtcga_protein) <- gtcga_protein$Description

gtcga_protein <- gtcga_protein[, -c(1:2, 183:184, 606:608)] #19193   601
saveRDS(gtcga_protein, "gtcga_proteins.RDS")

genesPC = fread("~/rprojects/XENA/zz_gene.protein.coding.csv")
xena_proteins <- genesPC$Gene_Symbol
biomart_proteins <- gtcga_protein$Description
length(intersect(xena_proteins, biomart_proteins)) #18290 tiek pavadinimu sutampa tarp duombaziu

