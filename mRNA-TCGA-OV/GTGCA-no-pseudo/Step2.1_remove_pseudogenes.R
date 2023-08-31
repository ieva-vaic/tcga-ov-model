#Remove pseudogenes from TCGA/GTEX database data
#first get all names and their types form biomart
setwd("~/rprojects/TCGA-OV-data") #wsl
library(biomaRt)

#get all the counts
gtcga <- readRDS("joined_gtex_tcga.RDS")
# regain count data: genes in rows
counts_gtcga <- gtcga[, 48:35165] #nusimazinuinu tik countus, bet dar liko gtex
which( colnames(counts_gtcga)=="gtex" ) #35118 stulpelis
counts_gtcga <- counts_gtcga[, -35118]
counts_gtcga <- t(counts_gtcga) #large numeric matrix with rows as genes
counts_gtcga <- as.data.frame(counts_gtcga)
#get all the ensembl names
counts_gtcga$ensembl <- rownames(counts_gtcga)
counts_gtcga$ensembl_gene_id <- gsub("\\..*", "",counts_gtcga$ensembl)
length(counts_gtcga$ensembl_gene_id) # 35117 genes

#Biomart: cia uztenka karta padaryt
listEnsembl() #shows available databases
ensembl <- useEnsembl(biomart = "genes") #unsuportive gali but
datasets <- listDatasets(ensembl) #issirinkti human genes
ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl') #name of data base ir data set
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)
#get gene names and biotype
gtcga_genes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', "gene_biotype"), #values to retreve
                     filters = "ensembl_gene_id", #input on quary
                     values = counts_gtcga$ensembl_gene_id,
                     mart = ensembl.con) #sukurtas conection object
dim(gtcga_genes) #34966 = 151 pamesta

#prisijoinsiu dfs
counts_gtcga_with_gene_names <- left_join(counts_gtcga, gtcga_genes, by= "ensembl_gene_id")
dim(counts_gtcga_with_gene_names) #35117   604

#remove pseudogenes
table(counts_gtcga_with_gene_names$gene_biotype, useNA = "a") #viso 34966 genu
#omg  2421 +  73+ 616 +81 +406+33 +13 +998 +10001 +187 +1+3+9+11+481 = 15334 random genu
#tec = to be experimentaly confirmed
unvanted <- c("unprocessed_pseudogene", "unitary_pseudogene", 
"transcribed_unprocessed_pseudogene", "transcribed_unitary_pseudogene",
"transcribed_processed_pseudogene",  "TR_V_pseudogene", "TR_J_pseudogene",
"pseudogene", "TEC", "processed_pseudogene", 'IG_V_pseudogene', 'IG_J_pseudogene',
"IG_pseudogene", "IG_C_pseudogene", "artifact", "rRNA_pseudogene", "unitary_pseudogene")

gtgca_names_sub <- counts_gtcga_with_gene_names[!(counts_gtcga_with_gene_names$gene_biotype %in% unvanted), ]
dim(gtgca_names_sub) #liko 19766 genu is 34966     
table(gtgca_names_sub$gene_biotype, useNA = "a") #lieka 124 NA
gtgca_names_sub <- gtgca_names_sub %>% drop_na(gene_biotype) #lieka 19642 genu

# negaliu pakeisti gene names � rownames: non-unique values when setting 'row.names
#todel grazinu transcriptu vardus
rownames(gtgca_names_sub) <- gtgca_names_sub$ensembl
repeating_genes <- c("", "5_8S_rRNA", "5S_rRNA", "7SK",
                     "ASMTL-AS1", "DHRSX-IT1", "DNAJC9-AS1",
                     "GOLGA8M", "GPR84-AS1", "GTPBP6", "Y_RNA",
                     "IL9R", "LINC00102", "LINC00106", "LINC00685",
                     "LINC01238", "LINC03023", "Metazoa_SRP",
                     "MIR3690", "MIR6089", "NPIPA9", "RAET1E-AS1", 
                     "SNORA62", "SNORA63", "SNORA70", "SNORA71", "SNORA72",
                     "SNORA73", "SNORA74", "SNORA75", "SNORD115", "SNORD116",
                     "SNORD18", "SNORD22", "SNORD27", "SNORD30",
                     "SNORD33", "SNORD39", "SNORD42", "SNORD63",
                     "SNORD81", "U1", "U2", "U3", "U4", "U6", "U7", 
                     "U8", "Vault", "WASH6P", "WASIR1")

repeating_counts <- gtgca_names_sub[(gtgca_names_sub$external_gene_name %in% repeating_genes), ]
table(repeating_counts$gene_biotype, useNA = "a") #pasalina 8717 transcriptu

final_counts_w_names <- gtgca_names_sub[!(gtgca_names_sub$external_gene_name %in% repeating_genes), ]
table(final_counts_w_names$gene_biotype, useNA = "a") #lieka 10925 transcriptu

rownames(final_counts_w_names) <- final_counts_w_names$external_gene_name
filtered_gtgca_counts <- final_counts_w_names[, -c(602:605)]
saveRDS(filtered_gtgca_counts, "filtered_gtgca_counts.RDS") 
#xx
################################################################################
