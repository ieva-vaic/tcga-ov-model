Kiekvienas relevant failas iš kiekvieno scripto:
DUOMENŲ PARSISIUNTIMAS:
STEP0.gtex: "gtex_counts.RDS"
STEP0.1.xena: "00_ClinTraits.csv"
STEP1.tcga: "tcga_data.RDS"; "pheno.RDS"
DUOMENŲ TVARKYMAS
STEP2.1.PHENO: "pheno_no_empty_data.RDS" -> TCGA PHENO BE NEREIKALINGŲ DUOMENŲ
 "joinedTCGA_XENA_clinical.RDS" -> APJUNGTAS TCGA IR XENA PHENO
STEP2.2. REMOVED WEIRD CASES FROM TCGA:
7 Recurrent Solid Tumor
1 Prior treatment
5 Kitos ligos nei serous cystadenocarcinoma
"tcga_no_weird_full.RDS"-> 416 žmonės (tcga) 60727 (60660 genai + 67 clinical) 
"tcga_no_weird_counts.RDS" -> 416 tcga žmonės  67 pheno duomenys
"tcga_no_weird_pheno_XENA_TCGA.RDS" -> 416 60660
STEP2.3.biomart. prideda genu vardus ir funckijas 
"tcga_with_names_all.RDS"
STEP3.splice.all.data and removes non-protein coding genes. 
"gtcga_final.RDS" -> 55469 genai  603 tcga ir gtex + external genu vardai (basically raw apjuntas)
"gtcga_final_counts.RDS" -> 19197 genai 596 tcga ir gtex žmonės, eilutėse genų pavadinimai readable
STEP4.normalization. otliers, htrees, gdc tools normalization, filtered, with model including the study. 
"mrna_voom_protein.RDS"
STEP5.train-test-split:
"train_gtcga_normcounts_prot.RDS" #153gtex  336tcga
"test_gtcga_normcounts_prot.RDS" #27gtex   79 tcga
STEP6.Elastic_net: 
"gtcga_elastic.RDS" -> pasirinktų ~200 genų vardų.
"elastic_net_model_gtex.RDS" -> elastic_net glm object
SKIP SIDE STUFF. TIK STEP7 nufiltruoja pheno kad būtų tik test data.
STEP10.COX: coxnet ir tada survroc ant train data
"pheno_survival_only.RDS" -> pheno tik su survival, tcga train žmonėms
"res_coef_coxnet_names.csv" -> coxnet išsirinktų 10 genų vardai
"coxnet_fit.RDS" -> glmnet(surv_counts, y2, family="cox")
"coxnet_cvfit.RDS" -> cv.glmnet(surv_counts, y2, family = "cox", type.measure = "C") -> cia lamdas nubraižys
STEP12.testai: multiple survrocs saved
