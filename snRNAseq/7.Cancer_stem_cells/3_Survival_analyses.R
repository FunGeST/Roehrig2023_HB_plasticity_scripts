
# choose if Aurore or Theo  -----------

user <- c("Aurore", "Theo")[1] # choisir ici

# libraries and functions -----------
library(stringr)
library(survival)
library(rms) # for Kapplan Meier
library(magrittr)
library(dplyr)
library(readxl)
library(prismatic) # for color




LoadRdata <- function(file.name, 
                      print.dim = T, 
                      print.real.name = T) {
  file.real.name = load(file.name, verbose = print.real.name)
  if (print.dim) cat("dimensions: ", dim(get(file.real.name)))
  return(get(file.real.name))
}

factoall <- function (d, ncmax = 10) {
  n <- ncol(d)
  for (i in 1:n) {
    if (is.factor(d[, i])) {
      d[, i] <- as.character(d[, i])
      na <- which(is.na(d[, i]))
      num <- suppressWarnings(as.numeric(d[, i]))
      nanum <- which(is.na(num))
      if (length(nanum) == length(na)) {
        d[, i] <- num
      }
    }
  }
  d
}

ChangeDateFormat <- function(date.string.vec, as.date = F, american = F) {
  month.rank = if (american) 1:2 else 4:5
  day.rank = setdiff(c(1, 2, 4, 5), month.rank)
  date.string.vec = paste0(substring(date.string.vec, 7), "-",
                           substring(date.string.vec, month.rank[1], month.rank[2]), "-",
                           substring(date.string.vec, day.rank[1], day.rank[2]))
  if (as.date) date.string.vec = as.Date(date.string.vec, "%Y-%m-%d")
  
  date.string.vec
}

FormatPValue = function(pval) {
  if (pval >= 0.01) return(round(pval, 3))
  return(formatC(pval, format = "e", digits = 2))
}

GivePatientsSamples <- function(annot) {
  print(glue::glue("{length(annot$CHCID)} samples from {length(unique(annot$Patient.identification))} patients"))
}

CreateDirIfNecessary <- function(dir_path) {
  if (!file.exists(dir_path)) dir.create(dir_path)
}
# 1.2 Directories ------------
GEPELIN_dropbox_dir <- c("C:/Users/Aurore/Dropbox/labo/GEPELIN_GENOMIC_ANALYSES",
                         "D:/TZH/Dropbox/GEPELIN_GENOMIC_ANALYSES")[1 + (user == "Theo")] # ne pas changer

resdir <- file.path(GEPELIN_dropbox_dir, "ANALYSES/Correlations_integrated_table")
data_dir <- file.path(GEPELIN_dropbox_dir, "DATA")
computer <- c("Mac", "PC")[2]
hepato_dir <- ifelse(computer == "Mac",
                     "/Volumes/HEPATO PARTAGE",
                     "//10.93.23.19/hepato partage")
fluidigm_dir = file.path(hepato_dir, 
                         "TRANSCRIPTOME-RNA&miRNA/FLUIDIGM-NANOSTRING/FLUIDIGM/RESULTS_FLUIDIGM/2-DeltaDeltaCT")

# 1.3 Parameters ------------
# 1.3.1 Samples to test / filter ----------------
keep_only_one_tumor_by_patient <- c(T, F)[1]
tumor_to_select <- c("worst", "representative")[2] # choose how to select one tumor per patient
include_sample_without_chemo <- c(T, F)[1]
keep_only_fetal_majority <- c(T, F)[2]
keep_only_E_inf_30 <- c(T, F)[2]
if (keep_only_E_inf_30) {
  print("Warning: the condition about E < 30 overrides the condition about fetal majority")
  keep_only_fetal_majority <- F
}
series_to_work_on <- c("NGS", "Validation", "NGS_and_validation")[3]

#annot %>% 
# filter(Histological.Diagnosis != "HB") %>% 
#select(To_keep_for_survival_analysis_worst, CHCID)

# for (series_to_work_on in c("NGS", "Validation", "NGS_and_validation")) {

# for (keep_only_one_tumor_by_patient in c(T, F)) {

# for (tumor_to_select in c("worst", "representative")) {

# for (include_sample_without_chemo in c(T, F)) {

# for (keep_only_fetal_majority in c(T, F)) {

subfolder_name <- "multi_samples_par_patient"

if (keep_only_one_tumor_by_patient) {
  subfolder_name <- glue::glue("only_the_{gsub('representative', 'most_representative', tumor_to_select)}_tumor_per_patient")
} 
if (!include_sample_without_chemo) {
  subfolder_name <- paste0(subfolder_name, "_only_chemo")
} 
if (keep_only_fetal_majority) {
  subfolder_name <- paste0(subfolder_name, "_only_fetal_majority")
} 
if (keep_only_E_inf_30) {
  subfolder_name <- paste0(subfolder_name, "_only_E_inf_30")
} 


subfolder_name <- paste(series_to_work_on, subfolder_name, sep = "_")

print(subfolder_name)
resdir_parameters <- file.path(resdir, subfolder_name) %T>% 
  CreateDirIfNecessary()

# 1.3.2 annot to test -----------
genes_to_test <- c("CTNNB1", "NFE2L2", "TERT", "APC", "ARID1A", "AXIN1",
                   "IRF2", "CCND1", "MDM4", "RPS6KA3",
                   "GPC3", "CDKN1C")

fluidigm_to_test <- c("IGF2", "H19")

to_test <- c("Gender", "Type", "Histological.Diagnosis", "Age", "Age_at_diag","Sup_5_yrs","Sup_3_yrs", "Inf_3_yrs","Etiology", 
              "Death", "RC1_recode", # "RC1", "Last_Outcome"
             "F_by_Tumor_sup30", "E_by_Tumor_sup30", "M_by_Tumor_sup30", "MT_by_Tumor_sup20",
             "C_by_Tumor_sup30", "SCUD_by_Tumor_sup0", "FMA_by_Tumor_sup30", "FP_by_Tumor_sup30", 
             "Majority_by_Sample", "Liver_progenitor","Majority_by_Tumor",
             "F_by_Sample_sup30", "E_by_Sample_sup30", "M_by_Sample_sup30", "C_by_Sample_sup30",
             "FMA_by_Sample_sup30", "FP_by_Sample_sup30", "MT_by_Sample_sup30", "SCUD_by_Sample_sup30",
             paste0(rep(c("WISP.F", "WISP.E", "WISP.M", "WISP.F_sup50", "WISP.E_sup50", "WISP.M_sup50", "WISP_Majority"), each = 2),
                    c("", "_Fluidigm")),
             "Chemotherapy..yes.no.ND.", "Tumor.purity", "Vascular.Embols.invasion", "Vascular.Embols.viable",
             "nmut", "nmut.clonal", "nmut.subclonal", "cisplat.sig.group", "ROS.sig.group",
             paste0(fluidigm_to_test, "_fluidigm"),
             paste0(genes_to_test, "_yes_no"), "CTNNB1", 
             "second_drivers", 
             "CTNNB1_simple", "CTNNB1_simple_recode", "CTNNB1_recode",
             "AXIN1_recode",
             "locus_11p15_merged", "locus_11p15_recode_merged", "GOM_IC1_merged", "LOM_IC2_merged", "GOM_IC1_meca_recode_merged",
             "status_11p15_recode", "status_11p15",
             "RNAseq_clustering_3", "RNAseq_clustering_4", "RRBS_clustering", "C1C2.RNAseq", "C1.C2.fluidigm", 
             "immune_HB", "TLS", "TLS.max", "B.cat.N",
             "AFP_pre_chemo", "AFP_post_chemo", "log10_AFP_reduction", 
             "AFP_reduced_1.5_log", "AFP_pre_chemo_sup_1.000.000", "AFP_post_chemo_sup_10.000",
             "Response_to_neoadj_therapy", "response_to_neoadj_therapy_sd","reascension_AFP",
             "PRETEXT_Stade","PRETEXT_M", "PRETEXT_V","PRETEXT_P","PRETEXT_E","PRETEXT_R","PRETEXT_F", "PRETEXT_IV",
             "Refractory","Relapse_progression", 
             "POSTEXT_M","Viable_M_at_histology")

to_test_surv <- c("Gender", "Sup_5_yrs","Sup_3_yrs", "Inf_3_yrs","Age_at_diag",
                  "Majority_by_Tumor", "Majority_by_Sample", "F_by_Tumor_sup30","E_by_Tumor_sup30", "M_by_Tumor_sup30", 
                  "MT_by_Tumor_sup20", "E_by_Tumor_sup20",
                  "RC1_recode", #"Last_Outcome",
                  "SCUD_by_Tumor_sup0", "C_by_Tumor_sup30", "FMA_by_Tumor_sup30", "FP_by_Tumor_sup30", 
                  "Vascular.Embols.invasion", "Vascular.Embols.viable",
                  paste0(rep(c("WISP.F_sup50", "WISP.E_sup50", "WISP.M_sup50", "WISP_Majority"), each = 2),
                         c("", "_Fluidigm")),
                  "WISP.E_max_sup_50",
                  "Chemotherapy..yes.no.ND.", "cisplat.sig.group", "ROS.sig.group", "any_cisplat.sig.group",
                  paste0(genes_to_test, "_yes_no"), 
                  "locus_11p15_merged", "locus_11p15_recode_merged", "GOM_IC1_merged", "LOM_IC2_merged", "GOM_IC1_meca_recode_merged",
                  "status_11p15_recode", "status_11p15", 
                  "second_drivers",
                  "status_11p15_ext_recode", "status_IC1_recode",
                  "CTNNB1_simple", "CTNNB1_simple_recode","CTNNB1_recode",
                  "AXIN1_recode",
                  "RNAseq_clustering_3", "RNAseq_clustering_4","Liver_progenitor",
                  "RRBS_clustering", "C1C2.RNAseq", "C1.C2.fluidigm",
                  "immune_HB", "TLS", "TLS.max",
                  "Hepatocytic",
                  "B.cat.N", "AFP_reduced_1.5_log", "AFP_reduced_1_log", "AFP_pre_chemo_sup_1.000.000", 
                  "AFP_post_chemo_sup_10.000", "reascension_AFP",
                  "Response_to_neoadj_therapy", "response_to_neoadj_therapy_sd",
                  "PRETEXT_Stade","PRETEXT_M","PRETEXT_V","PRETEXT_P","PRETEXT_E","PRETEXT_R","PRETEXT_F", "PRETEXT_IV",
                  "Refractory","Relapse_progression", 
                  "POSTEXT_M","Viable_M_at_histology", "Death")
                  

mycol <- c("cold" = "blue", "mixed" = "orange", "hot" = "indianred", "mixed_hot" = "gold", # immune
           "yes" = "#B3381D", "no" = "black", # yes no features
           "yes (no viable cells)" = "blue", # yes / yes (no viable cells) / no 
           "FEW" = "lightblue", # B.cat.N  
           "Mut" = "#B3381D", "NM" = "black", # mutations
           "RC1" = "gold", "no_RC1" = "black", # RC1
           "cn-LOH" = "royalblue", "DEL" = "forestgreen", "DEL-LOH" = "skyblue3", "wt" = "black", # 11p15
           "GOM_IC1" = "royalblue", "LOM_IC2" = "skyblue3", # 11p15 recode
           "LOH" = "darkblue",  "no-LOH" = "grey90", # 11p15 recode
           "0" = "#FCFDBFFF","sup0" = "#FEAF77FF", "sup10" = "#F1605DFF","sup20" = "#B63679FF","sup30" = "#721F81FF", "sup60" = "#2D1160FF", "sup90" = "#000004FF", # histo features
           "F" = "#DE81C0","E" = "#500787","M" = "#A84325","SCUD" = "#0066E8", # histo majority
           "Hep" = "#A6CEE3", "H-hot" = "#A6CEE3","H-cold" = "#B7A7E3", "LP" = "#500787", # RNAseq cluster 3
           "F1" = "#A6CEE3","F2" = "#B7A7E3","old" = "#33A02C", # RNAseq 4 clusters and RRBS cluster
           # paste0("E", 1:9),viridis::viridis(9), # clustering
           "C1" = "forestgreen", "C2" = "red", # C1 C2
           "GOM_by_LOH" = "darkblue","GOM_no_LOH" = "blue", "no_GOM" = "grey", # GOM meca
           "alt" = "#B3381D", # 11p15 recode
           "no.TLS" = "grey90", "Aggregate" = "lightblue", "Primary.follicle" = "darkblue", "Secondary.follicle" = "black", # TLS max
           "Male" = "blue", "Female" = "gold", # male female
           "other" = "grey", # other in E vs other RNAseq clustering
           "Pretext_I" = "lightblue", "Pretext_II" = "darkblue", "Pretext_III" = "gold", "Pretext_IV" = "red", # Pretext
           "PR" = "green","SD" = "gold","PD" = "red", "no_neoadjchemo" = "blue", # response_chemo
           "relapse" = "gold","progression" = "red", #relapse
           "young" = "lightblue", "old" = "darkblue",
           "Refractory" = "#B3381D", #Refractory
           "Large_deletion" = "green", "Missense" = "gold", "Small_deletion" = "red", "deleted" = "red", #CTNNB1_simple
NA)



# 2. Load annot, split anapath columns & add AFP response -------------
# 2.1 Load annot table --------------

annot_full <- read_excel(file.path(data_dir,"2022_05_04_integrated_pediatric_table.xlsx"))   

# annot_full <- grep("integrated_pediatric_table.Rdata", list.files(data_dir),
#                    value = T) %T>% {print(.)} %>%  
#   first() %>% 
#   file.path(data_dir, .) %>% 
#   LoadRdata()      

patients_to_work_with <- annot_full %>% filter(To_keep_for_survival_analysis_2022 == "yes") %>% pull(Patient.identification) %>% unique()

length(unique(annot_full$Patient.identification))
setdiff(patients_to_work_with, annot_full$Patient.identification)

annot_full %>% count(cisplat.sig.group)



annot_cisplat_any <- annot_full %>% 
  filter(!is.na(cisplat.sig.group), Histological.Diagnosis == "HB") %>% 
  group_by(Patient.identification) %>% 
  summarise(any_cisplat.sig.group = ifelse(any(cisplat.sig.group == "yes"), "yes", "no"))



length(unique(annot_cisplat_any$Patient.identification))
setdiff(patients_to_work_with, annot_cisplat_any$Patient.identification)

annot <- annot_full %>% 
  left_join(annot_cisplat_any) %>% 
  rowwise() %>% 
  mutate(CHCID = gsub("CHC", "#", CHCID),
         To_keep_for_survival_analysis_worst = gsub("yes\\?", "yes", To_keep_for_survival_analysis_worst),
         To_keep_for_survival_analysis_representative = gsub("yes\\?", "yes", To_keep_for_survival_analysis_representative),
         AFP_pre_chemo = as.numeric(AFP.at.diagnosis..ng.ml.),
         AFP_pre_chemo= ifelse(AFP_pre_chemo == "NA", NA_character_, AFP_pre_chemo),
         AFP_post_chemo = as.numeric(AFP.at.surgery.post.chemotherapy..ng.ml.),
         AFP_post_chemo= ifelse(AFP_post_chemo == "NA", NA_character_, AFP_post_chemo),
         log10_AFP_reduction = log10(AFP_pre_chemo/AFP_post_chemo),
         AFP_reduced_1.5_log = ifelse(log10_AFP_reduction >= 1.5, "yes", "no"),
         AFP_reduced_1_log = ifelse(log10_AFP_reduction >= 1, "yes", "no"),
         reascension_AFP = ifelse(sup1_AFP_augmentation_durant_chimioneoadj == "NA", 
                                  NA_character_, sup1_AFP_augmentation_durant_chimioneoadj),
         Age = as.numeric(Age.at.surgery.recode.months),
         Age_at_diag = as.numeric(Age_at_diag),
         Sup_5_yrs = ifelse(Age_at_diag > 60, "yes", "no"),
         Sup_3_yrs = ifelse(Age_at_diag > 36, "yes", "no"),
         Inf_3_yrs = ifelse(Age_at_diag <= 36, "yes", "no"),
         Hepatocytic = ifelse(RNAseq_clustering_4 %in% c("F1", "F2"), RNAseq_clustering_4, NA_character_),
         Liver_progenitor = ifelse(RNAseq_clustering_4 == "E", "yes", "no"),
         GOM_IC1_merged = ifelse(GOM_IC1_merged == "na", NA_character_, GOM_IC1_merged),
         AFP_pre_chemo_sup_1.000.000 = ifelse(AFP_pre_chemo >= 10^6, "yes", "no"),
         AFP_post_chemo_sup_10.000 = ifelse(AFP_post_chemo >= 10^4, "yes", "no"),
         Response_to_neoadj_therapy = ifelse(Response_to_neoadj_therapy %in% c("PD", "PR", "SD", "no_neoadjchemo"), 
                                                    Response_to_neoadj_therapy, NA_character_),
         Refractory = ifelse(Refractory %in% c("Refractory", "no"),Refractory, NA_character_),
         Second_drivers = ifelse(second_drivers %in% c("yes", "no"),second_drivers, NA_character_),
         response_to_neoadj_therapy_sd = ifelse(Recode_Response_to_neoadj_therapy_Stable_disease %in% c("PD", "PR", "SD","no_neoadjchemo"), 
                                                          Recode_Response_to_neoadj_therapy_Stable_disease, NA_character_),
         PRETEXT_IV = ifelse(PRETEXT_Stade == "Pretext_IV", "yes", "no"),
         # RC1 = ifelse(Last_Outcome == "RC1", "RC1", "no_RC1"),
         RC1_recode = ifelse(Last_Outcome == "RC1", "yes", "no"),
         WISP_Majority = ifelse(is.na(WISP.F), NA_character_, c("F", "E", "M")[order(c(WISP.F, WISP.E, WISP.M))[3]]),
         WISP_Majority_Fluidigm = ifelse(is.na(WISP.F_Fluidigm), NA_character_, c("F", "E", "M")[order(c(WISP.F_Fluidigm, WISP.E_Fluidigm, WISP.M_Fluidigm))[3]]),
         F_by_Sample_sup30 = ifelse(F_by_Sample %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         E_by_Sample_sup30 = ifelse(E_by_Sample %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         M_by_Sample_sup30 = ifelse(M_by_Sample %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         FP_by_Sample_sup30 = ifelse(FP_by_Sample %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         FMA_by_Sample_sup30 = ifelse(FMA_by_Sample %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         C_by_Sample_sup30 = ifelse(C_by_Sample %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         MT_by_Sample_sup30 = ifelse(MT_by_Sample %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         SCUD_by_Sample_sup30 = ifelse(SCUD_by_Sample %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         F_by_Tumor_sup30 = ifelse(F_by_Tumor %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         E_by_Tumor_sup30 = ifelse(E_by_Tumor %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         E_by_Tumor_sup20 = ifelse(E_by_Tumor %in% paste0("sup", c(20, 30, 60, 90)), "yes", "no"),
         M_by_Tumor_sup30 = ifelse(M_by_Tumor %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         FP_by_Tumor_sup30 = ifelse(FP_by_Tumor %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         FMA_by_Tumor_sup30 = ifelse(FMA_by_Tumor %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         C_by_Tumor_sup30 = ifelse(C_by_Tumor %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         MT_by_Tumor_sup20 = ifelse(MT_by_Tumor %in% paste0("sup", c(20, 30, 60, 90)), "yes", "no"),
         SCUD_by_Tumor_sup0 = ifelse(SCUD_by_Tumor == "0", "no", "yes"),
         WISP.F_sup50 = ifelse(WISP.F > .5, "yes", "no"),
         WISP.E_sup50 = ifelse(WISP.E > .5, "yes", "no"),
         WISP.M_sup50 = ifelse(WISP.M > .5, "yes", "no"),
         WISP.F_sup50_Fluidigm = ifelse(WISP.F_Fluidigm > .5, "yes", "no"),
         WISP.E_sup50_Fluidigm = ifelse(WISP.E_Fluidigm > .5, "yes", "no"),
         WISP.M_sup50_Fluidigm = ifelse(WISP.M_Fluidigm > .5, "yes", "no"),
         Vascular.Embols.invasion = case_when(Vascular.Embols.invasion %in% c("yes (sans cell viable)", "Thrombose sans reliquat tumoral viable",
                                                                              "yes (no viable cell)", "Emboles non viables") ~ "yes (no viable cells)",
                                              Vascular.Embols.invasion == "yes" ~ "yes",
                                              Vascular.Embols.invasion == "no" ~ "no",
                                              T ~ NA_character_),
         Vascular.Embols.viable = ifelse(Vascular.Embols.invasion == "yes (no viable cells)", "no", Vascular.Embols.invasion),
         Gender = gsub("M", "Male", gsub("F", "Female", Gender)),
         TLS = gsub("TLS", "yes", gsub("no.TLS", "no", TLS)),
         Viable_M_at_histology = ifelse(Viable_M_at_histology == "?", NA_character_, Viable_M_at_histology),
         locus_11p15_recode = ifelse(grepl("LOH", locus_11p15), "LOH", "no-LOH"),
         locus_11p15_recode_merged = ifelse(grepl("LOH", locus_11p15_merged), "LOH", "no-LOH"),
         locus_11p15_recode_MLPA = ifelse(grepl("LOH", locus_11p15_MLPA), "LOH", "no-LOH"),
         CTNNB1_simple_recode = ifelse(CTNNB1_simple == "Large_deletion" | CTNNB1_simple == "Small_deletion" , 
                                       "deleted", CTNNB1_simple),
         CTNNB1_recode = ifelse(CTNNB1_simple_recode != "wt", "Alt", "wt"),
         AXIN1_recode = ifelse(AXIN1 != "wt", "Alt", "wt" ),
         GOM_IC1_meca_recode_merged = ifelse(GOM_IC1_merged == "yes", ifelse(locus_11p15_recode_merged == "LOH", "GOM_by_LOH", "GOM_no_LOH"),
                                             "no_GOM"),
         GOM_IC1_meca_recode_MLPA = ifelse(GOM_IC1_MLPA == "yes", ifelse(locus_11p15_recode_MLPA == "LOH", "GOM_by_LOH", "GOM_no_LOH"),
                                           "no_GOM"),
         GOM_IC1_meca_recode = ifelse(GOM_IC1 == "yes", ifelse(locus_11p15_recode == "LOH", "GOM_by_LOH", "GOM_no_LOH"),
                                      "no_GOM"),
         status_11p15_recode = ifelse(status_11p15 != "wt", 
                                      "alt",
                                      "wt"),
         status_11p15_ext_recode = ifelse(status_11p15_ext != "wt", 
                                          "alt",
                                          "wt"),
         status_IC1_recode = ifelse(status_IC1 == "wt" | status_IC1 == "GAIN", 
                                    "wt",
                                    "alt"),
         RNAseq_clustering_3 = sub("F", "Hep", sub("E", "LP", RNAseq_clustering_3)), 
         RNAseq_clustering_4 = sub("F1", "H-hot", sub("F2", "H-cold", sub("E", "LP", RNAseq_clustering_4))), 
         C1C2_E = ifelse(is.na(C1C2.RNAseq) | is.na(E_by_Sample),
                         NA_character_,
                         ifelse(C1C2.RNAseq == "C1", 
                                "C1",
                                ifelse(E_by_Sample %in% paste0("sup", c(30, 60, 90)), 
                                       "C2_E",
                                       "C2_noE")))) %>% 
  set_colnames(gsub(paste0("To_keep_for_survival_analysis_", tumor_to_select), "to_keep", colnames(.))) %>% 
  as.data.frame()


# 2. bis load fluidigm data --------------

fluidigm_2DDCt <- grep("^F.+xlsx", list.files(fluidigm_dir), value = T) %T>% print(.) %>% 
  first() %>% 
  file.path(fluidigm_dir, .) %>% 
  read_xlsx(skip = 7) %>% 
  select(CHCID, 
         H19_fluidigm = H19.Hs00399294_g1,
         IGF2_fluidigm = IGF2.Hs00171254_m1)

fluidigm_DDCt <- cbind(CHCID = gsub("CHC", "#", fluidigm_2DDCt[, 1]),
                       log2(fluidigm_2DDCt[, -1]))

annot <- annot %>% 
  left_join(fluidigm_DDCt)


# 2.2 reformat mutations -------------
# remove "NA" for TERT and recode mutations as M / NM

for(g in genes_to_test) {
  ind <- which(annot[, g] == "NA")
  if(length(ind)) annot[ind, g] <- NA
  annot[, str_interp("${g}_yes_no")] <- c("Mut", "NM")[as.numeric(annot[, g] == "wt") + 1]
}
table(annot$TERT)
table(annot$TERT_yes_no)
table(annot$CTNNB1_simple_recode)
table(annot$AXIN1_recode)
table(annot$CTNNB1_recode)

# 2.3 Restrict annot according to parameters chosen in 1.3.1 ------------------
print(subfolder_name)
annot <- annot %T>% GivePatientsSamples() %>% 
  filter(grepl("T", CHCID))  %T>% GivePatientsSamples() %>% # only tumors
  filter(grepl("HB", Histological.Diagnosis))  %T>% GivePatientsSamples() %>%  # only HB diag
  as.data.frame() 

length(unique(annot$Patient.identification))
  # setdiff(patients_to_work_with, annot$Patient.identification)

NGS_alias <- c("NGS", "NGS_Validation", "NGS_redondant")
validation_alias <- c("NGS_Validation", "NGS_exclu_Validation", "Validation", "Validation_redondant")
if (series_to_work_on == "NGS_and_validation") annot <- annot %>% filter(Serie_2022 %in% c(NGS_alias, validation_alias)) %T>% GivePatientsSamples() 
if (series_to_work_on == "NGS") annot <- annot %>% filter(Serie_2022 %in% NGS_alias) %T>% GivePatientsSamples() 
if (series_to_work_on == "Validation") annot <- annot %>% filter(Serie_2022 %in% validation_alias) %T>% GivePatientsSamples() 

length(unique(annot$Patient.identification))
setdiff(patients_to_work_with, annot$Patient.identification)

annot_integrated_multi <- annot %>% 
  filter(Type == "T", !is.na(WISP.E)) %>% 
  group_by(Patient.identification) %>% 
  summarise(WISP.E_max_patient = max(WISP.E)) %>% 
  mutate(WISP.E_max_sup_50 = ifelse(WISP.E_max_patient > .5, "yes", "no")) 
annot <- annot %>% 
  left_join(annot_integrated_multi)


if (keep_only_one_tumor_by_patient) annot <- annot %>% filter(To_keep_for_survival_analysis_2022 == "yes", CHCID != "#3513T") %T>% GivePatientsSamples() 
if (!include_sample_without_chemo) annot <- annot %>% filter(Chemotherapy..yes.no.ND. == "yes") %T>% GivePatientsSamples()
if (keep_only_fetal_majority) annot <- annot %>% filter(Majority_by_Tumor == "F") %T>% GivePatientsSamples() 
if (keep_only_E_inf_30) annot <- annot %>% filter(E_by_Tumor_sup30 == "no") %T>% GivePatientsSamples() 

length(unique(annot$Patient.identification))
setdiff(patients_to_work_with, annot$Patient.identification)

#table(annot$Chemotherapy..yes.no.ND)
#table(annot$To_keep_for_survival_analysis_2022)
#table(annot$Majority_by_Tumor)
#table(annot$Type) # only primary tumors thanks to the to_keep
#count(annot, Patient.identification) %>% .[["n"]] %>% table() # only 1 CHCID per patient thanks to the to_keep
#table(annot$WISP.E_sup50, annot$WISP.E_max_sup_50)
#table(annot$Refractory)
#table(annot$PRETEXT_Stade)
#table(annot$CHCID)
#table(annot$patients_to_work_with)

#setdiff(patients_to_work_with, annot$To_keep_for_survival_analysis_2022)

# Analysis -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --------------

# 5. Survival vs annot -----------
# 5.1 Restrict survival -----------
stop_survival_data <- 60 # choose how many months you want to represent
# subfolder_name_survival = paste(subfolder_name, stop_survival_data, sep = "_")

resdir_parameters_surv <- file.path(resdir_parameters, glue::glue("Survival_trunctaed_to_{stop_survival_data}_months")) %T>% 
  CreateDirIfNecessary()

annot_surv <- annot %>% 
  filter(!is.na(DSS.status),
         !is.na(PFS.status),
         !is.na(TTD)) %>% 
  mutate(OSS.status_truncated = ifelse(TTD > stop_survival_data,
                                       0,
                                       OSS.status),
         DSS.status_truncated = ifelse(TTD > stop_survival_data,
                                       0,
                                       DSS.status),
         TTD_truncated = ifelse(TTD > stop_survival_data,
                                stop_survival_data,
                                TTD),
         PFS.status_truncated = ifelse(TTP > stop_survival_data,
                                       0,
                                       PFS.status),
         TTP_truncated = ifelse(TTP > stop_survival_data,
                                stop_survival_data,
                                TTP)) %>% 
  as.data.frame() %T>% GivePatientsSamples()

# table(annot_surv$DSS.status_truncated)
# table(annot_surv$PFS.status_truncated)
# 
# table(annot_surv$TTP_truncated)

  # annot %>% count(AFP_reduced_1.5_log)


  #setdiff(annot_surv$Patient.identification,annot$Patient.identification)

# 5.2 Create survival objects -----------
annot_surv$DSS_object = with(annot_surv, 
                             Surv(TTD_truncated, DSS.status_truncated == 1))


annot_surv$OS_object = with(annot_surv, 
                            Surv(TTD_truncated, OSS.status_truncated == 1))


annot_surv$PFS_object = with(annot_surv, 
                             Surv(TTP_truncated , PFS.status_truncated == 1))
annot_surv

# table (annot_surv$TTP_truncated)
# table (annot_surv$PFS.status_truncated)
# table (annot_surv$PFS_object)
# 
# 
# table (annot_surv$TTD_truncated)
# table (annot_surv$OSS.status_truncated)
# table (annot_surv$OS_object)

  # PFS <- survfit(annot_surv$PFS_object ~ 1)    #Intervalle de confiance 0.95
  # plot(PFS, conf.int = TRUE, col = c("blue", "red"))

  # PFS99 <- survfit(annot_surv$PFS_object ~ 1, conf.int = 0.99) #Intervalle de confiance 0.99
  # plot(PFS99, conf.int = TRUE, col = c("blue", "red"))

  # setdiff(categ_survival, colnames(annot_surv))


# 5.3 Spot significant association -----------
surv_table <- data.frame(categ = "test", DSS = 1, PFS = 1, effectif_size = "100/100", big_enough = T) %>% factoall()
for (categ_survival in to_test_surv) {
  print(categ_survival)
  if (length(unique(na.omit(annot_surv[, categ_survival]))) < 2) {
    print("cannot compared, all are the same")
    surv_table <- rbind(surv_table, c(categ_survival, NA_real_, NA_real_, NA_real_, NA_real_))
    next
  } 
  
  surv_test_DSS <- survdiff(DSS_object ~ annot_surv[, categ_survival], data = annot_surv, rho = 0) # rho = 0 for logrank test
  log_rank_pval_DSS <- 1 - pchisq(surv_test_DSS$chisq, length(surv_test_DSS$n) - 1) 
  print(glue::glue("DSS p-value = {log_rank_pval_DSS}"))
  
  surv_test_PFS <- survdiff(PFS_object ~ annot_surv[, categ_survival], data = annot_surv, rho = 0) # rho = 0 for logrank test
  log_rank_pval_PFS <- 1 - pchisq(surv_test_PFS$chisq, length(surv_test_PFS$n) - 1)
  print(glue::glue("PFS p-value = {log_rank_pval_PFS}"))
  
  effectif_s <- table(annot_surv[, categ_survival])
  
  surv_table <- rbind(surv_table, 
                      c(categ_survival, log_rank_pval_DSS, log_rank_pval_PFS,
                        paste(effectif_s, collapse = "/"), sum(effectif_s > 3) > 1))
  
} # end of the loop on to_test_surv


surv_table <- surv_table %>% 
  rowwise() %>% 
  mutate(DSS = as.numeric(DSS),
         PFS = as.numeric(PFS),
         big_enough = big_enough == "TRUE",
         min_p_val = min(DSS, PFS)) %>% 
  arrange(min_p_val) %>% 
  as.data.frame()

WriteXLS::WriteXLS("surv_table",
                   ExcelFileName = file.path(resdir_parameters_surv,
                                             glue::glue("Survival_trunctaed_to_{stop_survival_data}_months.xlsx")),
                   AutoFilter = T, BoldHeaderRow = T, na = "", FreezeRow = 1, FreezeCol = 1)

surv_signif <- surv_table %>% 
  arrange(-big_enough) %>% 
  filter(min_p_val < .1) %>% 
  pull(categ)

# 5.4 Plot survival vs different annotations -----------

{
  pdf(file.path(resdir_parameters_surv, glue::glue("Survival_trunctaed_to_{stop_survival_data}_months.pdf")),
      width = 8, height = 8)
  
  for (categ_survival in setdiff(union(surv_signif, "RNAseq_clustering_4"), "test")) {
    print(categ_survival)
    annot_var = annot_surv %T>% GivePatientsSamples() %>% 
      .[!is.na(annot_surv[, categ_survival]), ] %T>% GivePatientsSamples() # remove samples where the variable is unknown or not part of the groups you're interested in
    
    myvar_surv_factors = sort(unique(annot_var[, categ_survival]))
    if (length(myvar_surv_factors) < 2) {
      print("cannot compared, all are the same")
      next
    } 
    
    DSS_my_categ = npsurv(DSS_object ~ annot_var[, categ_survival], data = annot_var)
    #OS_my_categ = npsurv(OS_object ~ annot_var[, categ_survival], data = annot_var)
    PFS_my_categ = npsurv(PFS_object ~ annot_var[, categ_survival], data = annot_var)
    
    myvar_surv_factors <- annot_var %>% filter(!is.na(DSS.status_truncated)) %>% pull(categ_survival) %>% unique() %>% sort()
    
    gg_DSS <- survminer::ggsurvplot(DSS_my_categ, risk.table = T, pval = T, pval.method = T, legend.title = categ_survival,
                                    xlab = "Survival time (months)", ylab = "Disease related survival", break.time.by = 12,
                                    palette = mycol[myvar_surv_factors], legend.labs = myvar_surv_factors)
    #gg_OS <- survminer::ggsurvplot(OS_my_categ, risk.table = T, pval = T, pval.method = T, legend.title = categ_survival,
    # xlab = "Survival time (months)", ylab = "Overall survival", break.time.by = 12,
    # palette = mycol[myvar_surv_factors], legend.labs = myvar_surv_factors)
    
    myvar_surv_factors <- annot_var %>% filter(!is.na(PFS.status_truncated)) %>% pull(categ_survival) %>% unique() %>% sort()
    
    gg_PFS <- survminer::ggsurvplot(PFS_my_categ, risk.table = T, pval = T, pval.method = T, legend.title = categ_survival,
                                    xlab = "Survival time (months)", ylab = "Progression free survival", break.time.by = 12,
                                    palette = mycol[myvar_surv_factors], legend.labs = myvar_surv_factors)
    print(gg_DSS) # for pdf
    print(gg_PFS) # for pdf
    png(file.path(resdir_parameters_surv, glue::glue("{categ_survival}_DSS.png")),
        width = 1000, height = 1000, res = 200)
    print(gg_DSS)
    dev.off() # dev off for png DSS
    png(file.path(resdir_parameters_surv, glue::glue("{categ_survival}_PFS.png")),
        width = 1000, height = 1000, res = 200)
    print(gg_PFS)
    dev.off() # dev off for png PFS
    
  } # end of the loop on significant parameters for survival
  
  
  dev.off() # dev off for pdf survival
}



dev.off()







#annot_var$tot <- "yes"

#PFS_my_categ = npsurv(PFS_object ~ annot_var [,"tot"], data = annot_var)

#gg_PFS <- survminer::ggsurvplot(PFS_my_categ, risk.table = T, pval = T, pval.method = T, legend.title = categ_survival,
# xlab = "Survival time (months)", ylab = "Progression free survival", break.time.by = 12)
 
                                
                             
                                   
                                
                                                               


OS <- survfit(annot_surv$OS_object ~ 1)   #Intervalle de confiance 0.95
plot(OS, conf.int = TRUE, col = c("blue", "red"))



OS99 <- survfit(annot_surv$OS_object ~ 1, conf.int = 0.99)  #Intervalle de confiance 0.99
plot(OS99, conf.int = TRUE, col = c("blue", "red"))


PFS <- survfit(annot_surv$PFS_object ~ 1)    #Intervalle de confiance 0.95
plot(PFS, conf.int = TRUE, col = c("blue", "red"))

PFS99 <- survfit(annot_surv$PFS_object ~ 1, conf.int = 0.99) #Intervalle de confiance 0.99
plot(PFS99, conf.int = TRUE, col = c("blue", "red"))



