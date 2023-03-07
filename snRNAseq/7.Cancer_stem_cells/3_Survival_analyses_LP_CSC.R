
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(ggpubr)
library(stringr)
library(survival)
library(rms) # for Kapplan Meier
library(magrittr)
library(dplyr)
library(readxl)
library(prismatic) # for color


geco.load <- function (filename) # File to be loaded
{
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error - function geco.load : file ", filename, 
            " not found"))
  NULL
}

factoall <- function (d, ncmax = 10) 
{
  n <- ncol(d)
  for (i in 1:n) {
    if (is.factor(d[, i])) {
      d[, i] <- as.character(d[, i])
      na <- which(is.na(d[, i]))
      num <- suppressWarnings(as.numeric(d[, i]))
      nanum <- which(is.na(num))
      if (length(nanum) == length(na)) {
        #                int <- suppressWarnings(as.integer(d[, i]))
        #                naint <- which(is.na(int))
        #                nc <- nchar(num)
        #                if (length(naint) == length(nanum) & all(nc < ncmax)) {
        #                  d[, i] <- int
        #                }
        #                else {
        d[, i] <- num
        #                }
      }
    }
  }
  d
}

GivePatientsSamples <- function(annot) {
  print(glue::glue("{length(annot$CHCID)} samples from {length(unique(annot$Patient.identification))} patients"))
}

######################
##### PARAMETERS #####
######################

### WARNING here everything is based on the study of the LP cluster: LP with CSC, LP without CSC, other tumors

datadir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/10.Cancer_stem_cells/Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/"
resdir <- file.path(datadir, "Survival_analysis_LP_csc")

### Which genes to use to separate LP tumors with/without CSC?
genes_for_LP_csc <- c("PROM1 and EPCAM", "PROM1 EPCAM CSC signature")[1] # "PROM1 and EPCAM" means bulk average expression of PROM1 and EPCAM, "PROM1 EPCAM CSC signature" means the genes overexpressed in LP CSC PROM1 and/or EPCAM in the differential expression analysis ~100 genes 
threshold_PROM1_EPCAM = 9 # select mean expression threshold to separate LP with/without CSC if we rely on PROM1 and EPCAM only
threshold_PROM1_EPCAM_csc_sig = 9.3 # select mean expression threshold to separate LP with/without CSC if we rely on the CSC signature otained by differential expression for PROM1 and EPCAM CSC

### Choose number of months to display on the survival plot
stop_survival_data <- 60 


###########################
##### 0. DATA LOADING #####
###########################

##### 1. Load mean bulk expression for each sample to distinguish LP with and without CSC #####
#----------------------------------------------------------------------------------------------

if(genes_for_LP_csc == "PROM1 and EPCAM"){ # if we rely only on PROM1 and EPCAM
  bulk_mean_exp <- geco.load(file.path(datadir, "PROM1_EPCAM_bulk_mean_exp.RData"))
  bulk_mean_exp <- bulk_mean_exp %>% filter(cluster != "NT") # keep only tumors 
  
  ### Create the 3 groups for survival analysis
  LP_with_CSC <- bulk_mean_exp[which(bulk_mean_exp$cluster=="LP" & bulk_mean_exp$mean_exp >= threshold_PROM1_EPCAM), "sample"]
  LP_without_CSC <- bulk_mean_exp[which(bulk_mean_exp$cluster=="LP" & bulk_mean_exp$mean_exp < threshold_PROM1_EPCAM), "sample"]
  other_tumors <- bulk_mean_exp[which(bulk_mean_exp$cluster!="LP"), "sample"]
  
}else if(genes_for_LP_csc == "PROM1 EPCAM CSC signature"){ # if we rely on the whole PROM1 and EPCAM csc signature (diff exp)
  bulk_mean_exp <- geco.load(file.path(datadir, "PROM1_EPCAM_csc_signature_bulk_mean_exp.RData"))
  bulk_mean_exp <- bulk_mean_exp %>% filter(cluster != "NT") # keep only tumors 
  
  ### Create the3 groups for survival analysis
  LP_with_CSC <- bulk_mean_exp[which(bulk_mean_exp$cluster=="LP" & bulk_mean_exp$mean_exp >= threshold_PROM1_EPCAM_csc_sig), "sample"]
  LP_without_CSC <- bulk_mean_exp[which(bulk_mean_exp$cluster=="LP" & bulk_mean_exp$mean_exp < threshold_PROM1_EPCAM_csc_sig), "sample"]
  other_tumors <- bulk_mean_exp[which(bulk_mean_exp$cluster!="LP"), "sample"]
}


##### 2. Load and prepare patient/sample bulk annotation table #####
#--------------------------------------------------------------------

annot <- geco.load("F:/Amelie-Datas/Annotations/2022_05_04_integrated_pediatric_table.RData"); dim(annot)


#############################################################
##### 1. PREPARE ANNOTATION TABLE FOR SURVIVAL ANALYSIS #####
#############################################################

### Restrict to HB tumors in the expression table
annot <- annot %>% filter(grepl("T", CHCID)) %>% filter(grepl("HB", Histological.Diagnosis)) %>% filter(CHCID %in% bulk_mean_exp$sample) %>% as.data.frame() # 100 samples

### Define patients to keep for survival analyses
patients_to_work_with <- annot %>% filter(To_keep_for_survival_analysis_2022 == "yes") %>% pull(Patient.identification) %>% unique(); length(patients_to_work_with) # 59 patients

### change colname of interest to "to_keep"
annot <-  annot %>% set_colnames(gsub("To_keep_for_survival_analysis_representative", "to_keep", colnames(.))) %>% 
  as.data.frame()

### Restrict to series of NGS and validation
NGS_alias <- c("NGS", "NGS_Validation", "NGS_redondant")
validation_alias <- c("NGS_Validation", "NGS_exclu_Validation", "Validation", "Validation_redondant")
annot <- annot %>% filter(Serie_2022 %in% c(NGS_alias, validation_alias)) %T>% GivePatientsSamples()  # 100 samples, 64 patients

### Keep only samples annotated for survival analysis
annot <- annot %>% filter(To_keep_for_survival_analysis_2022 == "yes", CHCID != "#3513T") %T>% GivePatientsSamples() # 59 samples, 59 patients

### add info of the 3 groups for survival analysis
annot$groups_for_survival_analysis <- ifelse(annot$CHCID %in% LP_with_CSC, "LP with CSC",
                                             ifelse(annot$CHCID %in% LP_without_CSC, "LP without CSC",
                                                    ifelse(annot$CHCID %in% other_tumors, "other tumors", "")))


################################
##### 2. SURVIVAL ANALYSIS #####
################################

##### 1. Prepare variables for survival analysis #####
#-----------------------------------------------------

categ_variables = list(groups_for_survival_analysis = list(variable_name = "groups_for_survival_analysis",
                                                  variable_factors = c("LP with CSC", "LP without CSC", "other tumors"),
                                                  variable_factors_names = c("LP with CSC", "LP without CSC", "other tumors"),
                                                  variable_factors_colors = c("red", "blue", "grey")))
categ <- categ_variables[[1]]
my_categ <- categ$variable_name


##### 2. Prepare survival annotation table #####
#-----------------------------------------------

annot$TTP <- as.numeric(annot$TTP) # seemingly this column was set to character...

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
  as.data.frame() %T>% GivePatientsSamples() # 59 samples


##### 3. Create survival objects #####
#-------------------------------------

annot_surv$DSS_object = with(annot_surv, 
                             Surv(TTD_truncated, DSS.status_truncated == 1))


annot_surv$OS_object = with(annot_surv, 
                            Surv(TTD_truncated, OSS.status_truncated == 1))


annot_surv$PFS_object = with(annot_surv, 
                             Surv(TTP_truncated , PFS.status_truncated == 1))
annot_surv


##### 4. Plot survival curves #####
#----------------------------------

categ_survival = my_categ

### Retrieve survival group names
myvar_surv_factors = sort(unique(annot_surv[, categ_survival]))

### Create objects comparing survival to the category of interest
DSS_my_categ = npsurv(DSS_object ~ annot_surv[, categ_survival], data = annot_surv)
PFS_my_categ = npsurv(PFS_object ~ annot_surv[, categ_survival], data = annot_surv)
myvar_surv_factors <- annot_surv %>% filter(!is.na(DSS.status_truncated)) %>% pull(categ_survival) %>% unique() %>% sort()

gg_DSS <- survminer::ggsurvplot(DSS_my_categ, pval = T, pval.method = T, legend.title = categ_survival,
                                  xlab = "Survival time (months)", ylab = "Disease related survival", break.time.by = 12,
                                  palette = c("LP with CSC" = "red", "LP without CSC" = "blue", "other tumors" = "grey"), legend.labs = myvar_surv_factors)
 
myvar_surv_factors <- annot_surv %>% filter(!is.na(PFS.status_truncated)) %>% pull(categ_survival) %>% unique() %>% sort()

gg_PFS <- survminer::ggsurvplot(PFS_my_categ, risk.table = T, pval = T, pval.method = T, legend.title = categ_survival,
                                  xlab = "Survival time (months)", ylab = "Progression free survival", break.time.by = 12,
                                  palette = c("LP with CSC" = "red", "LP without CSC" = "blue", "other tumors" = "grey"), legend.labs = myvar_surv_factors)
print(gg_DSS) 
print(gg_PFS) 
annot_surv %>% count(groups_for_survival_analysis, OSS.status)

png(file.path(resdir_parameters_surv, glue::glue("{categ_survival}_DSS.png")),
      width = 1000, height = 1000, res = 200)
print(gg_DSS)
dev.off() 
png(file.path(resdir_parameters_surv, glue::glue("{categ_survival}_PFS.png")),
    width = 1000, height = 1000, res = 200)
print(gg_PFS)
dev.off() 
  
















