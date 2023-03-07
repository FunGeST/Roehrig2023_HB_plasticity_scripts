

library(dplyr)
library(tidyr)
#library(Signac)
library(Seurat)
#library(ArchR)
library(ggplot2)
library(cowplot)
library(Rmagic)
library(SummarizedExperiment)

set.seed(1234)
memory.limit(30000000000)

geco.load <- function (filename) # File to be loaded
{
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error - function geco.load : file ", filename, 
            " not found"))
  NULL
}

geco.changeRange <- function (v, newmin = 1, newmax = 10) 
{
  oldmin <- min(v, na.rm = TRUE)
  oldmax <- max(v, na.rm = TRUE)
  newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin))
}


#addArchRGenome("hg38") # load the genome for each R session
# addArchRThreads(threads = 16) # takes automatically half of the threads available

#####################
##### FUNCTIONS #####
#####################

cor_function <- function(x, y){
  cor_test <- cor.test(x, y)
  return(cor_test$estimate)
}

cor_pval_function <- function(x, y){
  cor_test <- cor.test(x, y)
  return(cor_test$p.value)
}


######################
##### PARAMETERS #####
######################

name_samples <- c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3377T", "CHC3610T", "CHC3662T")
rna_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/"
rna_dir_imputation <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/3.Denoising/Imputation_in_house_tumor_cells/"
atac_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/3.snATACseq/Final_analyses/All_samples_merged/"
bulk_dir <- "F:/Amelie-Datas/Pediatric tumors/2.Expression/RNAseq/Without new NL samples/Supervised"

resdir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/4.Integration/2.Gene_regulatory_network/Final_GRN/"

# Define qvalue for significance
qval_sig=0.05

# Define thresholds for snRNAseq detection
#prop_cells = 0.005 # proportion of cells that must have >= nb-counts
#nb_counts = 1 # minimum nb of counts in prop_cells


###########################
##### 0. DATA LOADING #####
###########################

##### Trancription factors #####
#-------------------------------
TF_table <- geco.load("F:/Amelie-Datas/Data-hepatopartage/human.TF.rda") # 1639 TF identified by Lambert et al 2018
TF <- TF_table$external_gene_name

##### snRNA expression data #####
#--------------------------------
# Expression tables
rna_sample_all <- geco.load(file.path(rna_dir, "Seurat_object_analyses.RData"))
rna_exp <- as.matrix(GetAssayData(rna_sample_all, assay="RNA", slot="data")) # take RNA data assay (log normalized) for expression visualization
rna_counts <- as.matrix(GetAssayData(rna_sample_all, assay="RNA", slot="counts")) # take raw counts to identify which genes have sufficient expression
#keep_genes. <- apply(rna_counts, 1, function(x) length(x[x >= nb_counts]) >= (prop_cells * ncol(rna_counts))) # keep genes with sufficient expression
#keep_genes <- names(keep_genes.[which(keep_genes.==T)])
CHC3610T_cells <- grep("CHC3610T", colnames(rna_exp), value=T) # isolate name of mesenchymal cells to remove for some correlations (ex PC2 vs TF) 

# Imputated expression (cluster-based)
rna_exp_imputation <- geco.load(file.path(rna_dir_imputation, "exp_mat_clusters.RData")) 
imputation_clusters_annot <- geco.load(file.path(rna_dir_imputation, "cluster_annot.RData")); CHC3610T_clusters <- imputation_clusters_annot[grep("CHC3610T", imputation_clusters_annot$sample),"cluster"]
cells_imputation <- geco.load(file.path(rna_dir_imputation, "cells_in_clusters.RData")) 
keep_genes <- apply(rna_exp_imputation, 1, sum); keep_genes <- names(keep_genes[which(keep_genes!=0)]) # remove genes that have 0 everywgere in clusters

# PC1 and PC2 values
pc1_values <- rna_sample_all@reductions[["pca"]]@cell.embeddings[,"PC_1"] # extract PC1 values 
pc2_values <- rna_sample_all@reductions[["pca"]]@cell.embeddings[,"PC_2"] # extract PC2 values 

##### snATAC accessibility data #####
#------------------------------------
# TF motif deviations, cisbp v1
TF_motif_dev_v1 <- geco.load(file.path(atac_dir, "Motif_TF/MotifMatrix_cisbp.RData"))
TF_motif_dev_v1 <- TF_motif_dev_v1@assays@data[["z"]]
colnames(TF_motif_dev_v1) <- gsub("#", "_", gsub("-1", "", colnames(TF_motif_dev_v1))) # change to seurat nomenclature
rownames(TF_motif_dev_v1) <- paste0(gsub("\\_.*", "", rownames(TF_motif_dev_v1)), "_v1") # change nomenclature of TF names: remove the number after the underscore, add the version of cisbp
# compute cluster imputation for TF motif deviation (from snRNAseq clusters)
TF_motif_dev_imputation_v1 <- as.data.frame(matrix(NA, nrow=nrow(TF_motif_dev_v1), ncol=ncol(rna_exp_imputation))) 
rownames(TF_motif_dev_imputation_v1) <- rownames(TF_motif_dev_v1); colnames(TF_motif_dev_imputation_v1) <- colnames(rna_exp_imputation)
for(cluster in colnames(TF_motif_dev_imputation_v1)){
  cells_in_cluster <- intersect(names(cells_imputation[which(cells_imputation==cluster)]), colnames(TF_motif_dev_v1)) 
  TF_motif_dev_imputation_v1[,cluster] <- apply(TF_motif_dev_v1[,cells_in_cluster], 1, mean)
}
# TF motif deviations, cisbp v2
TF_motif_dev_v2 <- geco.load(file.path(atac_dir, "Motif_TF/MotifMatrix_cisbp_v2.00.RData"))
TF_motif_dev_v2 <- TF_motif_dev_v2@assays@data[["z"]]
colnames(TF_motif_dev_v2) <- gsub("#", "_", gsub("-1", "", colnames(TF_motif_dev_v2))) # change to seurat nomenclature
rownames(TF_motif_dev_v2) <- paste0(gsub("\\_.*", "", rownames(TF_motif_dev_v2)), "_v2") # change nomenclature of TF names: remove the number after the underscore, add the version of cisbp
# compute cluster imputation for TF motif deviation (from snRNAseq clusters)
TF_motif_dev_imputation_v2 <- as.data.frame(matrix(NA, nrow=nrow(TF_motif_dev_v2), ncol=ncol(rna_exp_imputation))) 
rownames(TF_motif_dev_imputation_v2) <- rownames(TF_motif_dev_v2); colnames(TF_motif_dev_imputation_v2) <- colnames(rna_exp_imputation)
for(cluster in colnames(TF_motif_dev_imputation_v2)){
  cells_in_cluster <- intersect(names(cells_imputation[which(cells_imputation==cluster)]), colnames(TF_motif_dev_v2)) 
  TF_motif_dev_imputation_v2[,cluster] <- apply(TF_motif_dev_v2[,cells_in_cluster], 1, mean)
}

# comprises 1244 TF in total (691 common v1/v2, 179 exclusive v1, 374 exclusive v2)
# ~8% of "TF" present in ArchR are not in the TF list of Lambert et al (but no sign of interesting TF in this difference)


##### Bulk expression data #####
#-------------------------------
bulk_counts <- geco.load("F:/Amelie-Datas/Data-hepatopartage/Expression_matrix_new_2020/count_matrix.Rdata") # load bulk counts to have genes
# Differential tables
diff_bulk_H_vs_LP <- geco.load(file.path(bulk_dir, "cluster_F_vs_cluster_E/all_resullts_limma.RData")); diff_bulk_H_vs_LP$Gene <- rownames(diff_bulk_H_vs_LP)
diff_bulk_H_vs_LP$logFC <- -diff_bulk_H_vs_LP$logFC # We want LP genes to have positive logFC
diff_bulk_M_vs_HB <- geco.load(file.path(bulk_dir, "cluster_M_vs_HB/all_resullts_limma.RData")); diff_bulk_M_vs_HB$Gene <- rownames(diff_bulk_M_vs_HB)


###########################################################################
##### 1. CREATE TABLE WITH BULK/sNRNASEQ/SNATACSEQ INFORMATION FOR TF #####
###########################################################################

df_TF <- data.frame(TF=TF, stringsAsFactors = F); rownames(df_TF) <- df_TF$TF

##### 1. Add Bulk RNA info #####
#-------------------------------
df_TF$bulk_M_vs_H_LP_qval <- df_TF$bulk_M_vs_H_LP_logFC <- df_TF$bulk_LP_vs_H_qval <- df_TF$bulk_LP_vs_H_logFC <-  NA 
df_TF[intersect(df_TF$TF, rownames(diff_bulk_H_vs_LP)),"bulk_LP_vs_H_logFC"] <- diff_bulk_H_vs_LP[intersect(df_TF$TF, rownames(diff_bulk_H_vs_LP)),"logFC"]
df_TF[intersect(df_TF$TF, rownames(diff_bulk_H_vs_LP)),"bulk_LP_vs_H_qval"] <- diff_bulk_H_vs_LP[intersect(df_TF$TF, rownames(diff_bulk_H_vs_LP)),"adj.P.Val"]
df_TF[intersect(df_TF$TF, rownames(diff_bulk_M_vs_HB)),"bulk_M_vs_H_LP_logFC"] <- diff_bulk_M_vs_HB[intersect(df_TF$TF, rownames(diff_bulk_M_vs_HB)),"logFC"]
df_TF[intersect(df_TF$TF, rownames(diff_bulk_M_vs_HB)),"bulk_M_vs_H_LP_qval"] <- diff_bulk_M_vs_HB[intersect(df_TF$TF, rownames(diff_bulk_M_vs_HB)),"adj.P.Val"]

df_TF <- df_TF[intersect(rownames(df_TF), keep_genes),] # restrict to genes without 0 everywhere in snRNAseq


##### 2. Add snRNAseq info #####
#-------------------------------
df_TF$snRNAseq_cor_exp_PC2_pval <- df_TF$snRNAseq_cor_exp_PC2 <- df_TF$snRNAseq_cor_exp_PC1_pval <- df_TF$snRNAseq_cor_exp_PC1 <- NA 

### PC1 vs TF expression
y1=imputation_clusters_annot$mean_PC1; names(y1) <- imputation_clusters_annot$cluster
df_TF[intersect(df_TF$TF, rownames(rna_exp_imputation)),"snRNAseq_cor_exp_PC1"] <- apply(rna_exp_imputation[intersect(df_TF$TF, rownames(rna_exp_imputation)),names(y1)], 1, function(x) cor_function(x, y=y1))
df_TF[intersect(df_TF$TF, rownames(rna_exp_imputation)),"snRNAseq_cor_exp_PC1_pval"] <- apply(rna_exp_imputation[intersect(df_TF$TF, rownames(rna_exp_imputation)),names(y1)], 1, function(x) cor_pval_function(x, y=y1))
df_TF$snRNAseq_cor_exp_PC1_qval <- p.adjust(df_TF$snRNAseq_cor_exp_PC1_pval, method="BH") # adjust the p-values
df_TF$snRNAseq_cor_exp_PC1_pval <- NULL

### PC2 vs TF expression
y2=imputation_clusters_annot$mean_PC2; names(y2) <- imputation_clusters_annot$cluster; y2 <- y2[which(!names(y2) %in% CHC3610T_clusters)] # remove mesenchymal cells for PC2
df_TF[intersect(df_TF$TF, rownames(rna_exp_imputation)),"snRNAseq_cor_exp_PC2"] <- apply(rna_exp_imputation[intersect(df_TF$TF, rownames(rna_exp_imputation)),names(y2)], 1, function(x) cor_function(x, y=y2))
df_TF[intersect(df_TF$TF, rownames(rna_exp_imputation)),"snRNAseq_cor_exp_PC2_pval"] <- apply(rna_exp_imputation[intersect(df_TF$TF, rownames(rna_exp_imputation)),names(y2)], 1, function(x) cor_pval_function(x, y=y2))
df_TF$snRNAseq_cor_exp_PC2_qval <- p.adjust(df_TF$snRNAseq_cor_exp_PC2_pval, method="BH") # adjust the p-values
df_TF$snRNAseq_cor_exp_PC2_pval <- NULL

# restrict to TF with snRNAseq info
df_TF <- df_TF[which(df_TF$TF %in% rownames(rna_exp)),] # 1514/1639 TF = 93%

##### 3. Add snATACseq info ####
#-------------------------------
df_TF$snATACseq_cor_dev_PC2_pval <- df_TF$snATACseq_cor_dev_PC2 <- df_TF$snATACseq_cor_dev_PC1_pval <- df_TF$snATACseq_cor_dev_PC1 <- df_TF$cisbp_version <- NA 

### PC1 vs TF motif deviation
y1=imputation_clusters_annot$mean_PC1; names(y1) <- imputation_clusters_annot$cluster
y2=imputation_clusters_annot$mean_PC2; names(y2) <- imputation_clusters_annot$cluster; y2 <- y2[which(!names(y2) %in% CHC3610T_clusters)] # remove mesenchymal cells for PC2
for(g in df_TF$TF){
  if(paste0(g, "_v1") %in% rownames(TF_motif_dev_imputation_v1) | paste0(g, "_v2") %in% rownames(TF_motif_dev_imputation_v2)){
    # Goal: identify which of cisbp v1 or v2 fits best to the snRNAseq info
    # to do so: keep version with best correlation between motif deviations and  gene expression
    if(paste0(g, "_v1") %in% rownames(TF_motif_dev_imputation_v1)){
      ATAC_X_RNA_v1 <- cor.test(as.numeric(TF_motif_dev_imputation_v1[paste0(g, "_v1"), colnames(rna_exp_imputation)]), as.numeric(rna_exp_imputation[g,]))$estimate
    }else{ATAC_X_RNA_v1 <- -10} # set to -10 and not NA because if we set NA the comparison will not be possible afterwards, while we know that -10 will never be the best version
    if(paste0(g, "_v2") %in% rownames(TF_motif_dev_imputation_v2)){
      ATAC_X_RNA_v2 <- cor.test(as.numeric(TF_motif_dev_imputation_v2[paste0(g, "_v2"), colnames(rna_exp_imputation)]), as.numeric(rna_exp_imputation[g,]))$estimate
    }else{ATAC_X_RNA_v2 <- -10}
    # keep the version that gives the highest positive correlation to snRNAseq
    # if there is only one version available, keep that version
    if(max(ATAC_X_RNA_v1, ATAC_X_RNA_v2, na.rm=T)==ATAC_X_RNA_v1){
      best_version <- TF_motif_dev_imputation_v1; version="_v1"
      df_TF[g, "snATACseq_cor_dev_PC1"] <- cor.test(as.numeric(best_version[paste0(g, version), names(y1)]), y1)$estimate
      df_TF[g, "snATACseq_cor_dev_PC1_pval"] <- cor.test(as.numeric(best_version[paste0(g, version), names(y1)]), y1)$p.value
      df_TF[g, "snATACseq_cor_dev_PC2"] <- cor.test(as.numeric(best_version[paste0(g, version), names(y2)]), y2)$estimate
      df_TF[g, "snATACseq_cor_dev_PC2_pval"] <- cor.test(as.numeric(best_version[paste0(g, version), names(y2)]), y2)$p.value
      df_TF[g, "cisbp_version"] <- "v1"
    }else if(max(ATAC_X_RNA_v1, ATAC_X_RNA_v2, na.rm=T)==ATAC_X_RNA_v2){
      best_version <- TF_motif_dev_imputation_v2; version="_v2"
      df_TF[g, "snATACseq_cor_dev_PC1"] <- cor.test(as.numeric(best_version[paste0(g, version), names(y1)]), y1)$estimate
      df_TF[g, "snATACseq_cor_dev_PC1_pval"] <- cor.test(as.numeric(best_version[paste0(g, version), names(y1)]), y1)$p.value
      df_TF[g, "snATACseq_cor_dev_PC2"] <- cor.test(as.numeric(best_version[paste0(g, version), names(y2)]), y2)$estimate
      df_TF[g, "snATACseq_cor_dev_PC2_pval"] <- cor.test(as.numeric(best_version[paste0(g, version), names(y2)]), y2)$p.value
      df_TF[g, "cisbp_version"] <- "v2"
    }
  }else{df_TF[g, "snATACseq_cor_dev_PC1"] <- df_TF[g, "snATACseq_cor_dev_PC1_pval"] <- df_TF[g, "snATACseq_cor_dev_PC2"] <- df_TF[g, "snATACseq_cor_dev_PC2_pval"] <- NA} # if no ATAC info: set NA
  ATAC_X_RNA_v1 <- ATAC_X_RNA_v2 <- NA # reset just for safety
}

df_TF$snATACseq_cor_dev_PC1_qval <- p.adjust(df_TF$snATACseq_cor_dev_PC1_pval, method="BH") # adjust the p-values
df_TF$snATACseq_cor_dev_PC1_pval<- NULL
df_TF$snATACseq_cor_dev_PC2_qval <- p.adjust(df_TF$snATACseq_cor_dev_PC2_pval, method="BH") # adjust the p-values
df_TF$snATACseq_cor_dev_PC2_pval <- NULL

##### 4. Save #####
#------------------

save(df_TF, file=file.path(resdir, "TF_characteristics_bulk_snRNAseq_snATACseq.RData"))

# here NA = no information 

#######################################
##### 2. COMPUTE COMPOSITE SCORES #####
#######################################

##### 1. Compute quantiles for each info #####
#---------------------------------------------
# We directly remove non significant values to compute percentiles and NAs due to absence of information are not considered to compute percentiles

### Bulk
bulk_LP_vs_H <- df_TF[which(df_TF$bulk_LP_vs_H_qval <= qval_sig), "bulk_LP_vs_H_logFC"]; names(bulk_LP_vs_H) <- df_TF[which(df_TF$bulk_LP_vs_H_qval <= qval_sig), "TF"]
df_TF$percentile_bulk_LP_vs_H <- NA; df_TF[names(bulk_LP_vs_H),"percentile_bulk_LP_vs_H"]  <- cut(bulk_LP_vs_H,quantile(bulk_LP_vs_H, seq(0, by=0.01), na.rm=T),include.lowest=TRUE,labels=FALSE)
bulk_M_vs_H_LP <- df_TF[which(df_TF$bulk_M_vs_H_LP_qval <= qval_sig), "bulk_M_vs_H_LP_logFC"]; names(bulk_M_vs_H_LP) <- df_TF[which(df_TF$bulk_M_vs_H_LP_qval <= qval_sig), "TF"]
df_TF$percentile_bulk_M_vs_H_LP <- NA; df_TF[names(bulk_M_vs_H_LP),"percentile_bulk_M_vs_H_LP"]  <- cut(bulk_M_vs_H_LP,quantile(bulk_M_vs_H_LP, seq(0, by=0.01), na.rm=T),include.lowest=TRUE,labels=FALSE)

### snRNAseq
snRNAseq_LP_vs_H <- df_TF[which(df_TF$snRNAseq_cor_exp_PC2_qval <= qval_sig), "snRNAseq_cor_exp_PC2"]; names(snRNAseq_LP_vs_H) <- df_TF[which(df_TF$snRNAseq_cor_exp_PC2_qval <= qval_sig), "TF"]
df_TF$percentile_snRNAseq_LP_vs_H  <- NA; df_TF[names(snRNAseq_LP_vs_H),"percentile_snRNAseq_LP_vs_H"]  <- cut(snRNAseq_LP_vs_H, quantile(snRNAseq_LP_vs_H, seq(0, by=0.01), na.rm=T),include.lowest=TRUE,labels=FALSE)
snRNAseq_M_vs_H_LP <- df_TF[which(df_TF$snRNAseq_cor_exp_PC1_qval <= qval_sig), "snRNAseq_cor_exp_PC1"]; names(snRNAseq_M_vs_H_LP) <- df_TF[which(df_TF$snRNAseq_cor_exp_PC1_qval <= qval_sig), "TF"]
df_TF$percentile_snRNAseq_M_vs_H_LP  <- NA; df_TF[names(snRNAseq_M_vs_H_LP),"percentile_snRNAseq_M_vs_H_LP"]  <- cut(snRNAseq_M_vs_H_LP, quantile(snRNAseq_M_vs_H_LP, seq(0, by=0.01), na.rm=T),include.lowest=TRUE,labels=FALSE)

### snATACseq
snATACseq_LP_vs_H <- df_TF[which(df_TF$snATACseq_cor_dev_PC2_qval <= qval_sig), "snATACseq_cor_dev_PC2"]; names(snATACseq_LP_vs_H) <- df_TF[which(df_TF$snATACseq_cor_dev_PC2_qval <= qval_sig), "TF"]
df_TF$percentile_snATACseq_LP_vs_H  <- NA; df_TF[names(snATACseq_LP_vs_H),"percentile_snATACseq_LP_vs_H"]  <- cut(snATACseq_LP_vs_H, quantile(snATACseq_LP_vs_H, seq(0, by=0.01), na.rm=T),include.lowest=TRUE,labels=FALSE)
snATACseq_M_vs_H_LP <- df_TF[which(df_TF$snATACseq_cor_dev_PC1_qval <= qval_sig), "snATACseq_cor_dev_PC1"]; names(snATACseq_M_vs_H_LP) <- df_TF[which(df_TF$snATACseq_cor_dev_PC1_qval <= qval_sig), "TF"]
df_TF$percentile_snATACseq_M_vs_H_LP  <- NA; df_TF[names(snATACseq_M_vs_H_LP),"percentile_snATACseq_M_vs_H_LP"]  <- cut(snATACseq_M_vs_H_LP, quantile(snATACseq_M_vs_H_LP, seq(0, by=0.01), na.rm=T),include.lowest=TRUE,labels=FALSE)


##### 3. Compute composite scores #####
#--------------------------------------
# sum the percentiles of bulk, snRNAseq and snATACseq. this has to be done for each subtype separately
# WARNING: the problem of summing is that when we have NA the score will be much lower while it could be very good marker in RNA but the motif is not available --> loss of important info!!
# to balance this problem: divide by the nb of values (3 if no NAs, 2 if 1 NAS, 1 if 2 NAs)
# REMINDER: here NAs correspond to missing info or non significant data

df_TF$composite_score_H_LP <- df_TF$composite_score_M  <- df_TF$composite_score_LP  <- df_TF$composite_score_H  <- NA

### H score ###
for(i in 1:nrow(df_TF)){
  # Extract percentiles necessary to compute composite score
  # WARNING here percentile=100 = max in LP and min in H --> we have to reverse 
  p_bulk <- 100-df_TF[i,"percentile_bulk_LP_vs_H"]+1
  p_snRNA <- 100-df_TF[i,"percentile_snRNAseq_LP_vs_H"]+1
  p_snATAC <- 100-df_TF[i,"percentile_snATACseq_LP_vs_H"]+1
  if(length(which(is.na(c(p_bulk, p_snRNA, p_snATAC))))>=2){df_TF[i, "composite_score_H"]  <- NA # if at least 2 missing values we can't keep the TF
  }else{df_TF[i, "composite_score_H"]  <- sum(p_bulk, p_snRNA, p_snATAC, na.rm=T)/(3-length(which(is.na(c(p_bulk, p_snRNA, p_snATAC))))) # we do not want to take NA as info --> we can't divide by 3 but by 3-(nb of NAs) to have accurate info
  } 
}

### LP score ###
for(i in 1:nrow(df_TF)){
  # Extract percentiles necessary to compute composite score
  p_bulk <- df_TF[i,"percentile_bulk_LP_vs_H"] 
  p_snRNA <- df_TF[i,"percentile_snRNAseq_LP_vs_H"]
  p_snATAC <- df_TF[i,"percentile_snATACseq_LP_vs_H"]
  if(length(which(is.na(c(p_bulk, p_snRNA, p_snATAC))))>=2){df_TF[i, "composite_score_LP"]  <- NA # if at least 2 missing values we can't keep the TF
  }else{df_TF[i, "composite_score_LP"]  <- sum(p_bulk, p_snRNA, p_snATAC, na.rm=T)/(3-length(which(is.na(c(p_bulk, p_snRNA, p_snATAC))))) # we do not want to take NA as info --> we can't divide by 3 but by 3-(nb of NAs) to have accurate info
  } 
}

### M score ###
for(i in 1:nrow(df_TF)){
  # Extract percentiles necessary to compute composite score
  p_bulk <- df_TF[i,"percentile_bulk_M_vs_H_LP"] 
  p_snRNA <- df_TF[i,"percentile_snRNAseq_M_vs_H_LP"]
  p_snATAC <- df_TF[i,"percentile_snATACseq_M_vs_H_LP"]
  if(length(which(is.na(c(p_bulk, p_snRNA, p_snATAC))))>=2){df_TF[i, "composite_score_M"]  <- NA # if at least 2 missing values we can't keep the TF
  }else{df_TF[i, "composite_score_M"]  <- sum(p_bulk, p_snRNA, p_snATAC, na.rm=T)/(3-length(which(is.na(c(p_bulk, p_snRNA, p_snATAC))))) # we do not want to take NA as info --> we can't divide by 3 but by 3-(nb of NAs) to have accurate info
  } 
}

### H/LP score ###
for(i in 1:nrow(df_TF)){
  # Extract percentiles necessary to compute composite score
  # WARNING here percentile=100 = max in M and min in H LP --> we have to reverse 
  p_bulk <- 100-df_TF[i,"percentile_bulk_M_vs_H_LP"]+1
  p_snRNA <- 100-df_TF[i,"percentile_snRNAseq_M_vs_H_LP"]+1
  p_snATAC <- 100-df_TF[i,"percentile_snATACseq_M_vs_H_LP"]+1
  if(length(which(is.na(c(p_bulk, p_snRNA, p_snATAC))))>=2){df_TF[i, "composite_score_H_LP"]  <- NA # if at least 2 missing values we can't keep the TF
  }else{df_TF[i, "composite_score_H_LP"]  <- sum(p_bulk, p_snRNA, p_snATAC, na.rm=T)/(3-length(which(is.na(c(p_bulk, p_snRNA, p_snATAC))))) # we do not want to take NA as info --> we can't divide by 3 but by 3-(nb of NAs) to have accurate info
  } 
}

##### 3. Save #####
#------------------

save(df_TF, file=file.path(resdir, "TF_characteristics_bulk_snRNAseq_snATACseq_composite_scores.RData"))
openxlsx::write.xlsx(df_TF, file.path(resdir, "TF_characteristics_bulk_snRNAseq_snATACseq_composite_scores.xlsx"))



