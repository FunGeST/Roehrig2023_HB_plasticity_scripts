
# strategy based on https://www.biorxiv.org/content/10.1101/2022.02.26.482041v1.full.pdf (The chromatin landscape of Th17cells reveals mechanisms of diversification of regulatory and pro inflammatory states, Thakore et al, 2022)

library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(corrplot)
library(dendextend)
set.seed(1234)

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

col_fun <- function(column, subtype){
  # Define color to make gradient
  if(subtype=="H"){col <- "#E8AAD6"}else if(subtype=="LP"){col <- "#510787"}else if(subtype=="M"){col <- "#6F3323"} 
  circlize::colorRamp2(c(min(column, na.rm=T), max(column, na.rm=T)), c("white", col))
}


######################
##### PARAMETERS #####
######################

rna_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/"
rna_dir_imputation <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/3.Denoising/Imputation_in_house_tumor_cells/"
archrdir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/3.snATACseq/Final_analyses/All_samples_merged/"
bulk_diff_dir <- "F:/Amelie-Datas/Pediatric tumors/2.Expression/RNAseq/Without new NL samples/Supervised"
bulk_dir <- "F:/Amelie-Datas/Data-hepatopartage/Expression_matrix_new_2020/"

resdir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/4.Integration/2.Gene_regulatory_network/Final_GRN"

# Define thresholds for snRNAseq detection
prop_cells = 0.01 # proportion of cells that must have >= nb-counts
nb_counts = 1 # minimum nb of counts in prop_cells

# Thresholds for peak2gene
FDR_threshold_peaksgene = 0.01 # for correlation > 0.4 all FDR are far below 0.01 anyway
varCutOffATAC= varCutOffRNA = 0.25 # default ArchR tfhreshold
cor_threshold=0.4

# Differential expression thresholds
FDR_threshold_bulk = 0.05 
logFC_threshold_bulk = log2(1.5) 


###########################
##### 0. DATA LOADING #####
###########################

### Multiome ###
# Expression tables
rna_sample_all <- geco.load(file.path(rna_dir, "Seurat_object_analyses.RData"))
rna_exp <- as.matrix(GetAssayData(rna_sample_all, assay="RNA", slot="data")) # take RNA data assay (log normalized) for expression visualization
rna_counts <- as.matrix(GetAssayData(rna_sample_all, assay="RNA", slot="counts")) # take raw counts to identify which genes have sufficient expression
CHC3610T_cells <- grep("CHC3610T", colnames(rna_exp), value=T) # isolate name of mesenchymal cells to remove for some correlations (ex PC2 vs TF) 
cells_H <- rownames(rna_sample_all@meta.data[which(rna_sample_all$H_LP_M_clean_groups_1=="H"),])
cells_LP <- rownames(rna_sample_all@meta.data[which(rna_sample_all$H_LP_M_clean_groups_1=="LP"),])
cells_M <- rownames(rna_sample_all@meta.data[which(rna_sample_all$H_LP_M_clean_groups_1=="M"),])
keep_genes. <- apply(rna_counts, 1, function(x) length(x[x >= nb_counts]) >= (prop_cells * ncol(rna_counts))) # keep genes with sufficient expression in H, HLP, LP, M tumor cells
keep_genes <- names(keep_genes.[which(keep_genes.==T)])

### TF defined by bulk/Multiome composite scores
df_TF_all <- geco.load(file.path(resdir, "TF_characteristics_bulk_snRNAseq_snATACseq_composite_scores.RData"))
df_TF <- df_TF_all[which(df_TF_all$TF %in% keep_genes),] # keep TF sufficiently detected in snRNA
df_TF <- df_TF[which((!is.na(df_TF$percentile_bulk_LP_vs_H) | !is.na(df_TF$percentile_bulk_M_vs_H_LP))),] # keep TF with significant bulk info, removes 31% TF
# Take top 20 TF for each subtype
top_TF_H <- df_TF %>% arrange(desc(composite_score_H)) %>% pull(TF) %>% head(20)
top_TF_LP <- df_TF %>% arrange(desc(composite_score_LP)) %>% pull(TF) %>% head(20)
top_TF_M <- df_TF %>% arrange(desc(composite_score_M)) %>% pull(TF) %>% head(20)
top_TF_H_LP <- df_TF %>% arrange(desc(composite_score_H_LP)) %>% pull(TF) %>% head(20) # we do not select the same nb of TF for each subtype because LP and M TF tend to have less motif info than H and H+LP
# adjust the initial nb of TF to have ~15 TF at the end after intersecting with the motif matrix
TF_GEPELIN <- c("HNF1A", "HNF4A", "MYCN", "MIXL1", "TWIST1", "TBX5") # Add the TF taken as example for the GEPELIN paper

# Imputated expression (cluster-based)
rna_exp_imputation <- geco.load(file.path(rna_dir_imputation, "exp_mat_clusters.RData")) 
imputation_clusters_annot <- geco.load(file.path(rna_dir_imputation, "cluster_annot.RData"))
cells_imputation <- geco.load(file.path(rna_dir_imputation, "cells_in_clusters.RData")) 

# Load ArchR motif match matrix
Motifs_in_peaks_v1_SE <- readRDS(file.path(archrdir, "/Annotations/Motif_cisbp-Matches-In-Peaks.rds")) # cisbp v1
r1 <- rowRanges(Motifs_in_peaks_v1_SE)# identify the row names
peak_names <- paste0(seqnames(r1), ":", start(r1), "-", end(r1))
Motifs_in_peaks_v1 <- Motifs_in_peaks_v1_SE@assays@data[["matches"]] # extract matrix part
rownames(Motifs_in_peaks_v1) <- peak_names
colnames(Motifs_in_peaks_v1) <- paste0(gsub("_.*", "", colnames(Motifs_in_peaks_v1)), "_v1")
Motifs_in_peaks_v2_SE <- readRDS(file.path(archrdir, "/Annotations/Motif_cisbp_v2.00-Matches-In-Peaks.rds")) # cisbp v2
r2 <- rowRanges(Motifs_in_peaks_v2_SE)# identify the row names
peak_names_v2 <- paste0(seqnames(r2), ":", start(r2), "-", end(r2))
Motifs_in_peaks_v2 <- Motifs_in_peaks_v2_SE@assays@data[["matches"]] # extract matrix part
rownames(Motifs_in_peaks_v2) <- peak_names_v2
colnames(Motifs_in_peaks_v2) <- paste0(gsub("_.*", "", colnames(Motifs_in_peaks_v2)), "_v2")
Motifs_in_peaks <- cbind(Motifs_in_peaks_v1, Motifs_in_peaks_v2) # combine cisbp v1 and v2
  
# Load TF deviation scores
TF_motif_dev_v1 <- geco.load(file.path(archrdir, "Motif_TF/MotifMatrix_cisbp.RData")) # cisbp v1
TF_motif_dev_v1 <- TF_motif_dev_v1@assays@data[["z"]]
colnames(TF_motif_dev_v1) <- gsub("#", "_", gsub("-1", "", colnames(TF_motif_dev_v1))) # change to seurat nomenclature
rownames(TF_motif_dev_v1) <- paste0(gsub("_.*", "", rownames(TF_motif_dev_v1)), "_v1")
TF_motif_dev_v2 <- geco.load(file.path(archrdir, "Motif_TF/MotifMatrix_cisbp_v2.00.RData")) # cisbp v2
TF_motif_dev_v2 <- TF_motif_dev_v2@assays@data[["z"]]
colnames(TF_motif_dev_v2) <- gsub("#", "_", gsub("-1", "", colnames(TF_motif_dev_v2))) # change to seurat nomenclature
rownames(TF_motif_dev_v2) <- paste0(gsub("_.*", "", rownames(TF_motif_dev_v2)), "_v2")
TF_motif_dev <- rbind(TF_motif_dev_v1, TF_motif_dev_v2) # combine cisbp v1 and v2

# Define top TF motifs
TF_with_ATAC_scores <- intersect(unique(c(top_TF_H, top_TF_LP, top_TF_M, top_TF_H_LP, TF_GEPELIN)), gsub("_.*", "", colnames(Motifs_in_peaks)))
TF_with_ATAC_scores <- intersect(TF_with_ATAC_scores, rownames(rna_exp))
top_TF_motifs <- c()
for(g in TF_with_ATAC_scores){ # add the corresponding version to the TF name
  if(df_TF_all[g, "cisbp_version"] == "v1"){top_TF_motifs <- c(top_TF_motifs, paste0(g, "_v1"))}
  if(df_TF_all[g, "cisbp_version"] == "v2"){top_TF_motifs <- c(top_TF_motifs, paste0(g, "_v2"))}
}

# Load peak2gene info
peaks2genes <- geco.load(file.path(archrdir, "Peaks_to_genes/peaks2genes_all_correlations_df.RData"))
peaks2genes_sig <- peaks2genes[which(peaks2genes$FDR <= FDR_threshold_peaksgene & peaks2genes$VarQATAC >= varCutOffATAC & peaks2genes$VarQRNA >= varCutOffRNA & peaks2genes$Correlation >= cor_threshold),]


###########################################
##### 1. IDENTIFY PEAKS WITH TF MOTIF #####
###########################################

##### 1. Look for TF motifs in peaks linked to genes #####
#---------------------------------------------------------
### Keep peaks linked to genes sufficiently expressed
#peaks2genes_sig <- peaks2genes_sig %>% filter(gene %in% keep_genes) # 8815 unique genes

### Search TF motifs in these peaks (for TF identified in 1st part)
Motifs_in_peaks_peak2genelink <- Motifs_in_peaks[which(rownames(Motifs_in_peaks) %in% peaks2genes_sig$peak), top_TF_motifs]
colnames(Motifs_in_peaks_peak2genelink) <- gsub("\\_.*", "", colnames(Motifs_in_peaks_peak2genelink)) # remove the cisbp version name in TF names 
a <- apply(Motifs_in_peaks_peak2genelink, 1, sum)
Motifs_in_peaks_peak2genelink <- Motifs_in_peaks_peak2genelink[names(a)[which(a>0)],] # keep peaks that contain the motifs of the regulator TF 

pdf(file.path(resdir, "Nb_TF_with_motif_in_peaks_peak2genelinkage.pdf"))
hist(a, breaks=40, main = "Nb of TF with motif in peaks from peak2genelinkage", xlab="Nb of TF motifs per peak", ylab="nb of peaks (peak2genelinkage)")
dev.off()

### Annotate TF motifs found in each peak
df_motifs_in_peaks <- data.frame(peaks=rownames(Motifs_in_peaks_peak2genelink), Motifs_in_peaks_peak2genelink=NA, stringsAsFactors = F)
for(i in 1:nrow(Motifs_in_peaks_peak2genelink)){ # for each peak
  df_motifs_in_peaks[i, "Motifs_in_peaks_peak2genelink"] <- paste(colnames(Motifs_in_peaks_peak2genelink)[which(Motifs_in_peaks_peak2genelink[i,]==1)], collapse=", ")
}

save(df_motifs_in_peaks, file=file.path(resdir, "H_LP_M/top_TF_motifs_in_peaks_linked_to_genes.RData"))


##################################################################################
##### 2. COMPUTE EXPRESSION CORRELATION BETWEEN TF AND GENES LINKED TO PEAKS #####
##################################################################################
# genes that are not linked to peaks are probably not related to the TF binding

##### 1. Create dataframe for expression correlations #####
#----------------------------------------------------------
TF_vs_target_genes_cor <- TF_vs_target_genes_pval <- matrix(NA, nrow=length(unique(peaks2genes_sig$gene)), ncol=length(TF_with_ATAC_scores))
colnames(TF_vs_target_genes_cor) <- colnames(TF_vs_target_genes_pval) <- TF_with_ATAC_scores; rownames(TF_vs_target_genes_cor) <- rownames(TF_vs_target_genes_pval) <- unique(peaks2genes_sig$gene)


##### 2.Compute correlation between TF and genes #####
#-------------------------------------------------------
# p values are adjusted for each TF over all potential target genes
for(TF in TF_with_ATAC_scores){
  print(TF)
  y1=as.numeric(rna_exp_imputation[TF,]) # restrict to common cells snRNA/snATAC
  names(y1) <- colnames(rna_exp_imputation)
  TF_vs_target_genes_cor[,TF] <- apply(rna_exp_imputation[rownames(TF_vs_target_genes_cor), names(y1)], 1, function(x) cor_function(x, y=y1))
  TF_vs_target_genes_pval[,TF] <-  p.adjust(apply(rna_exp_imputation[rownames(TF_vs_target_genes_pval), names(y1)], 1, function(x) cor_pval_function(x, y=y1)), method="BH") # adjust the p-values
}

save(TF_vs_target_genes_cor, file=file.path(resdir, "H_LP_M/Correlation_expression_TF_vs_target_genes.RData"))
save(TF_vs_target_genes_pval, file=file.path(resdir, "H_LP_M/Correlation_expression_TF_vs_target_genes_qval.RData"))


###############################################
##### 3. IDENTIFY THE NETWORK FOR EACH TF #####
###############################################
list_GRN <- list() # create empty list to contain GRN

##### 1. Restrict to the genes linked to peaks that contain the TF motif #####
#-----------------------------------------------------------------------------
for(TF in TF_with_ATAC_scores){
  b <- Motifs_in_peaks_peak2genelink[,TF]
  peaks_with_motif <- names(b)[which(b==TRUE)]
  potential_target_genes <- unique(peaks2genes_sig[which(peaks2genes_sig$peak %in% peaks_with_motif), "gene"])
  potential_target_genes <- intersect(potential_target_genes, keep_genes) # restrict to sufficiently detected genes
  list_GRN[[which(TF_with_ATAC_scores==TF)]] <- potential_target_genes
}
names(list_GRN) <- TF_with_ATAC_scores

##### 2. Remove TF from target genes #####
#-----------------------------------------
for(TF in TF_with_ATAC_scores){
  list_GRN[[TF]] <- list_GRN[[TF]][which(!list_GRN[[TF]] %in% names(list_GRN))]
}

##### 3. Restrict to the genes that are correlated with TF in snRNAseq #####
#---------------------------------------------------------------------------
for(TF in TF_with_ATAC_scores){
  # For whole visualization
 c <- TF_vs_target_genes_cor[list_GRN[[TF]],TF]
 d <- TF_vs_target_genes_pval[list_GRN[[TF]],TF]
 genes_cor_with_TF <- names(c[which(c>0.5)]) # keep genes with high correlation
 genes_cor_with_TF <- intersect(genes_cor_with_TF, names(d[which(d <= 0.001)])) # keep genes with significant correlation
 list_GRN[[TF]] <- intersect(genes_cor_with_TF, list_GRN[[TF]])
} # ~5000 genes

##### 4. Select targets for visualization #####
#----------------------------------------------
### Big  heatmap: restrict on nb of TF that detect the target to decrease nb of target genes
list_GRN_big <- list_GRN
for(g in unique(unlist(list_GRN_big))){ # for each gene in all possible targets in GRN
  if(length(grep(g, list_GRN_big)) < 5){ # if the gene is associated with > 5 TF
    for(TF in TF_with_ATAC_scores){
      list_GRN_big[[TF]] <- list_GRN_big[[TF]][which(list_GRN_big[[TF]] != g)]
    }
  }
} # ~ 2000 genes
### Short heatmap: restrict on the nb of targets per module to decrease nb of target genes
list_GRN_short1 <- list_GRN_short2 <- list_GRN
for(TF in TF_with_ATAC_scores){
  c <- TF_vs_target_genes_cor[list_GRN[[TF]],TF] # the targets linked to TF are already filtered for correlation and significance
  c <- c[which(names(c) %in% keep_genes)] # restrict to sufficiently expressed
  list_GRN_short1[[TF]] <- intersect(names(c[order(c, decreasing=T)][1:2]), list_GRN_short1[[TF]]) # Based on nb of targets per TF
  list_GRN_short2[[TF]] <- intersect(names(c[order(c, decreasing=T)][1:round(0.01*length(c))]), list_GRN_short2[[TF]]) # based on % of targets per TF
} # approximately 100-200 genes for visualization


##### 5. Save lists #####
#------------------------
save(list_GRN, file=file.path(resdir, "H_LP_M/GRN_list.RData"))
save(list_GRN_big, file=file.path(resdir, "H_LP_M/GRN_list_visu_big_heatmap.RData"))
save(list_GRN_short1, file=file.path(resdir, "H_LP_M/GRN_list_for_heatmap_top_2_targets_per_TF.RData"))
save(list_GRN_short2, file=file.path(resdir, "H_LP_M/GRN_list_for_heatmap_top_1_percent_targets_per_TF.RData"))

##### 6. Display nb of targets per TF #####
#------------------------------------------
pdf(file.path(resdir, "H_LP_M/Nb_of_target_genes_per_TF.pdf"))
nb_targets_per_TF <- rep(NA, times=length(TF_with_ATAC_scores)); names(nb_targets_per_TF) <- TF_with_ATAC_scores
for(TF in names(nb_targets_per_TF)){nb_targets_per_TF[TF] <- length(list_GRN[[TF]])}
barplot(nb_targets_per_TF[order(nb_targets_per_TF)], horiz = T, las=2, xlab="Nb of targets genes per TF", cex.names=0.4) # visualize results

nb_targets_per_TF <- rep(NA, times=length(TF_with_ATAC_scores)); names(nb_targets_per_TF) <- TF_with_ATAC_scores
for(TF in names(nb_targets_per_TF)){nb_targets_per_TF[TF] <- length(list_GRN_short1[[TF]])}
barplot(nb_targets_per_TF[order(nb_targets_per_TF)], horiz = T, las=2, xlab="Nb of targets genes per TF (visualization nb targets)", cex.names=0.4) # visualize results

nb_targets_per_TF <- rep(NA, times=length(TF_with_ATAC_scores)); names(nb_targets_per_TF) <- TF_with_ATAC_scores
for(TF in names(nb_targets_per_TF)){nb_targets_per_TF[TF] <- length(list_GRN_short2[[TF]])}
barplot(nb_targets_per_TF[order(nb_targets_per_TF)], horiz = T, las=2, xlab="Nb of targets genes per TF (visualization % targets)", cex.names=0.4) # visualize results
dev.off()


######################################################
##### 4. TF/TARGET GENES MODULES WITH CLUSTERING #####
######################################################
# Display snRNAseq correlations between TF and target gene

list_to_use_name = c("list_GRN_big", "list_GRN_short1", "list_GRN_short2")[2]
if(list_to_use_name=="list_GRN_big"){list_to_use <- list_GRN_big; heatmap_name = "cor target/TF >0.5, targets in >=5 TF"; file_name = "big"}
if(list_to_use_name=="list_GRN_short1"){list_to_use <- list_GRN_short1; heatmap_name = "top 2 targets per TF"; file_name = "short_nb_targets"}
if(list_to_use_name=="list_GRN_short2"){list_to_use <- list_GRN_short2; heatmap_name = "top 1% targets per TF"; file_name = "short_percent_targets"}
length(unique(unlist(list_to_use))) # how many genes in total


### Determine text size for labels
size_col_names=8; size_row_names=8; size_annot_labels=6; width_annot=0.2

##### 1. Prepare table for heatmap #####
#---------------------------------------
mat <- as.data.frame(TF_vs_target_genes_cor[unique(unlist(list_to_use)), names(list_to_use)])
mat <- mat[which(!rownames(mat) %in% colnames(mat)), ] # Restrict to target genes that are not among the TF under study
# keep all correlations because the non significant ones will be near 0 anyway

##### 2. Heatmap row annotations (TF) #####
#------------------------------------------
### Compute mean on H/LP/M subtype for each gene
row_annot <- as.data.frame(matrix(NA, nrow=ncol(mat), ncol=6)); rownames(row_annot) <- colnames(mat); colnames(row_annot) <- c("expH_sn", "expLP_sn", "expM_sn", "TFdevH_sn", "TFdevLP_sn", "TFdevM_sn")
row_annot$expH_sn <- apply(rna_exp[match(colnames(mat), rownames(rna_exp)),cells_H], 1, mean)
row_annot$expLP_sn <- apply(rna_exp[match(colnames(mat), rownames(rna_exp)),cells_LP], 1, mean)
row_annot$expM_sn <- apply(rna_exp[match(colnames(mat), rownames(rna_exp)),cells_M], 1, mean)
row_annot$TFdevH_sn <- apply(TF_motif_dev[top_TF_motifs, intersect(cells_H, colnames(TF_motif_dev))], 1, mean)
row_annot$TFdevLP_sn <- apply(TF_motif_dev[top_TF_motifs,intersect(cells_LP, colnames(TF_motif_dev))], 1, mean)
row_annot$TFdevM_sn <- apply(TF_motif_dev[top_TF_motifs,intersect(cells_M, colnames(TF_motif_dev))], 1, mean)

### set values between 0 and 1 to compare H/LP/M subtype for genes
for(i in 1:nrow(row_annot)){ 
  row_annot[i,c("expH_sn", "expLP_sn", "expM_sn")] <- geco.changeRange(row_annot[i,c("expH_sn", "expLP_sn", "expM_sn")], newmin=0, newmax=1)
} # no change in range for ATAC because those are scores that can be negative or positive where 0 = not enriched, not expression between 0 and X

### Create annotation
# For ATAC we only show the positive TF deviation scores 
row_ha = rowAnnotation(expH_sn =  row_annot$expH_sn, expLP_sn = row_annot$expLP_sn, expM_sn = row_annot$expM_sn, 
                       TFdevH_sn = row_annot$TFdevH_sn, TFdevLP_sn = row_annot$TFdevLP_sn, TFdevM_sn = row_annot$TFdevM_sn, 
                       col=list(expH_sn=col_fun(row_annot$expH_sn, "H"), expLP_sn=col_fun(row_annot$expLP_sn, "LP"), expM_sn=col_fun(row_annot$expM_sn, "M"), 
                                TFdevH_sn=circlize::colorRamp2(c(min(row_annot$TFdevH_sn, na.rm=T), 0, max(row_annot$TFdevH_sn, na.rm=T)), c("white", "white", "#E8AAD6")), 
                                TFdevLP_sn=circlize::colorRamp2(c(min(row_annot$TFdevLP_sn, na.rm=T), 0, max(row_annot$TFdevLP_sn, na.rm=T)), c("white", "white", "#510787")), 
                                TFdevM_sn=circlize::colorRamp2(c(min(row_annot$TFdevM_sn, na.rm=T), 0, max(row_annot$TFdevM_sn, na.rm=T)), c("white", "white", "#6F3323"))), 
                       na_col = "grey50", simple_anno_size=unit(width_annot, "cm"), annotation_name_gp= gpar(fontsize = size_annot_labels))


##### 3. Heatmap column annotations (target genes) ##### TO DO
#-------------------------------------------------------
### Compute mean on H/LP/M subtype for each gene
col_annot <- as.data.frame(matrix(NA, nrow=nrow(mat), ncol=3)); rownames(col_annot) <- rownames(mat); colnames(col_annot) <- c("expH_sn", "expLP_sn", "expM_sn")
col_annot$expH_sn <- apply(rna_exp[match(rownames(mat), rownames(rna_exp)),cells_H], 1, mean)
col_annot$expLP_sn <- apply(rna_exp[match(rownames(mat), rownames(rna_exp)),cells_LP], 1, mean)
col_annot$expM_sn <- apply(rna_exp[match(rownames(mat), rownames(rna_exp)),cells_M], 1, mean)

### set values between 0 and 1 to compare H/LP/M subtype for genes
for(i in 1:nrow(col_annot)){ 
  col_annot[i,c("expH_sn", "expLP_sn", "expM_sn")] <- geco.changeRange(col_annot[i,c("expH_sn", "expLP_sn", "expM_sn")], newmin=0, newmax=1)
} 

### Create annotation
col_ha = HeatmapAnnotation(expH_sn =  col_annot$expH_sn, expLP_sn = col_annot$expLP_sn, expM_sn = col_annot$expM_sn, 
                           col=list(expH_sn=col_fun(col_annot$expH_sn, "H"), expLP_sn=col_fun(col_annot$expLP_sn, "LP"), expM_sn=col_fun(col_annot$expM_sn, "M"), 
                                    na_col = "grey50", simple_anno_size=unit(width_annot, "cm"), annotation_name_gp= gpar(fontsize = size_annot_labels)))

##### 3. Heatmap #####
#---------------------
ht <- draw(ComplexHeatmap::Heatmap(t(mat), show_row_names=T, show_column_names = T, column_dend_height = unit(4, "cm"), clustering_distance_columns = "pearson", clustering_method_columns = "ward.D", 
                                   clustering_distance_rows = "pearson", clustering_method_rows = "ward.D", 
                                   column_title = paste(heatmap_name, "\n", "Target genes"), column_names_gp = gpar(fontsize = size_col_names), row_names_gp = gpar(fontsize = size_row_names), right_annotation = row_ha, top_annotation=col_ha,
                                   heatmap_legend_param = list(legend_direction="horizontal"), row_title="TF", na_col = "white"), 
           heatmap_legend_side="bottom")

# Re order clusters for nice visualization
if(list_to_use_name == "list_GRN_short1"){
  TF_dendro <- row_dend(ht); TF_current.order <- cutree(TF_dendro,k = 4)
  TF_new.order <- c(TF_current.order[c(which(TF_current.order==2), which(TF_current.order==4), which(TF_current.order==3), which(TF_current.order==1))]); TF_dendro.reordered <- rotate(TF_dendro,order = names(TF_new.order))
  targets_dendro <- column_dend(ht); targets_current.order <- cutree(targets_dendro,k = 5)
  targets_new.order <- c(targets_current.order[c(which(targets_current.order==2), which(targets_current.order==5), which(targets_current.order==4), which(targets_current.order==3), which(targets_current.order==1))])
  targets_dendro.reordered <- rotate(targets_dendro,order = names(targets_new.order))

ht_ordered <- draw(ComplexHeatmap::Heatmap(t(mat), show_row_names=T, show_column_names = T, column_dend_height = unit(4, "cm"), cluster_rows=TF_dendro.reordered, cluster_columns = targets_dendro.reordered,
                                   column_title = paste(heatmap_name, ", clusters reordered", "\n", "Target genes"), column_names_gp = gpar(fontsize = size_col_names), row_names_gp = gpar(fontsize = size_row_names), right_annotation = row_ha, top_annotation=col_ha,
                                   heatmap_legend_param = list(legend_direction="horizontal"), row_title="TF", na_col = "white"), heatmap_legend_side="bottom")
}

save(TF_new.order, file=file.path(resdir, paste0("H_LP_M/Heatmap_TF_targets_", file_name, "_TF_ordered.RData")))
save(targets_new.order, file=file.path(resdir, paste0("H_LP_M/Heatmap_TF_targets_", file_name, "_targets_ordered.RData")))


##### 4. Map of targets with TF binding site #####
#-------------------------------------------------
# The previous GRN heatmap is based on snRNAseq for clustering but that does not mean that all target genes in the module have the corresponding TF binding sites in the related enhancer peaks

### Recover TF and target order in heatmap
TF_order <- colnames(mat[,row_order(ht)]) 
targets_order <- rownames(mat[column_order(ht),])

mat_TFBS <- matrix(0, nrow=nrow(mat), ncol=ncol(mat)); colnames(mat_TFBS) <-TF_order; rownames(mat_TFBS) <- targets_order # initialize matrix of TF binding site presence in target genes
for(TF in colnames(mat_TFBS)){
  mat_TFBS[intersect(names(mat_TFBS[,TF]), list_GRN[[TF]]),TF] <- 1 # set 1 for target genes of the snRNAseq correlation module that contain the TF binding motif in their related enhancer peaks (peak2genelinkage)
}

if(list_to_use_name == "list_GRN_short1"){
  TF_order_reordered <- colnames(mat[,row_order(ht_ordered)]) 
  targets_order_reordered <- rownames(mat[column_order(ht_ordered),])
  
  mat_TFBS_reordered <- matrix(0, nrow=nrow(mat), ncol=ncol(mat)); colnames(mat_TFBS_reordered) <-TF_order_reordered; rownames(mat_TFBS_reordered) <- targets_order_reordered # initialize matrix of TF binding site presence in target genes
  for(TF in colnames(mat_TFBS_reordered)){
    mat_TFBS_reordered[intersect(names(mat_TFBS_reordered[,TF]), list_GRN[[TF]]),TF] <- 1 # set 1 for target genes of the snRNAseq correlation module that contain the TF binding motif in their related enhancer peaks (peak2genelinkage)
  }
}


ht2 <- draw(ComplexHeatmap::Heatmap(t(mat_TFBS), show_row_names=T, show_column_names = T, cluster_rows = F, cluster_columns = F, column_dend_height = unit(4, "cm"), col=c("white", "darkblue"),
                                    column_names_gp = gpar(fontsize = size_col_names), row_names_gp = gpar(fontsize = size_row_names), 
                                    heatmap_legend_param = list(legend_direction="horizontal"), na_col = "white"), heatmap_legend_side="bottom")

if(list_to_use_name == "list_GRN_short1"){
ht2_ordered <- draw(ComplexHeatmap::Heatmap(t(mat_TFBS_reordered), show_row_names=T, show_column_names = T, cluster_rows = F, cluster_columns = F, column_dend_height = unit(4, "cm"), col=c("white", "darkblue"),
                                    column_names_gp = gpar(fontsize = size_col_names), row_names_gp = gpar(fontsize = size_row_names), 
                                    heatmap_legend_param = list(legend_direction="horizontal"), na_col = "white"), heatmap_legend_side="bottom")
}

save(TF_order_reordered, file=file.path(resdir, paste0("H_LP_M/Heatmap_TF_targets_", file_name, "_TF_ordered.RData")))
save(targets_order_reordered, file=file.path(resdir, paste0("H_LP_M/Heatmap_TF_targets_", file_name, "_targets_ordered.RData")))


##### 5. Save the heatmaps #####
#-------------------------------

pdf(file.path(resdir, paste0("H_LP_M/Heatmap_TF_targets_", file_name, ".pdf")), width=14, height=9)
print(ht)
print(ht2)
print(ht_ordered)
print(ht2_ordered)
dev.off()

png(file.path(resdir, paste0("H_LP_M/Heatmap_TF_targets_", file_name, ".png")), width=2600, height=2100, res=200)
print(ht)
dev.off()
png(file.path(resdir, paste0("H_LP_M/Heatmap_TF_targets_", file_name, "_motifs.png")), width=2600, height=2100, res=200)
print(ht2)
dev.off()
if(list_to_use_name == "list_GRN_short1"){
png(file.path(resdir, paste0("H_LP_M/Heatmap_TF_targets_", file_name, "_ordered.png")), width=2600, height=2100, res=200)
print(ht_ordered)
dev.off()
png(file.path(resdir, paste0("H_LP_M/Heatmap_TF_targets_", file_name, "_motifs_ordered.png")), width=2600, height=2100, res=200)
print(ht2_ordered)
dev.off()
}



################################################################
##### 5. MEAN EXPRESSION IN CELL POPULATIONS FOR ALL GENES #####
################################################################

df_exp_H_LP_M <- data.frame(gene=unique(c(TF_with_ATAC_scores, unique(peaks2genes_sig$gene))), mean_exp_H=NA, mean_exp_LP=NA, mean_exp_M=NA, stringsAsFactors = F)
if(length(setdiff(df_exp_H_LP_M$gene, rownames(rna_exp)))==0){ # check that all genes are in expression matrix (should be)
  df_exp_H_LP_M$mean_exp_H <- apply(rna_exp[df_exp_H_LP_M$gene, cells_H], 1, mean)
  df_exp_H_LP_M$mean_exp_LP <- apply(rna_exp[df_exp_H_LP_M$gene, cells_LP], 1, mean)
  df_exp_H_LP_M$mean_exp_M <- apply(rna_exp[df_exp_H_LP_M$gene, cells_M], 1, mean)
}

write.table(df_exp_H_LP_M, file.path(resdir, "H_LP_M/mean_exp_H_LP_M_per_gene.csv"), sep=";", col.names=T, row.names=F)
