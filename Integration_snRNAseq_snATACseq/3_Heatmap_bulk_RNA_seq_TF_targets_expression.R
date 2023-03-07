
# strategy based on https://www.biorxiv.org/content/10.1101/2022.02.26.482041v1.full.pdf (The chromatin landscape of Th17cells reveals mechanisms of diversification of regulatory and pro inflammatory states, Thakore et al, 2022)

### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ArchR)
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


######################
##### PARAMETERS #####
######################

rna_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/"
bulk_dir <- "F:/Amelie-Datas/Data-hepatopartage/Expression_matrix_new_2020/"
grn_dir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/4.Integration/2.Gene_regulatory_network/Final_GRN"

resdir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/4.Integration/2.Gene_regulatory_network/Final_GRN/H_LP_M"


###########################
##### 0. DATA LOADING #####
###########################

### Load Multiome GRN
TF_ordered <- geco.load(file.path(resdir, "Heatmaps_TF/TF_ordered_by_PC2_position.RData")) # load TF ordered by position on PC2 using bins of 100 cells
list_GRN <- geco.load(file.path(resdir, "GRN_list.RData")) # list of targets per TF

### Load single cell markers
sample_seurat <- geco.load(file.path(rna_dir, "Seurat_object_analyses.RData")) # load object with all tumor cells from all samples to get PC2 values
pc2_values <- sample_seurat@reductions$pca@feature.loadings[,"PC_2"] # load pc2 values
#Define H and LP markers from PC2 by takin top 120 genes of PC2 extremities
LP_markers_sn <- names(sort(pc2_values, decreasing = T)[1:50])
H_markers_sn <- names(sort(pc2_values, decreasing = F)[1:50]) 

### Load bulk expression
bulk_exp <- geco.load("F:/Amelie-Datas/Data-hepatopartage/Expression_matrix_new_2020/exp.Rdata")
bulk_order <- geco.load("F:/Amelie-Datas/Pediatric tumors/2.Expression/RNAseq/Without new NL samples/Unsupervised/papier_GEPELIN/HB_NL_3000nmostvar/my.RNA.order_4clusters.RData")
bulk_order <- bulk_order[grep("T", names(bulk_order))] # remove non tumor livers
bulk_order <- sort(bulk_order, decreasing = T) # order as such: M, LP, H
bulk_exp <- bulk_exp[,names(bulk_order)]
H_tumors <- names(bulk_order[which(bulk_order %in% c(1,2))])
LP_tumors <- names(bulk_order[which(bulk_order %in% c(3))])

### Order tumors by scH/LP expression
bulk_exp_scH <- apply(bulk_exp[intersect(H_markers_sn, rownames(bulk_exp)), c(LP_tumors, H_tumors)], 2, mean); tumors_order_H <- names(sort(bulk_exp_scH, decreasing = F))
bulk_exp_scLP <- apply(bulk_exp[intersect(LP_markers_sn, rownames(bulk_exp)), c(LP_tumors, H_tumors)], 2, mean); tumors_order_LP <- names(sort(bulk_exp_scLP, decreasing = T))
tumors_order <- c(intersect(tumors_order_LP, LP_tumors), intersect(tumors_order_H, H_tumors))
  

####################################################
##### 1. MULTIOME TF EXPRESSION IN BULK RNASEQ #####
####################################################

exp. <- as.matrix(t(bulk_exp[TF_ordered,tumors_order])) # restrict bulk expression matrix to TF under study

for(j in 1:ncol(exp.)){
  exp.[,j] <- geco.changeRange( exp.[,j], newmin = 0, newmax = 1) # for each TF set values between 0 and 1
}

hmColors <- colorRampPalette(c("royalblue","white","indianred"))(256)# select heatmap gradient

pdf(file.path(resdir, "Heatmaps_TF/Heatmap_bulk_TF_targets_short_nb_targets_ordered_TF.pdf"))
par(mar=c(3,4,2,2))
image(exp.[,ncol(exp.):1], xaxt= "n", yaxt= "n", col = hmColors)
axis( 2, at=seq(0,1,length.out=ncol(exp.[,ncol(exp.):1]) ), labels=colnames(exp.[,ncol(exp.):1]), las= 2, cex.axis=0.6)
axis( 1, at=seq(0,1,length.out=nrow(exp.) ), labels=rownames(exp.), las= 2, cex.axis=0.6)
dev.off()

########################################################
##### 2. MULTIOME TARGET EXPRESSION IN BULK RNASEQ #####
########################################################

# create mean expression matrix for targets for each TF
exp_targets <- matrix(NA, nrow=length(TF_ordered), ncol=length(tumors_order)); rownames(exp_targets) <- TF_ordered; colnames(exp_targets) <- tumors_order
for(TF in TF_ordered){
  targets <- list_GRN[[TF]]
  exp_targets[TF, tumors_order] <- apply(bulk_exp[intersect(targets, rownames(bulk_exp)), tumors_order], 2, mean)
}

exp. <- as.matrix(t(exp_targets)) 

for(j in 1:ncol(exp.)){
  exp.[,j] <- geco.changeRange( exp.[,j], newmin = 0, newmax = 1) # for each TF set values between 0 and 1
}

hmColors <- colorRampPalette(c("royalblue","white","indianred"))(256)# select heatmap gradient

pdf(file.path(resdir, "Heatmaps_TF/Heatmap_bulk_TF_targets_short_nb_targets_ordered_targets.pdf"))
par(mar=c(3,4,2,2))
image(exp.[,ncol(exp.):1], xaxt= "n", yaxt= "n", col = hmColors)
axis( 2, at=seq(0,1,length.out=ncol(exp.[,ncol(exp.):1]) ), labels=colnames(exp.[,ncol(exp.):1]), las= 2, cex.axis=0.6)
axis( 1, at=seq(0,1,length.out=nrow(exp.) ), labels=rownames(exp.), las= 2, cex.axis=0.6)
dev.off()

