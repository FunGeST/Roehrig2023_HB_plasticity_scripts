### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(Seurat)
library(reticulate)

geco.load <- function (filename) # File to be loaded
{
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error - function geco.load : file ", filename, 
            " not found"))
  NULL
}

######################
##### PARAMETERS #####
######################

datadir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict"
resdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/3.Denoising/Imputation_in_house_tumor_cells"

nbcomp=40 # nb of components to use for clustering


###########################
##### 0. DATA LOADING #####
###########################

sample <- geco.load(file.path(datadir, "Seurat_object_analyses.RData"))
exp_mat <- as.matrix(GetAssayData(sample, assay="RNA", slot="data")) # retrieve matrix of normalized expression data
pca_values <- sample@reductions[["pca"]]@cell.embeddings # retrieve PCA cell coordinates


##################################################
##### 1. CLUSTER-BASED MEAN EXPRESSION TABLE #####
##################################################

##### 1. Identify cell clusters #####
#------------------------------------
# Separate the data in numerous small clusters = like mini pseudo-bulk samples
DefaultAssay(sample) <- "SCT"
sample <- FindNeighbors(sample, dims = 1:nbcomp)
sample <- FindClusters(sample, resolution = 20) # identify clusters on which aggregate the expression data
clusters_order <- paste0("C", levels(sample$seurat_clusters))
sample$seurat_clusters <- paste0("C", as.character(sample$seurat_clusters)) # clusters are automatically set in the "seurat_clusters" slot
cells_in_cluster <- sample$seurat_clusters; names(cells_in_cluster) <- rownames(sample@meta.data) # identify cells in each cluster
save(cells_in_cluster, file=file.path(resdir, "cells_in_clusters.RData"))

DimPlot(sample, group.by="seurat_clusters", label=T)+NoLegend() # project on UMAP

##### 2. Create expression table #####
#-------------------------------------
exp_mat_clusters <- as.data.frame(matrix(NA, nrow=nrow(exp_mat), ncol=length(clusters_order)))
rownames(exp_mat_clusters) <- rownames(exp_mat); colnames(exp_mat_clusters) <- clusters_order

for(cluster in colnames(exp_mat_clusters)){ # for each cell cluster compute the average expression for each gene
  cells_in_cluster <- rownames(sample@meta.data[which(sample$seurat_clusters == cluster),])
  exp_mat_clusters[,cluster] <- apply(exp_mat[,cells_in_cluster], 1, mean)
}
save(exp_mat_clusters, file=file.path(resdir, "exp_mat_clusters.RData"))

##### 3. Add cluster annotation #####
#------------------------------------
cluster_annot <- data.frame(cluster=clusters_order, sample=NA, mean_PC1=NA, mean_PC2=NA, stringsAsFactors = F)
rownames(cluster_annot) <- cluster_annot$cluster

for(cluster in cluster_annot$cluster){ 
  # Add sample of origin
  cluster_annot[cluster,"sample"] <- paste(unique(sample@meta.data[which(sample$seurat_clusters==cluster),"orig.ident"]), collapse=", ")
  # Add mean PC
  cluster_annot[cluster,"mean_PC1"] <- mean(pca_values[rownames(sample@meta.data[which(sample$seurat_clusters==cluster),]),"PC_1"])
  cluster_annot[cluster,"mean_PC2"] <- mean(pca_values[rownames(sample@meta.data[which(sample$seurat_clusters==cluster),]),"PC_2"])
}

save(cluster_annot, file=file.path(resdir, "cluster_annot.RData"))


##################################################################
##### 2. PROJECT CLUSTER MEAN EXPRESSION ON INDIVIDUAL CELLS #####
##################################################################

exp_mat_clusters_per_cell <- exp_mat
for(cluster in colnames(exp_mat_clusters)){
  cells_in_cluster <- rownames(sample@meta.data[which(sample$seurat_clusters == cluster),])
  exp_mat_clusters_per_cell[,cells_in_cluster] <- exp_mat_clusters[,cluster]
}

save(exp_mat_clusters_per_cell, file=file.path(resdir, "exp_mat_clusters_per_cell.RData"))

## Test on Featureplots
sample$test <- exp_mat_clusters_per_cell["LHX1",]
FeaturePlot(sample, features="test")
FeaturePlot(sample, features="LHX1")


