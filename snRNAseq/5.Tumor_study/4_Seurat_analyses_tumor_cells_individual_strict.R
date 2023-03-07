library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)

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


######################
##### PARAMETERS #####
######################

# select sample to study
sample_name="CHC3133T"

datadir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict"
resdir <- file.path(datadir, sample_name); if(!dir.exists(resdir)){dir.create(resdir)}

# Select number of components for UMAP etc
nbcomp=30

# Define colors for visualization (H, H+LP, LP, M subtypes)
cols1 = c("#E8AAD6", "#BC53CD", "#510787", "#6F3323")
cols2 = c("#E8AAD6", "#E088C3", "#BC53CD", "#8C3FC1", "#510787", "#6F3323")

if(sample_name=="CHC2959T"){cols1 <- cols1[c(1,2,3)]; cols2 <- cols2[c(2,3,4,5)]}
if(sample_name=="CHC2960T"){cols1 <- cols1[c(1,2,3)]; cols2 <- cols2[c(1,2,3,4,5)]}
if(sample_name=="CHC3133T"){cols1 <- cols1[c(1,2,3)]; cols2 <- cols2[c(2,3,4,5)]}
if(sample_name=="CHC3377T"){cols1 <- cols1[c(1,2,3)]; cols2 <- cols2[c(1,2,3,4)]}
if(sample_name=="CHC3662T"){cols1 <- cols1[c(1,2,3)]; cols2 <- cols2[c(1,2,3,4,5)]}
if(sample_name=="CHC3610T"){cols1 <- cols1[c(4)]; cols2 <- cols2[c(6)]}


###########################
##### 0. DATA LOADING #####
###########################

sample_all <- geco.load(file.path(datadir, "Seurat_object_analyses.RData")) # seurat object with all 14448 tumor cells
sample <- subset(sample_all, subset = orig.ident == sample_name) # keep only the cells from the sample of interest
DefaultAssay(sample) <- "SCT" # use SCT assay to redo normalization UMAP 

sample$mean_H <- sample$mean_LP <- sample$mean_M <- sample$SCT_snn_res.0.3 <- sample$SCT_snn_res.5 <- sample$clusters_Seurat_res_5 <- NULL # remove unnecessary metadata to avoid mistakes


############################
##### 1. PREPROCESSING #####
############################

### Normalization with SCT
sample <- SCTransform(sample, verbose = FALSE, variable.features.n=3000) # sctransform = par défaut 3000 nmostvar et pas 2000 comme lognormalize
# SCTransform works by default on RNA counts assay
DefaultAssay(sample) <- "RNA"
sample <- NormalizeData(sample) # normalize with classical pipeline for gene visualization (better to use "RNA" than "SCT" when averaging genes)
sample <- ScaleData(sample, features = rownames(sample)) # scale the data for potential heatmap projections
DefaultAssay(sample) <- "SCT" # return to SCT assay for next steps

### Add cell cycle (must be done after normalization because uses the data slot)
s.genes <- cc.genes$s.genes;s.genes <- s.genes[which(s.genes%in%rownames(sample))]
g2m.genes <- cc.genes$g2m.genes;g2m.genes <- g2m.genes[which(g2m.genes%in%rownames(sample))]
sample <- CellCycleScoring(sample, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

### Dimensionality reduction
sample <- RunPCA(sample, features = VariableFeatures(sample), npcs=100) # 100 components

pdf(file.path(resdir, "PCA.pdf"))
DimPlot(sample, reduction="pca")
FeaturePlot(sample, features=c("nCount_RNA"), reduction="pca") # on vérifie qu'il n'y a pas de biais par rapport aux UMI counts
ElbowPlot(sample, ndims=100)
ElbowPlot(sample, ndims=50)
DimHeatmap(sample, dims = 1:12, cells = 500, balanced=T) # l'option "cells" permet de garder les cellules les plus extrêmes = facilite le calcul
DimHeatmap(sample, dims = 13:24, cells = 500, balanced=T) 
DimHeatmap(sample, dims = 25:36, cells = 500, balanced=T) 
DimHeatmap(sample, dims = 37:40, cells = 500, balanced=T) 
dev.off()


#################################
##### 2. UMAP VISUALIZATION #####
#################################

sample <- RunUMAP(sample, dims=1:nbcomp, verbose=F, umap.method = "uwot")
DimPlot(sample, label=T, reduction="umap") + NoLegend()


##### Clustering #####
#---------------------

sample <- FindNeighbors(sample, dims=1:nbcomp)
sample <- FindClusters(sample, resolution =3) # change resolution to change nb of clusters determined

pdf(file.path(resdir, "UMAP_clusters.pdf"))
DimPlot(sample, label=T, reduction="umap") + NoLegend() + ggtitle(paste0("res=3"))
dev.off()


##### Display H LP M clean clusters from merged UMAP of tumor cells #####
#------------------------------------------------------------------------

pdf(file.path(resdir, "UMAP_H_LP_M_clean_groups_from_merged.pdf"))
DimPlot(sample, label=T, group.by = "H_LP_M_clean_groups_1") + NoLegend() + scale_color_manual(values=cols1)
DimPlot(sample, label=T, group.by = "H_LP_M_clean_groups_2") + NoLegend() + scale_color_manual(values=cols2)
dev.off()


#########################################################
##### 3. CROSSING H/LP/M GROUPS AND SEURAT CLUSTERS ##### OPTIONAL
#########################################################

# Goal = to have final clusters that are based on Seurat results (unsupervised unlike our signatures), match the H/LP/M groups based on PCA and keep info of particular clusters (prolif...)

table(sample$SCT_snn_res.3, sample$H_LP_M_clean_groups_1)
table(sample$SCT_snn_res.3, sample$defined_clusters)

sample$H_LP_M_seurat_clusters <- NA
if(sample_name=="CHC2959T"){
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(1,5,7,10,19)),"H_LP_M_seurat_clusters"] <- "H+LP"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(0,3,4,6,8,9,11,12,13,16,17)),"H_LP_M_seurat_clusters"] <- "LP"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(2,15)),"H_LP_M_seurat_clusters"] <- "LP, prolif"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3==18),"H_LP_M_seurat_clusters"] <- "LP, low counts"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3==14),"H_LP_M_seurat_clusters"] <- "LP, high MT%"
  cols3 <- c("#BC53CD", rep("#510787", times=4))
}
if(sample_name=="CHC2960T"){
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(1,2,3,4,5,6,7,10,11,13,16,21)),"H_LP_M_seurat_clusters"] <- "H"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(0,8,9,15,18,19,23)),"H_LP_M_seurat_clusters"] <- "H+LP"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(14,17)),"H_LP_M_seurat_clusters"] <- "LP"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(12,20,22)),"H_LP_M_seurat_clusters"] <- "H+LP, prolif"
  cols3 <- c("#E8AAD6", "#BC53CD", "#BC53CD", "#510787")
}
if(sample_name=="CHC3133T"){
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(1,2,3,5,6,8,10,12,15,17,19,20,21,25,26,27)),"H_LP_M_seurat_clusters"] <- "H+LP"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(0,4,7,9,11,13,14,16,18,22)),"H_LP_M_seurat_clusters"] <- "LP"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(0,9)),"H_LP_M_seurat_clusters"] <- "LP, prolif"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(24)),"H_LP_M_seurat_clusters"] <- "LP, M?"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(23)),"H_LP_M_seurat_clusters"] <- "H+LP, high ALB"
  cols3 <- c("#BC53CD","#BC53CD", "#510787", "#510787", "#510787")
}
if(sample_name=="CHC3377T"){
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(5,8,12,15,28)),"H_LP_M_seurat_clusters"] <- "H"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(0,1,3,4,6,7,9,10,11,13,14,16,17,18,19,20,21,22,23,24,25,26,27,29,30)),"H_LP_M_seurat_clusters"] <- "H+LP"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(2)),"H_LP_M_seurat_clusters"] <- "H+LP, high MT%"
  cols3 <- c("#E8AAD6", "#BC53CD", "#BC53CD")
}
if(sample_name=="CHC3662T"){
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(0,1,4,6,8,9,17,18,20,21)),"H_LP_M_seurat_clusters"] <- "H"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(2,3,5,7,10,11,12,13,14,15,16,19)),"H_LP_M_seurat_clusters"] <- "H+LP"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(23)),"H_LP_M_seurat_clusters"] <- "LP"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(22)),"H_LP_M_seurat_clusters"] <- "H, high ALB"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(19)),"H_LP_M_seurat_clusters"] <- "H+LP, prolif"
  cols3 <- c("#E8AAD6", "#E8AAD6", "#BC53CD", "#BC53CD", "#BC53CD")
}
if(sample_name=="CHC3610T"){
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(0,1,2,3,4,5,6,7,9,10,11,12,13,14,17,22,24)),"H_LP_M_seurat_clusters"] <- "M"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(16,18)),"H_LP_M_seurat_clusters"] <- "M, GLI3"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(15)),"H_LP_M_seurat_clusters"] <- "M, prolif"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(8)),"H_LP_M_seurat_clusters"] <- "M, low counts"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(19)),"H_LP_M_seurat_clusters"] <- "M, low counts, high MT%"
  sample@meta.data[which(sample@meta.data$SCT_snn_res.3%in%c(20,21,23)),"H_LP_M_seurat_clusters"] <- "M, clustered with immune cells"
  cols3 <- c("#6F3323", "grey", rep("#6F3323", times=4))
}

pdf(file.path(resdir, "UMAP_H_LP_M_seurat_clusters.pdf"))
DimPlot(sample, group.by = "H_LP_M_seurat_clusters")
DimPlot(sample, group.by = "H_LP_M_seurat_clusters") + NoLegend() + scale_color_manual(values=cols3)
DimPlot(sample, label=T, group.by = "H_LP_M_seurat_clusters") + NoLegend() + scale_color_manual(values=cols3)
dev.off()


###################
##### 3. SAVE #####
###################

sample@meta.data <- sample@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.RPS", "ratio_UMI_genes", "log_ratio_genes_UMI", "S.Score", "G2M.Score", "Phase", "doublet_clusters_potential", "H_LP_M_clean_groups_1", "H_LP_M_clean_groups_2", "SCT_snn_res.3", "H_LP_M_seurat_clusters")]
save(sample, file=file.path(resdir, "Seurat_object_analyses.RData"))

