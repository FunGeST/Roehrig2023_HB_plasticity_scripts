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

datadir_rna = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/"
resdir_visium = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/Visium"


###########################
##### 0. DATA LOADING #####
###########################

### snRNA-seq data 
sample_seurat <- geco.load(file.path(datadir_rna, "Seurat_object_analyses.RData")) # load object with all tumor cells from all samples to get PC2 values
pc1_values <- sample_seurat@reductions$pca@feature.loadings[,"PC_1"] # load pc1 values
pc2_values <- sample_seurat@reductions$pca@feature.loadings[,"PC_2"] # load pc2 values
#Define H and LP markers from PC2 by taking top 120 genes of PC2 extremities
LP_markers_sn <- names(sort(pc2_values, decreasing = T)[1:120])
H_markers_sn <- names(sort(pc2_values, decreasing = F)[1:120]) 
M_markers_sn <- names(sort(pc1_values, decreasing = T)[1:120]) 
df_markers <- data.frame(scH_markers=H_markers_sn, scLP_markers=LP_markers_sn, scM_markers=M_markers_sn, stringsAsFactors = F)

### Visium data
sample_visium <- geco.load(file.path(resdir_visium, "0.Input_data/seurat_SCT_sample_3131_T.Rdata"))
DefaultAssay(sample_visium) <- "SCT"

exp_visium <- sample_visium@assays$SCT@data # extract expression matrix

# Restrict H/LP markers to common genes with Visium (have approximately 100 genes)
LP_markers_sn <- intersect(LP_markers_sn, rownames(sample_visium)) #102 genes
H_markers_sn <- intersect(H_markers_sn, rownames(sample_visium)) #91 genes
M_markers_sn <- intersect(M_markers_sn, rownames(sample_visium)) #90 genes
df_markers$in_visium_H <- ifelse(df_markers$scH_markers %in% H_markers_sn, "yes", "no") # are the genes found in visium or not?
df_markers$in_visium_LP <- ifelse(df_markers$scLP_markers %in% LP_markers_sn, "yes", "no") # are the genes found in visium or not?
df_markers$in_visium_M <- ifelse(df_markers$scM_markers %in% M_markers_sn, "yes", "no") # are the genes found in visium or not?
write.table(df_markers, file.path(resdir_visium, "sc_markers_in_visium.csv"), col.names=T, row.names=F, sep=";")


###########################################################
##### 1. COMPUTE H/LP MEAN EXPRESSION PER SPATIAL DOT #####
###########################################################

##### 1. Create H and LP average expression #####
#------------------------------------------------
sample_visium$mean_H_sn <- apply(exp_visium[H_markers_sn,], 2, mean)
sample_visium$mean_LP_sn <- apply(exp_visium[LP_markers_sn,], 2, mean)

##### 2. Project on visium #####
#-------------------------------
pdf(file.path(resdir_visium, "H_LP_snRNAseq_makers_on_visium_map.pdf"))
SpatialFeaturePlot(sample_visium, features="mean_H_sn", pt.size.factor = 1.6)
SpatialFeaturePlot(sample_visium, features="mean_LP_sn", pt.size.factor = 1.6)
dev.off()

##### 3. 2D plot of H vs LP markers in visium #####
#--------------------------------------------------

ggplot(sample_visium@meta.data, aes(x=mean_H_sn, y=mean_LP_sn))+
  geom_point()+
  theme_classic()

