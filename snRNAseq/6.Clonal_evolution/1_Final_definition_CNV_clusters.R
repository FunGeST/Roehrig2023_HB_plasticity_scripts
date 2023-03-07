
library(dplyr)
library(ggplot2)
library(infercnv)
library(phytools)
library(dendextend)
library(Seurat)
library(ComplexHeatmap)
library(lsa)


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

geco.changeRange <- function (v, newmin = 1, newmax = 10) 
{
  oldmin <- min(v, na.rm = TRUE)
  oldmax <- max(v, na.rm = TRUE)
  newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin))
}

distCosine <- function(m)
{
  nsamp <- nrow(m)
  res <- matrix(NA,nrow=nsamp,ncol=nsamp)
  rownames(res) <- colnames(res) <- rownames(m)
  for(i in 1:nsamp)
  {
    for(j in 1:nsamp)
    {
      res[i,j] <- cosine(m[i,],m[j,])
    }
  }
  as.dist(1-res)
}




######################
##### PARAMETERS #####
######################

# Here we use the inferCNV values computed on the merged dataset of all 8 samples  = directory "All_samples_merged"
# for individual study we will take these values and project them on the individuel seurat objects
name_sample = c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3377T", "CHC3610T", "CHC3662T")[6] 

datadir <- paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/4.Pseudo_CNV/1.Results/All_samples_merged/") # we take the values for the heatmap containing all cells --> easier to compare
resdir <- paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/", name_sample, "/Final_clusters"); if(!dir.exists(resdir)){dir.create(resdir)}
seurat_dir <- paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/", name_sample)

preliminary_results <- c("yes", "no")[1] # check "yes" if the inferCNV algo has not finished yet and only the preliminary results are available (but it should not change the groups obtained, only the CNV values that should be clearer after filtering)
if(preliminary_results=="yes"){
  preliminary_suffix <- "preliminary."
}else{preliminary_suffix <- ""}


### choose clustering distance
clustering_distance = "euclidean"

### choose clustering method
clustering_method="ward.D"

# euclidean distance and ward.D method were chosen because they gave good results 

# how many clusters?
if(name_sample=="CHC2959T"){j=3}else if(name_sample=="CHC2960T"){j=9}else if(name_sample=="CHC3377T"){j=7}else if(name_sample=="CHC3133T"){j=4}else if(name_sample=="CHC3610T"){j=4}else if(name_sample=="CHC3662T"){j=7}


###########################
##### 0. DATA LOADING #####
###########################

### load the inferCNV object (it contains the gene order, the kept genes and the inferCNV values used for the heatmaps)
object <- readRDS(paste0(datadir, preliminary_suffix, "infercnv_obj"))
object_gene_order <- object@gene_order # extract the gene coordinates
object_matrix <- object@expr.data # extract the values displayed on the inferCNV heatmap 

### load seurat object (tumor cells)
seurat_sample <- geco.load(file.path(seurat_dir, "Seurat_object_analyses.RData"))

object_matrix <- object_matrix[,intersect(colnames(object_matrix), colnames(seurat_sample))] # restrict to tumor cells


##########################################################
##### 1. CLUSTERING ON INFERCNV VALUES (PRELIMINARY) #####
##########################################################

### Create gene annotation (chromosome) #####

if(identical(rownames(object_gene_order), rownames(object_matrix))){
  row_ha =rowAnnotation(chromosomes = object_gene_order$chr)
  } # check that the genes in matrix are ordered by position


### Print heatmap
# takes ~10 min

pdf(paste0(resdir, "/Complexheatmap_inferCNV_values_", clustering_distance, "_", clustering_method, ".pdf"))
ht <- draw(Heatmap(as.matrix(object_matrix), show_column_names = F,show_row_names=F, cluster_rows=F,  column_dend_height = unit(4, "cm"),
                   clustering_distance_columns = clustering_distance, clustering_method_columns = clustering_method,
                   column_title = paste0(clustering_distance, " distance, ", clustering_method, " method"),
                   left_annotation = row_ha, col=circlize::colorRamp2(c(0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3), c("#00008E", "#4848AB", "#B6B6DD", "white", "#DDB6B6", "#AB4848", "#630000"))))
dev.off()

dendro <- column_dend(ht) # save the whole dendrogram
save(dendro, file=paste0(resdir, "/dendrogram_", clustering_distance, "_", clustering_method, ".RData"))

### split the dendrogram by expected clusters

dendro_clusters = color_branches(dendro, k = j) 

pdf(paste0(resdir, "/dendrogram_clusters_k_", j, "_", clustering_distance, "_", clustering_method, ".pdf"))
plot(dendro_clusters)
dev.off()


##############################################
##### 2. ANNOTATE THE FINAL CNV CLUSTERS #####
##############################################

# Determiner order of cells in clustering
cell_order <- colnames(as.matrix(object_matrix))[ht@ht_list[[grep("matrix", names(ht@ht_list), value=T)]]@column_order]
final_cnv_clusters <- factoall(data.frame(cells=cell_order, inferCNV_cluster=NA, final_CNV_cluster=NA))

if(name_sample=="CHC2959T"){
  clusters <- cutree(dendro, k=j)[cell_order] # we order clusters in the same order as heatmap (otherwise it is from 1 to n)
  final_cnv_clusters$inferCNV_cluster <- clusters # retrieve inferCNV clusters for each cell
  # defined final clusters based on inferCNV clusters (complexheatmap, ward.D euclidean)
  final_cnv_clusters[which(final_cnv_clusters$inferCNV_cluster %in% c(1,2,3)),"final_CNV_cluster"] <- "+1q, 2, 7, 8, 16, 20, 21"
  
  save(final_cnv_clusters, file=file.path(resdir, "final_CNV_clusters.RData"))
}else if(name_sample=="CHC2960T"){
  clusters <- cutree(dendro, k=j)[cell_order] # we order clusters in the same order as heatmap (otherwise it is from 1 to n)
  final_cnv_clusters$inferCNV_cluster <- clusters # retrieve inferCNV clusters for each cell
  # defined final clusters based on inferCNV clusters (complexheatmap, ward.D euclidean)
  final_cnv_clusters[which(final_cnv_clusters$inferCNV_cluster %in% c(8)),"final_CNV_cluster"] <- "+1q, 2q, 5, 6, 8, 12, 20, 2p"
  final_cnv_clusters[which(final_cnv_clusters$inferCNV_cluster %in% c(1,2,3,4,5,6,7,9)),"final_CNV_cluster"] <- "+1q, 2q, 5, 6, 8, 12, 20"
  
  save(final_cnv_clusters, file=file.path(resdir, "final_CNV_clusters.RData"))
}else if(name_sample=="CHC3133T"){
  clusters <- cutree(dendro, k=j)[cell_order] # we order clusters in the same order as heatmap (otherwise it is from 1 to n)
  final_cnv_clusters$inferCNV_cluster <- clusters # retrieve inferCNV clusters for each cell
  # defined final clusters based on inferCNV clusters (complexheatmap, ward.D euclidean)
  final_cnv_clusters[which(final_cnv_clusters$inferCNV_cluster %in% c(1,2,3)),"final_CNV_cluster"] <- "+1q, 2, 3, 4, 5, 6, 8, 12, 14, 16, 17, 19, 20"
  final_cnv_clusters[which(final_cnv_clusters$inferCNV_cluster %in% c(4)),"final_CNV_cluster"] <- "+1q, 2, 4, 5, 6, 8, 12, 14, 16, 17, 19, 20"
  
  save(final_cnv_clusters, file=file.path(resdir, "final_CNV_clusters.RData"))
}else if(name_sample=="CHC3377T"){
  clusters <- cutree(dendro, k=j)[cell_order] # we order clusters in the same order as heatmap (otherwise it is from 1 to n)
  final_cnv_clusters$inferCNV_cluster <- clusters # retrieve inferCNV clusters for each cell
  # defined final clusters based on inferCNV clusters (complexheatmap, ward.D euclidean)
  final_cnv_clusters[which(final_cnv_clusters$inferCNV_cluster %in% c(1)),"final_CNV_cluster"] <- "+2q, 6, 12, 17, 20, X, 1q"
  final_cnv_clusters[which(final_cnv_clusters$inferCNV_cluster %in% c(7)),"final_CNV_cluster"] <- "+2q, 6, 12, 17, 20, X, 8"
  final_cnv_clusters[which(final_cnv_clusters$inferCNV_cluster %in% c(2,3,4,5,6)),"final_CNV_cluster"] <- "+2q, 6, 12, 17, 20, X"
  
  save(final_cnv_clusters, file=file.path(resdir, "final_CNV_clusters.RData"))
}else if(name_sample=="CHC3610T"){
  clusters <- cutree(dendro, k=j)[cell_order] # we order clusters in the same order as heatmap (otherwise it is from 1 to n)
  final_cnv_clusters$inferCNV_cluster <- clusters # retrieve inferCNV clusters for each cell
  # defined final clusters based on inferCNV clusters (complexheatmap, ward.D euclidean)
  final_cnv_clusters[which(final_cnv_clusters$inferCNV_cluster %in% c(1,2,3,4)),"final_CNV_cluster"] <- "+2, 7, 8, 19, 20"
  
  save(final_cnv_clusters, file=file.path(resdir, "final_CNV_clusters.RData"))
}else if(name_sample=="CHC3662T"){
  clusters <- cutree(dendro, k=j)[cell_order] # we order clusters in the same order as heatmap (otherwise it is from 1 to n)
  final_cnv_clusters$inferCNV_cluster <- clusters # retrieve inferCNV clusters for each cell
  # defined final clusters based on inferCNV clusters (complexheatmap, ward.D euclidean)
  final_cnv_clusters[which(final_cnv_clusters$inferCNV_cluster %in% c(5,2,6)),"final_CNV_cluster"] <- "+2q, 17, 20"
  final_cnv_clusters[which(final_cnv_clusters$inferCNV_cluster %in% c(4,1,7)),"final_CNV_cluster"] <- "+2q, 17, 20, 1q"
  final_cnv_clusters[which(final_cnv_clusters$inferCNV_cluster %in% c(3)),"final_CNV_cluster"] <- "2q, 17, 1q, 2p, 8, 12"

  save(final_cnv_clusters, file=file.path(resdir, "final_CNV_clusters.RData"))
}




