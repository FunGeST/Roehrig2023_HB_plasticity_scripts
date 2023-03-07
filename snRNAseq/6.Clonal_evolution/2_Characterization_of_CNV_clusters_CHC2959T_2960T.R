

library(dplyr)
library(ggplot2)
library(infercnv)
library(phytools)
library(dendextend)
library(Seurat)
library(ComplexHeatmap)


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

# We mix 2959T and 2960T because those are tumors from the same patient --> same ancester cells -->  they are handled in a different script

datadir <- paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/4.Pseudo_CNV/1.Results/All_samples_merged/") # we take the values for the heatmap containing all cells --> easier to compare
datadir_2959 <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/CHC2959T/Final_clusters/"
datadir_2960 <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/CHC2960T/Final_clusters/"

resdir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/CHC2959T_CHC2960T/Final_clusters/"

seurat_dir_tumor_cells_2959 <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/CHC2959T/"
seurat_dir_tumor_cells_2960 <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/CHC2960T/"

# load PCA coordinates for tumor cells on the merge UMAP
pca_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/PCA_components_study" # PCA on all tumor cells from all samples = not individual

# load somatic mutation info
mut_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/Somatic_mutations/Patient_2959/"

### To redo clustering on inferCNV results
preliminary_results <- c("yes", "no")[1] # check "yes" if the inferCNV algo has not finished yet and only the preliminary results are available (but it should not change the groups obtained, only the CNV values that should be clearer after filtering)
if(preliminary_results=="yes"){
  preliminary_suffix <- "preliminary."
}else{preliminary_suffix <- ""}

### choose clustering distance
clustering_distance = "euclidean"

### choose clustering method
clustering_method="ward.D"


###########################
##### 0. DATA LOADING #####
###########################

### load the inferCNV object (it contains the gene order, the kept genes and the inferCNV values used for the heatmaps)
object <- readRDS(paste0(datadir, preliminary_suffix, "infercnv_obj"))
object_gene_order <- object@gene_order # extract the gene coordinates
object_matrix <- object@expr.data # extract the values displayed on the inferCNV heatmap 

### Load final CNV clusters identified in script 4_Final_definition_CNV_clusters
clusters_2959 <- geco.load(file.path(datadir_2959, "final_CNV_clusters.RData"))
clusters_2960 <- geco.load(file.path(datadir_2960, "final_CNV_clusters.RData"))
# merge cluster tables
clusters <- rbind(clusters_2959, clusters_2960)

### Load seurat object for the sample
sample_tumor_cells_2959 <- geco.load(paste0(seurat_dir_tumor_cells_2959, "/Seurat_object_analyses.RData"))
sample_tumor_cells_2960 <- geco.load(paste0(seurat_dir_tumor_cells_2960, "/Seurat_object_analyses.RData"))
# merge the samples
sample_tumor_cells <- merge(sample_tumor_cells_2959, sample_tumor_cells_2960) # no need to add cell ID, already present
# re normalize the merge + PCA + UMAP
sample_tumor_cells <- SCTransform(sample_tumor_cells, verbose = FALSE, variable.features.n=3000) 
sample_tumor_cells <- RunPCA(sample_tumor_cells, features = VariableFeatures(sample_tumor_cells), npcs=100)
sample_tumor_cells <- RunUMAP(sample_tumor_cells, dims=1:40, verbose=F, umap.method = "uwot") # use 40 dimensions (will work fine in any case)

### Restrict inferCNV data to tumor cells from patient 2959
object_matrix <- object_matrix[,intersect(colnames(object_matrix), colnames(sample_tumor_cells))] # restrict to tumor cells

### Load PCA values (UMAP merged 6 tumor samples)
pca_coord_all <- geco.load(file.path(pca_dir, "PCA_coordinates.RData")) # load PCA coordinates for tumor cells on the merge UMAP
pca_coord <- pca_coord_all[grep("CHC2959T|CHC2960T", rownames(pca_coord_all)),] # restrict to tumor cells in the sample


########################################
##### 1. CLUSTERING ON BOTH TUMORS #####
########################################

### Create gene annotation (chromosome) ###
if(identical(rownames(object_gene_order), rownames(object_matrix))){
  row_ha =rowAnnotation(chromosomes = object_gene_order$chr)
} # check that the genes in matrix are ordered by position

### Create cell annotation (sample of origin) 
col_ha = columnAnnotation(sample = sample_tumor_cells@meta.data[colnames(as.matrix(object_matrix)),"orig.ident"],
                          col= list(sample = c("CHC2959T" = "#F8766D", "CHC2960T" = "#B79F00")))

### Print heatmap
# takes ~10 min

pdf(paste0(resdir, "/Complexheatmap_inferCNV_values_", clustering_distance, "_", clustering_method, ".pdf"))
ht <- draw(Heatmap(as.matrix(object_matrix), show_column_names = F,show_row_names=F, cluster_rows=F,  column_dend_height = unit(4, "cm"),
                   clustering_distance_columns = clustering_distance, clustering_method_columns = clustering_method,
                   column_title = paste0(clustering_distance, " distance, ", clustering_method, " method"),
                   left_annotation = row_ha, top_annotation = col_ha, col=circlize::colorRamp2(c(0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3), c("#00008E", "#4848AB", "#B6B6DD", "white", "#DDB6B6", "#AB4848", "#630000"))))
dev.off()

dendro <- column_dend(ht) # save the whole dendrogram
save(dendro, file=paste0(resdir, "/dendrogram_", clustering_distance, "_", clustering_method, ".RData"))

### Identify subclones 
clusters <- cutree(dendro, k=2) 
save(clusters, file=paste0(resdir, "/clusters_", clustering_distance, "_", clustering_method, ".RData"))


##############################
##### 2. ADD ANNOTATIONS #####
##############################

##### 1. Add CNV clusters information to seurat object #####
#-----------------------------------------------------------

sample_tumor_cells$CNV_clusters <- NA
if(length(setdiff(colnames(sample_tumor_cells), names(clusters)))==0){sample_tumor_cells@meta.data[names(clusters), "CNV_clusters"] <- clusters} # check that the cell names are the same in both objects
table(sample_tumor_cells$orig.ident, sample_tumor_cells$CNV_clusters)
sample_tumor_cells$CNV_clusters <- gsub(1, "CNV_CHC2959T", gsub(2, "CNV_CHC2960T", sample_tumor_cells$CNV_clusters))


##### 2. Add PC2 (H vs LP) values to seurat object #####
#-------------------------------------------------------

sample_tumor_cells$PC2_coord <- pca_coord[rownames(sample_tumor_cells@meta.data),"PC_2"]

##### 3.  Recompute G2M/S mean score #####
#-----------------------------------------

DefaultAssay(sample_tumor_cells) <- "RNA" # better to use RNA assay for mean on several genes
# get S and G2M genes (from Seurat)
s.genes <- cc.genes$s.genes;s.genes <- s.genes[which(s.genes%in%rownames(sample_tumor_cells))] # 42 genes
g2m.genes <- cc.genes$g2m.genes;g2m.genes <- g2m.genes[which(g2m.genes%in%rownames(sample_tumor_cells))] # 52 genes

sample_tumor_cells$G2M.Score <- sample_tumor_cells$S.Score <- sample_tumor_cells$Phase <- NULL # reinitialize cell cycle results to have the most recent ones
sample_tumor_cells <- CellCycleScoring(sample_tumor_cells, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(sample_tumor_cells) <- "SCT" # reset to SCT for visualizations etc


#################################
##### 2. UMAP VISUALIZATION #####
#################################

pdf(file.path(resdir, "UMAP_final_CNV_clusters_annotations.pdf"))
# sample
DimPlot(sample_tumor_cells, group.by = "orig.ident", label=T, pt.size = 1.5)+ ggtitle(paste0("CHC2959T+60T"))
# CNV clusters
DimPlot(sample_tumor_cells, group.by = "CNV_clusters", label=T, pt.size = 1.5, cols = c("#F6493C", "#B79F00"))+ ggtitle(paste0("CHC2959T+60T, final CNV clusters"))
DimPlot(sample_tumor_cells, group.by = "CNV_clusters", label=F, pt.size = 1.5, cols = c("#F6493C", "#B79F00"))+ ggtitle(paste0("CHC2959T+60T, final CNV clusters"))+NoLegend()
# tumor subtypes
DimPlot(sample_tumor_cells, group.by="H_LP_M_clean_groups_1", pt.size = 1.5, label=T, cols=c("#E8AAD6", "#BC53CD", "#510787", "#6F3323")) + ggtitle(paste0("CHC2959T+60T, tumor subtype"))
DimPlot(sample_tumor_cells, group.by="H_LP_M_clean_groups_1", pt.size = 1.5, label=F, cols=c("#E8AAD6", "#BC53CD", "#510787", "#6F3323")) + ggtitle(paste0("CHC2959T+60T, tumor subtype"))+NoLegend()
FeaturePlot(sample_tumor_cells, features="PC2_coord", pt.size = 1.5)+ scale_colour_gradientn(colours = c("#F3D5EA", "#E8AAD6", "#BC53CD", "#510787", "#1B002A"), values=scales::rescale(c(-50,-25,0,25,50)), limits=c(min(pca_coord_all[,"PC_2"]), max(pca_coord_all[,"PC_2"])))
FeaturePlot(sample_tumor_cells, features="PC2_coord", pt.size = 1.5)+ scale_colour_gradientn(colours = c("#F3D5EA", "#E8AAD6", "#BC53CD", "#510787", "#1B002A"), values=scales::rescale(c(-50,-25,0,25,50)), limits=c(min(pca_coord_all[,"PC_2"]), max(pca_coord_all[,"PC_2"]))) + NoLegend()
# proliferation cluster
sample_tumor_cells$cycling_cluster <- sample_tumor_cells$H_LP_M_seurat_clusters # we will retrieve the cells that belong to cycling cluster
sample_tumor_cells@meta.data[grep("prolif", sample_tumor_cells@meta.data$cycling_cluster),"cycling_cluster"] <- "cycling cluster"
sample_tumor_cells@meta.data[which(sample_tumor_cells$cycling_cluster != "cycling cluster"),"cycling_cluster"] <- NA # keep only cycling cluster in this column
DimPlot(sample_tumor_cells, group.by="cycling_cluster", label=T, pt.size = 1.5) + ggtitle(paste0("CHC2959T+60T, cycling cluster"))
DimPlot(sample_tumor_cells, group.by="cycling_cluster", label=F, pt.size = 1.5) + ggtitle(paste0("CHC2959T+60T, cycling cluster"))+NoLegend()
DimPlot(sample_tumor_cells, group.by="Phase", label=F, pt.size = 1.5) + ggtitle(paste0("CHC2959T+60T, phase"))
FeaturePlot(sample_tumor_cells, features="G2M.Score", pt.size = 1.5)+ scale_colour_gradientn(colours = c("grey", "darkorange"))
FeaturePlot(sample_tumor_cells, features="G2M.Score", pt.size = 1.5)+ scale_colour_gradientn(colours = c("grey", "darkorange")) + NoLegend()
dev.off()

dev.off()


#######################################
##### 3. STRIPCHARTS AND BARPLOTS #####
#######################################

# create table containing interesting info
df_stripchart <- factoall(sample_tumor_cells@meta.data[,c("CNV_clusters", "PC2_coord", "G2M.Score", "S.Score")])

pdf(file.path(resdir, "Stripcharts_final_CNV_clusters_annotations.pdf"))
### tumor subtype correlation (PC2 = H/LP gradient)
ggplot(df_stripchart, aes(x=CNV_clusters, y=PC2_coord, color=PC2_coord))+
  geom_jitter(width=length(unique(df_stripchart$CNV_clusters))/10)+
  theme_classic()+
  ggtitle("")+
  ylim(c(min(pca_coord_all[,"PC_2"]), max(pca_coord_all[,"PC_2"])))+ # to set all tumor samples on the same axis limits
  scale_colour_gradientn(colours = c("#E8AAD6", "#BC53CD", "#510787"), limits=c(min(pca_coord_all[,"PC_2"]), max(pca_coord_all[,"PC_2"])))+
  stat_summary(fun=median, geom="point", shape=18, color="red",size=5)

### perform Fisher test for cycling cluster
# we use fisher test and not chisquare because for some values we may have <5 (few observations) --> anyway if we did chisquare test we would have a warning "l'approximation du Chi-2 est peut-être incorrecte" meaning we should use fisher test instead (too small dataset)
sample_tumor_cells@meta.data[which(is.na(sample_tumor_cells$cycling_cluster)),"cycling_cluster"] <- "not cycling cluster" # we need to replace NA for Fisher test because otherwise we will have just 1 column
tt <- table(sample_tumor_cells$CNV_clusters, sample_tumor_cells$cycling_cluster); p_val <- format(fisher.test(tt)$p.value, digits=3)
ggplot(df_stripchart, aes(x=CNV_clusters, y=G2M.Score))+
  geom_jitter()+
  theme_classic()+
  ggtitle(paste0("Fisher test CNV clusters vs cycling cluster, pval= ", p_val))+
  geom_hline(yintercept = 0.1, linetype="dashed") # this threshold was defined using scores between cycling cluster and other cells + merged UMAP of tumor samples
ggplot(df_stripchart, aes(x=CNV_clusters, y=S.Score))+
  geom_jitter()+
  theme_classic()+
  ggtitle(paste0("Fisher test CNV clusters vs cycling cluster, pval= ", p_val))+
  geom_hline(yintercept = 0.1, linetype="dashed") # this threshold was defined using scores between cycling cluster and other cells + merged UMAP of tumor samples
### show distribution of cells with high or low cycling scores for each CNV cluster
df_stripchart$high_low_cycle_G2M_scores <- ifelse(df_stripchart$G2M.Score>0.1, "high (> 0.1)", "low (< 0.1)")
df_stripchart$high_low_cycle_S_scores <- ifelse(df_stripchart$S.Score>0.1, "high (> 0.1)", "low (< 0.1)")
# G2M scores
tt2 <- prop.table(table(df_stripchart$CNV_clusters, df_stripchart$high_low_cycle_G2M_scores), margin=1) # prop.table computes frequency from contingency table, margin=1 means the frequency is computed on rows = CNV clusters
par(mfrow = c(1, 2), mar=c(10,5,5,5))
barplot(t(tt2), col=c("darkorange", "grey90"), main="Cells with high/low G2M scores", las=2, ylab="% of cells in CNV cluster")
# S scores
tt3 <- prop.table(table(df_stripchart$CNV_clusters, df_stripchart$high_low_cycle_S_scores), margin=1) # prop.table computes frequency from contingency table, margin=1 means the frequency is computed on rows = CNV clusters
barplot(t(tt3), col=c("darkorange", "grey90"), main="Cells with high/low S scores", las=2, ylab="% of cells in CNV cluster")

dev.off()


################
##### SAVE #####
################

save(df_stripchart, file=file.path(resdir, "df_final_CNV_clusters_annotations.RData"))
