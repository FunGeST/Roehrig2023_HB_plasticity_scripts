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

keep_MT_genes=c("yes", "no")[1] # Keep MT genes for further analyses? 
regression=c("No_regression", "With_MT_regression", "With_cell_cycle_regression")[1] # Perform MT, cell cycle regression or not?

if(keep_MT_genes=="yes"){genes_MT_or_not <- "With_MT_genes"}else{genes_MT_or_not <- "Without_MT_genes"}

resdir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/No_regression/Merged"; if(!dir.exists(resdir)){dir.create(resdir)}
bulk_dir = "F:/Amelie-Datas/Pediatric tumors/2.Expression/RNAseq/Without new NL samples/Supervised" # directory containing H, LP, M differential bulk markers to project on single cell

nbcomp=40 # select number of components to use for UMAP


###################################
##### 0. DATA LOADING + MERGE ##### 
###################################

### Load and merge samples ###

sample_CHC2959N <- geco.load(paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/", regression, "/CHC2959N/", genes_MT_or_not, "/Seurat_object_analyses.RData"))
sample_CHC2959T <- geco.load(paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/", regression, "/CHC2959T/", genes_MT_or_not, "/Seurat_object_analyses.RData"))
sample_CHC2960T <- geco.load(paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/", regression, "/CHC2960T/", genes_MT_or_not, "/Seurat_object_analyses.RData"))
sample_CHC3133T <- geco.load(paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/", regression, "/CHC3133T/", genes_MT_or_not, "/Seurat_object_analyses.RData"))
sample_CHC3377N <- geco.load(paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/", regression, "/CHC3377N/", genes_MT_or_not, "/Seurat_object_analyses.RData"))
sample_CHC3377T <- geco.load(paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/", regression, "/CHC3377T/", genes_MT_or_not, "/Seurat_object_analyses.RData"))
sample_CHC3610T <- geco.load(paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/", regression, "/CHC3610T/", genes_MT_or_not, "/Seurat_object_analyses.RData"))
sample_CHC3662T <- geco.load(paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/", regression, "/CHC3662T/", genes_MT_or_not, "/Seurat_object_analyses.RData"))

# Merge samples on raw count matrix
sample_merged <- merge(sample_CHC2959N, list(sample_CHC2959T, sample_CHC2960T, sample_CHC3133T, sample_CHC3377N, sample_CHC3377T, sample_CHC3610T, sample_CHC3662T), add.cell.ids=c("CHC2959N", "CHC2959T","CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T"))
# merge merges RNA counts
#save(sample_merged, file=file.path(resdir, "Seurat_object_merged.RData"))

### Gene signatures ###

diff_H <- geco.load(paste0(bulk_dir, "/cluster_F_vs_HB/top_resullts_limma.RData"))
diff_LP <- geco.load(paste0(bulk_dir, "/cluster_E_vs_HB/top_resullts_limma.RData"))
diff_M <- geco.load(paste0(bulk_dir, "/cluster_M_vs_HB/top_resullts_limma.RData"))

LP_markers <- rownames(diff_LP)[which(diff_LP$logFC>3)] #133 genes
H_markers <- rownames(diff_H)[which(diff_H$logFC>3)] # 71 genes
M_markers <- rownames(diff_M)[which(diff_M$logFC>3)] # 304 genes
LP_markers <- LP_markers[which(LP_markers%in%rownames(sample_merged))]
H_markers <- H_markers[which(H_markers%in%rownames(sample_merged))]
M_markers <- M_markers[which(M_markers%in%rownames(sample_merged))]


############################
##### 1. PREPROCESSING #####
############################

sample <- sample_merged

### Normalization ###
sample <- SCTransform(sample, verbose = FALSE, variable.features.n=3000) # sctransform = par défaut 3000 nmostvar et pas 2000 comme lognormalize
# SCTransform always performs on RNA assay event if default assay is SCT (is written inside the funciton)

### Dimensionality reduction ###
sample <- RunPCA(sample, features = VariableFeatures(sample), npcs=100) # 100 components (sufficient to enable determination of the final number of components to use)

pdf(file.path(resdir, "PCA.pdf"))
DimPlot(sample, reduction="pca")
FeaturePlot(sample, features=c("nCount_RNA")) # on vérifie qu'il n'y a pas de biais par rapport aux UMI counts
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
DimPlot(sample, label=T, reduction="umap", group.by = "orig.ident") + NoLegend()


##### Check potential cluster-driving biases: cell cycle, UMI counts, MT%... #####
#---------------------------------------------------------------------------------

s.genes <- cc.genes$s.genes;s.genes <- s.genes[which(s.genes%in%rownames(sample))]
g2m.genes <- cc.genes$g2m.genes;g2m.genes <- g2m.genes[which(g2m.genes%in%rownames(sample))]
sample <- CellCycleScoring(sample, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf(file.path(resdir, "UMAP_potential_biases.pdf"))
DimPlot(sample, label=T, reduction="umap", group.by = "orig.ident") + NoLegend()
FeaturePlot(sample, features=c("percent.mt"))
FeaturePlot(sample, features=c("nCount_RNA"))
FeaturePlot(sample, features=c("ratio_UMI_genes"))
FeaturePlot(sample, features=c("log_ratio_genes_UMI"))
if(keep_MT_genes=="yes"){FeaturePlot(sample, features=c("MT-ND1", "MT-ND2", "MT-CO2", "MT-ATP8", "MT-ND3", "MT-CYB"), pt.size=0.1)} # pour vérifier si le fait de laisser les gènes MT induit un biais
DimPlot(sample, reduction = "umap", group.by= "Phase")
DimPlot(sample, reduction = "umap", group.by= "Phase", split.by = "Phase")
dev.off()

pdf(file.path(resdir, "UMAP_link_to_PCA.pdf"))
for(i in 1:10){print(FeaturePlot(sample, features=c(paste0("PC_", i))))}
dev.off()


##### Clustering #####
#---------------------

sample <- FindNeighbors(sample, dims=1:nbcomp)
sample <- FindClusters(sample, resolution =0.3) # change resolution to change nb of clusters determined

pdf(file.path(resdir, "UMAP_clusters.pdf"))
DimPlot(sample, label=T, reduction="umap") + NoLegend()
dev.off()

pdf(file.path(resdir, "Vlnplots_potential_biases_clusters.pdf"))
Idents(sample) <- "SCT_snn_res.0.3"
sample@meta.data %>%
  group_by(SCT_snn_res.0.3,orig.ident) %>%
  count() %>%
  group_by(SCT_snn_res.0.3) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=SCT_snn_res.0.3,y=percent, fill=orig.ident)) +
  geom_col() +
  ggtitle("Distribution of individual samples per cluster (res=0.3)")
sample@meta.data %>%
  group_by(SCT_snn_res.0.3,Phase) %>%
  count() %>%
  group_by(SCT_snn_res.0.3) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=SCT_snn_res.0.3,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster (res=0.3)")
VlnPlot(object = sample, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "ratio_UMI_genes"))
VlnPlot(object = sample, features = c("log_ratio_genes_UMI", "S.Score", "G2M.Score", "G2M.Score"))

dev.off()


##### VIsualization of normal/immune/tumor cells #####
#-----------------------------------------------------

sample$tumor_normal <- NA
for(s in c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")){
  
  infercnv_clusters <- geco.load(paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/4.Pseudo_CNV/1.Results/All_samples/UMAP_visualization/inferCNV_clusters_", s, ".RData"))
  names_infercnv_clusters <- paste0(s, "_", names(infercnv_clusters))
  
  if(s=="CHC2959N"){infercnv_clusters <- gsub("1", "normal/immune", gsub("2", "normal/immune", gsub("3", "normal/immune", infercnv_clusters)))}
  if(s=="CHC2959T"){infercnv_clusters <- gsub("1", "tumor", gsub("2", "tumor", gsub("3", "normal/immune", gsub("4", "tumor", infercnv_clusters))))}
  if(s=="CHC2960T"){infercnv_clusters <- gsub("1", "normal/immune", gsub("2", "tumor", gsub("3", "tumor", infercnv_clusters)))}
  if(s=="CHC3133T"){infercnv_clusters <- gsub("1", "tumor", gsub("2", "tumor", gsub("3", "tumor", gsub("4", "tumor", gsub("5", "normal/immune", infercnv_clusters)))))}
  if(s=="CHC3377N"){infercnv_clusters <- gsub("1", "normal/immune", gsub("2", "normal/immune", infercnv_clusters))}
  if(s=="CHC3377T"){infercnv_clusters <- gsub("1", "normal/immune", gsub("2", "tumor", gsub("3", "tumor", infercnv_clusters)))}
  if(s=="CHC3610T"){infercnv_clusters <- gsub("1", "normal/immune", gsub("2", "tumor", gsub("3", "tumor", gsub("4", "unknown", gsub("5", "tumor", infercnv_clusters)))))}
  if(s=="CHC3662T"){infercnv_clusters <- gsub("1", "normal/immune", gsub("2", "tumor", gsub("3", "tumor", gsub("4", "tumor", infercnv_clusters))))}
  
  names(infercnv_clusters) <- names_infercnv_clusters
  
  sample@meta.data[intersect(rownames(sample@meta.data), names(infercnv_clusters)),"tumor_normal"] <- infercnv_clusters[intersect(rownames(sample@meta.data), names(infercnv_clusters))]
}
sample@meta.data[grep("CHC2959N", rownames(sample@meta.data)),"tumor_normal"] <- "normal/immune"

pdf(file.path(resdir, "UMAP_tumor_vs_normal_immune.pdf"))
DimPlot(sample, group.by = "tumor_normal", label=T) + NoLegend()
DimPlot(sample, group.by = "orig.ident", label=T) + NoLegend()
dev.off()


##### VIsualization of the detailed individually annotated clusters #####
#------------------------------------------------------------------------

sample$defined_clusters <- NA
for(s in c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")){
  
  annotation_clusters <- geco.load(paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/No_regression/", s, "/With_MT_genes/cluster_refined_identification.RData"))
  annotation_clusters$cell <- paste0(s, "_", annotation_clusters$cell)
  
  sample@meta.data[annotation_clusters$cell,"defined_clusters"] <- annotation_clusters$cluster_definition
}

pdf(file.path(resdir, "UMAP_individual_cluster_definition.pdf"))
plot1 <- DimPlot(sample, group.by = "defined_clusters") + NoLegend()
LabelClusters(plot1, id = "defined_clusters", color = "black", size = 2, repel = T,  box.padding = 1)
DimPlot(sample, group.by = "orig.ident", label=T) + NoLegend()
dev.off()


##### H/LP/M signatures #####
#----------------------------

DefaultAssay(sample) <- "RNA"
sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor=10000) 

mean_LP <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[LP_markers,], 2, mean); sample$mean_LP <- mean_LP
mean_H <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[H_markers,], 2, mean); sample$mean_H <- mean_H
mean_M <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[M_markers,], 2, mean); sample$mean_M <- mean_M

### correlation between components and signatures 

pca_coord <- sample@reductions$pca@cell.embeddings
signature_exp <- sample@meta.data[,c("mean_H", "mean_LP", "mean_M")]
identical(rownames(pca_coord), rownames(signature_exp)) # must be TRUE

d <- as.matrix(cor(pca_coord, signature_exp, method = "pearson"))
save(d, file=file.path(resdir, "correlation_PCA_coordinates_H_LP_M_signatures.RData"))

pdf(file.path(resdir, "PCA_H_LP_M_signatures.pdf"))
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
p1 <- FeaturePlot(sample, reduction="pca", dims=1:2, features = c("mean_H", "mean_LP", "mean_M"), order=T, combine = FALSE)
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
p1 <- FeaturePlot(sample, reduction="pca", dims=c(1,3), features = c("mean_H", "mean_LP", "mean_M"), order=T, combine = FALSE)
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
p1 <- FeaturePlot(sample, reduction="pca", dims=c(1,5), features = c("mean_H", "mean_LP", "mean_M"), order=T, combine = FALSE)
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
p1 <- FeaturePlot(sample, reduction="pca", dims=c(3,5), features = c("mean_H", "mean_LP", "mean_M"), order=T, combine = FALSE)
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
dev.off()

pdf(file.path(resdir, "PCA_1_5_hepatocytes_location.pdf"))
sample$hepatocytes <- NA
sample@meta.data[grep("Hepatocyte", sample$defined_clusters),"hepatocytes"] <- "hepatocytes"
DimPlot(sample, reduction="pca", dims=c(1,5), group.by = "hepatocytes")
dev.off()



################
##### SAVE #####
################

save(sample, file=file.path(resdir, "Seurat_object_analyses.RData"))

