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

sample_name="CHC3662T"
keep_MT_genes="yes"


datadir <-paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/1.Preprocess/Filtered_genes_gtf_files_5/", sample_name)
resdir <- paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/No_regression/", sample_name)
if(!dir.exists(resdir)){dir.create(resdir)}

if(keep_MT_genes=="yes"){
  resdir <- paste0(resdir, "/With_MT_genes")
}else{
  resdir <- paste0(resdir, "/Without_MT_genes")
}

if(!dir.exists(resdir)){dir.create(resdir)}


###########################
##### 0. DATA LOADING #####
###########################


if(keep_MT_genes=="yes"){
  sample <- geco.load(file.path(datadir, "Seurat_object_QC_with_MT_genes.RData"))
}else{
  sample <- geco.load(file.path(datadir, "Seurat_object_QC_without_MT_genes.RData"))
}


############################
##### 1. PREPROCESSING #####
############################

sample <- SCTransform(sample, verbose = FALSE, variable.features.n=3000) # sctransform = par défaut 3000 nmostvar et pas 2000 comme lognormalize
DefaultAssay(sample) <- "RNA"
sample <- NormalizeData(sample) # pour avoir aussi le slot data normalisé dans sample au cas où
sample <- ScaleData(sample, features = rownames(sample)) 
DefaultAssay(sample) <- "SCT"

sample <- RunPCA(sample, features = VariableFeatures(sample), npcs=100) # 100 components

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

nbcomp=30


#################################
##### 2. UMAP VISUALIZATION #####
#################################

sample <- RunUMAP(sample, dims=1:nbcomp, verbose=F, umap.method = "uwot")
DimPlot(sample, label=T, reduction="umap") + NoLegend()


##### Check potential cluster-driving biases: cell cycle, UMI counts, MT%... #####
#---------------------------------------------------------------------------------

s.genes <- cc.genes$s.genes;s.genes <- s.genes[which(s.genes%in%rownames(sample))]
g2m.genes <- cc.genes$g2m.genes;g2m.genes <- g2m.genes[which(g2m.genes%in%rownames(sample))]
sample <- CellCycleScoring(sample, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

sample$mean_G2M <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[g2m.genes,], 2, mean)
sample$mean_S <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[s.genes,], 2, mean)

# just check ribosomal genes in case
sample[["percent.RPS"]] <- PercentageFeatureSet(sample, pattern="^RPS") 
sample[["percent.RPL"]] <- PercentageFeatureSet(sample, pattern="^RPL")  

pdf(file.path(resdir, "UMAP_potential_biases.pdf"))
# typical sources of variation
FeaturePlot(sample, features=c("percent.mt"))
FeaturePlot(sample, features=c("percent.RPS"))
FeaturePlot(sample, features=c("percent.RPL"))
FeaturePlot(sample, features=c("nCount_RNA"))
FeaturePlot(sample, features=c("nFeature_RNA"))
FeaturePlot(sample, features=c("ratio_UMI_genes"))
FeaturePlot(sample, features=c("log_ratio_genes_UMI"))
if(keep_MT_genes=="yes"){FeaturePlot(sample, features=c("MT-ND1", "MT-ND2", "MT-CO2", "MT-ATP8", "MT-ND3", "MT-CYB"), pt.size=0.1)} # pour vérifier si le fait de laisser les gènes MT induit un biais
# cell cycle
FeaturePlot(sample, features=c("mean_G2M"))
FeaturePlot(sample, features=c("mean_S"))
DimPlot(sample, reduction = "umap", group.by= "Phase")
DimPlot(sample, reduction = "umap", group.by= "Phase", split.by = "Phase")

sample_cc <- sample; DefaultAssay(sample_cc) <- "RNA"
sample_cc <- RunPCA(sample_cc, features = c(s.genes, g2m.genes), npcs=100) # 100 components
sample_cc$Phase <- factor(sample_cc$Phase, levels=c("G1", "G2M", "S"))
DimPlot(sample_cc, group.by = "Phase")
FeaturePlot(sample_cc, features="mean_G2M")
FeaturePlot(sample_cc, features="mean_S")

dev.off()

# First identification using PCA

pca_genes <- sample@reductions$pca@feature.loadings
df_pca_genes <- as.data.frame(matrix(NA, nrow=20, ncol=20))
a <- c()
   for(i in 1:10){
    a <- c(a, paste0("PC_", i, "_pos"), paste0("PC_", i, "_neg"))
   }
colnames(df_pca_genes) <- a
  
pdf(file.path(resdir, "UMAP_link_to_PCA.pdf"))
for(i in 0:9){
  print(FeaturePlot(sample, features=c(paste0("PC_", i+1))))
  genes_order <- names(sort(pca_genes[,i+1], decreasing = T)[1:20])
  df_pca_genes[,i*2+1] <- genes_order
  sample$mean_PC <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[genes_order,], 2, mean)
  print(FeaturePlot(sample, features="mean_PC")+ scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+ ggtitle(label=paste0("PC_", i+1, " pos (20 top genes)")))
  genes_order <- names(sort(pca_genes[,i+1])[1:20])
  df_pca_genes[,i*2+2] <- genes_order
  sample$mean_PC <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[genes_order,], 2, mean)
  print(FeaturePlot(sample, features="mean_PC")+ scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+ ggtitle(label=paste0("PC_", i+1, " neg (20 top genes)")))
}
dev.off()
sample$mean_PC <- NULL

save(df_pca_genes, file=file.path(resdir, "PCA_genes.RData"))


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
  group_by(SCT_snn_res.0.3,Phase) %>%
  count() %>%
  group_by(SCT_snn_res.0.3) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=SCT_snn_res.0.3,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster (res=0.3)")
VlnPlot(object = sample, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.RPS", "percent.RPL", "ratio_UMI_genes"))
VlnPlot(object = sample, features = c("log_ratio_genes_UMI", "S.Score", "G2M.Score", "mean_G2M", "mean_S"))
dev.off()

save(sample, file=file.path(resdir, "Seurat_object_analyses.RData"))

