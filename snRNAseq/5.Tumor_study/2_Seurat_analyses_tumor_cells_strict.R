library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)
library(corrplot)

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

datadir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune"
resdir <- file.path(datadir, "Tumor/2.Tumor_selection_strict"); if(!dir.exists(resdir)){dir.create(resdir)}

bulk_dir = "F:/Amelie-Datas/Pediatric tumors/2.Expression/RNAseq/Without new NL samples/Supervised"

nbcomp=40 # nb of components for UMAP

###########################
##### 0. DATA LOADING ##### 
###########################

sample_all <- geco.load(file.path(datadir, "Seurat_object_analyses_all_cells.RData"))
sample <- subset(sample_all, subset = selection_proposition_strict == "tumor (clean)") # 14448 cells


##### Gene signatures #####
#--------------------------

diff_H <- geco.load(paste0(bulk_dir, "/cluster_F_vs_HB/top_resullts_limma.RData"))
diff_LP <- geco.load(paste0(bulk_dir, "/cluster_E_vs_HB/top_resullts_limma.RData"))
diff_M <- geco.load(paste0(bulk_dir, "/cluster_M_vs_HB/top_resullts_limma.RData"))

LP_markers <- rownames(diff_LP)[which(diff_LP$logFC>3)] #133 genes
H_markers <- rownames(diff_H)[which(diff_H$logFC>3)] # 71 genes
M_markers <- rownames(diff_M)[which(diff_M$logFC>3)] # 304 genes
LP_markers <- LP_markers[which(LP_markers%in%rownames(sample))]
H_markers <- H_markers[which(H_markers%in%rownames(sample))]
M_markers <- M_markers[which(M_markers%in%rownames(sample))]
  

############################
##### 1. PREPROCESSING #####
############################

# Normalization with SCTransform
sample <- SCTransform(sample, verbose = FALSE, variable.features.n=3000) # sctransform = par défaut 3000 nmostvar et pas 2000 comme lognormalize

# Dimensionality reduction with PCA
sample <- RunPCA(sample, features = VariableFeatures(sample), npcs=100) # 100 components

###########################
##### 2. PCA ANALYSIS #####
###########################

# Extract PCA genes
pca_genes <- sample@reductions$pca@feature.loadings
save(pca_genes, file=file.path(resdir, "PCA_genes.RData"))

# Show PCA variance
pdf(file.path(resdir, "PCA_variance.pdf"))
ElbowPlot(sample, ndims=100)
dev.off()

# Project PCA main components together
pdf(file.path(resdir, "PCA_main_components.pdf"))
for(i in 1:4){
  for(j in 1:4){
    if(i!=j){
      print(DimPlot(sample, reduction="pca", group.by = "orig.ident", dims=c(i,j)))
    }
  }
}
dev.off()

# Project top genes per PCA component
for(i in 1:3){
  pdf(paste0(resdir, "/PCA_15_top_genes_PC_", i, ".pdf"))
  for(j in 1:15){
    gene <- names(sort(pca_genes[,i])[j])
    print(FeaturePlot(sample, features=gene, reduction="pca"))
  }
  for(j in 1:15){
    gene <- names(sort(pca_genes[,i], decreasing = T)[j])
    print(FeaturePlot(sample, features=gene, reduction = "pca"))
  }

  dev.off()
}

# Project sample location on PC1 and PC2
pdf(file.path(resdir, "PCA_1_2_sample_location.pdf"))
sample$orig.ident_2959T <- sample$orig.ident_2960T <- sample$orig.ident_3133T <- sample$orig.ident_3377T <- sample$orig.ident_3610T <- sample$orig.ident_3662T <- "other"
sample@meta.data[grep("CHC2959T", rownames(sample@meta.data)),"orig.ident_2959T"] <- "CHC2959T"
sample@meta.data[grep("CHC2960T", rownames(sample@meta.data)),"orig.ident_2960T"] <- "CHC2960T"
sample@meta.data[grep("CHC3133T", rownames(sample@meta.data)),"orig.ident_3133T"] <- "CHC3133T"
sample@meta.data[grep("CHC3377T", rownames(sample@meta.data)),"orig.ident_3377T"] <- "CHC3377T"
sample@meta.data[grep("CHC3610T", rownames(sample@meta.data)),"orig.ident_3610T"] <- "CHC3610T"
sample@meta.data[grep("CHC3662T", rownames(sample@meta.data)),"orig.ident_3662T"] <- "CHC3662T"
DimPlot(sample, reduction="pca", dims=c(1,2), group.by = "orig.ident_2959T", cols = c("grey", "#CD9600"), order=c("CHC2959T", "other"))
DimPlot(sample, reduction="pca", dims=c(1,2), group.by = "orig.ident_2960T", cols = c("grey", "#7CAE00"), order=c("CHC2960T", "other"))
DimPlot(sample, reduction="pca", dims=c(1,2), group.by = "orig.ident_3133T", cols = c("grey", "#00BE67"), order=c("CHC3133T", "other"))
DimPlot(sample, reduction="pca", dims=c(1,2), group.by = "orig.ident_3377T", cols = c("grey", "#00A9FF"), order=c("CHC3377T", "other"))
DimPlot(sample, reduction="pca", dims=c(1,2), group.by = "orig.ident_3610T", cols = c("grey", "#C77CFF"), order=c("CHC3610T", "other"))
DimPlot(sample, reduction="pca", dims=c(1,2), group.by = "orig.ident_3662T", cols = c("grey", "#FF61CC"), order=c("CHC3662T", "other"))
dev.off()

# Split PCA bt sample  
pdf(file=file.path(resdir, "PCA_samples.pdf"))
DimPlot(sample, reduction="pca", split.by  = "orig.ident", group.by = "orig.ident")
dev.off()


############################
##### 3. VISUALIZATION #####
############################

DefaultAssay(sample) <- "SCT"
sample <- RunUMAP(sample, dims=1:nbcomp, verbose=F, umap.method = "uwot") # add UMAP for visualization

pdf(file=file.path(resdir, "UMAP_samples.pdf")) # project sample position on UMAP
DimPlot(sample, label=T, reduction="umap", group.by = "orig.ident") + NoLegend()
dev.off()

pdf(file=file.path(resdir, "UMAP_PC1_PC2_location.pdf")) # project PC1 and PC2 on UMAP
FeaturePlot(sample, features="PC_1")
FeaturePlot(sample, features="PC_2")
dev.off()


##############################################
##### DISPLAY FEATURES ON PCA PC1 vs PC2 #####
##############################################

##### H/LP/M signatures #####
#----------------------------

DefaultAssay(sample) <- "RNA"
sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor=10000) # use RNA lognormalized counts and not SCT when averaging gene expressions

mean_LP <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[LP_markers,], 2, mean); sample$mean_LP <- mean_LP
mean_H <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[H_markers,], 2, mean); sample$mean_H <- mean_H
mean_M <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[M_markers,], 2, mean); sample$mean_M <- mean_M

pdf(file.path(resdir, "PCA_H_LP_M_signatures.pdf"))
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
p1 <- FeaturePlot(sample, reduction="pca", dims=c(1,2), features = c("mean_H", "mean_LP", "mean_M"), order=T, combine = FALSE)
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
dev.off()

pdf(file.path(resdir, "UMAP_H_LP_M_signatures.pdf"))
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
p1 <- FeaturePlot(sample, reduction="umap", features = c("mean_H", "mean_LP", "mean_M"), order=T, combine = FALSE)
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
dev.off()


### correlation between components and signatures 

pca_coord <- sample@reductions$pca@cell.embeddings
signature_exp <- sample@meta.data[,c("mean_H", "mean_LP", "mean_M")]
identical(rownames(pca_coord), rownames(signature_exp)) # must be TRUE

d <- as.matrix(cor(pca_coord, signature_exp, method = "pearson"))
save(d, file=file.path(resdir, "correlation_PCA_coordinates_H_LP_M_signatures.RData"))

corrplot(d, order="original", tl.col="black", cl.cex=1.3, tl.cex=1.3)


###################
##### 3. SAVE #####
###################

DefaultAssay(sample) <- "SCT"

save(sample, file=file.path(resdir, "Seurat_object_analyses.RData"))
