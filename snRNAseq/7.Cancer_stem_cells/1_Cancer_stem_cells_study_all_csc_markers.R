
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(ggpubr)
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

datadir = paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/")

resdir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/10.Cancer_stem_cells/"; if(!dir.exists(resdir)){dir.create(resdir)}

### Cancer stem cell markers (found in literature,liver specific)

#csc_markers <- unique(c("PROM1", "THY1", "EPCAM", "CD24", "ANPEP", "CD34", "SOX9", "ABCG2", "CD44", "ALDH1A1", "ALDH3A1", "KRT19", "SOX12", "CD47",
#                        "CD24", "CD44", "CD90", "CD133", "EPCAM", "AFP", "NANOG", "NOTCH", "OCT3", "OCT4", "SOX2",
#                        "CD133", "ALDH", "CD90", "CD44", "CD13", "EPCAM", "OC6", "1B50-1", "SALL4", "ICAM1",
#                        "POU5F1"))

csc_markers <- unique(c("PROM1", "THY1", "EPCAM", "CD24", "ANPEP", "CD34", "SOX9", "ABCG2", "CD44", "ALDH1A1", "ALDH3A1", "KRT19", "SOX12", "CD47",
                 "BMI1", "GLI1", "LGR5", "NANOG", "SHH", "SOX2",
                 "KIT", "MUC1", "DLK1", "ICAM1", "SALL4", "KLF4", "CD19", "MS4A1", "DPP4", "CD27", "CD38", "PTPRC", "ITGA6", "CEACAM6", "ALCAM", "NGFR", "ENG", "CD177", "CD151",
                 "NES", "MSI1", "HAVCR2", "POU5F1", "MYC", "CXCR4", "CXCL12", "PLAUR", "ITGB1", "ITGB3", "CD70", "LETM1", "GLI2", "LINGO2", "AFP", "NOTCH1", "MSI2",
                 "CD33", "IL3RA", "CLEC12A", "IL2RA", "CD36", "KIT", "IL1RAP", "JAK2", "STAT1", "PODXL", "FUT4", "PROCR")) # 67 markers


### Define UMI counts threshold to keep cells with sufficient expression
umi_threshold = 20


###########################
##### 0. DATA LOADING #####
###########################

### Load seurat object for merged tumor samples (tumor cells only)
sample_seurat <- geco.load(file.path(datadir, "Seurat_object_analyses.RData"))
DefaultAssay(sample_seurat) <- "RNA" # switch to "RNA" assay for visualization
sample_seurat <- ScaleData(sample_seurat, features=csc_markers) # scale for the csc markers for heatmap

### Retrieve PCA values for most variable genes
pca_values_all <-sample_seurat@reductions$pca@feature.loadings
pca_values <- pca_values_all[intersect(csc_markers, rownames(pca_values_all)),1:50] # 20 genes in variable features, take the top 50 components

### restrict to CSC markers present in seurat object
setdiff(csc_markers, rownames(sample_seurat)) # check which markers are not present (sometimes it is a different gene name used)

csc_markers <- intersect(csc_markers, rownames(sample_seurat@assays$RNA@counts)) # 67 genes


#############################################
##### 1. MERGED UMAP OF ALL CSC MARKERS #####
#############################################

pdf(file.path(resdir, "Cancer_stem_cell_markers_umap_pca.pdf"))
csc_markers <- intersect(csc_markers, rownames(sample_seurat)) #53 markers
DimPlot(sample_seurat, group.by = "H_LP_M_clean_groups_1")
for(g in csc_markers){
  print(FeaturePlot(sample_seurat, features=g, reduction="umap") + scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20')))
  print(FeaturePlot(sample_seurat, features=g, reduction="pca") + scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20')))
}
dev.off()


##############################
##### 2. ALL CSC MARKERS ##### 
##############################

##### 1. Restrict exp data to CSC markers #####
#----------------------------------------------

counts_csc <- sample_seurat@assays$SCT@counts[csc_markers,] # choose SCT corrected counts to avoid depth sequencing bias for cell selection
exp_csc <- sample_seurat@assays$RNA@data[csc_markers,] # but to display expression we stick to RNA assay


##### 2. Keep cells with a certain amount of UMI counts for > 1 CSC marker #####
#-------------------------------------------------------------------------------

is_marker_detected <- apply(counts_csc, 2, function(x) length(x[x>=umi_threshold])>=1)
table(is_marker_detected) # for umi_threshold=3: ~6000 potential cancer stem cells, for umi_threshold=5: ~2400 potential cancer stem cells


##### 3. Visualize those cells on UMAP to check their consistency #####
#----------------------------------------------------------------------

### Create column in seurat object
sample_seurat$potential_csc <- is_marker_detected
potential_csc <- colnames(sample_seurat)[which(sample_seurat$potential_csc==T)] # retrieve potential cancer stem cells

### Visualize
DimPlot(sample_seurat, group.by = "potential_csc", reduction="umap")
DimPlot(sample_seurat, group.by = "potential_csc", reduction="pca")
DimPlot(sample_seurat, group.by = "orig.ident", reduction="pca")

# compare with depth distribution
sample_seurat$log_ncounts_RNA <- log10(sample_seurat$nCount_RNA)
FeaturePlot(sample_seurat, features="log_ncounts_RNA", reduction="umap")
FeaturePlot(sample_seurat, features="log_ncounts_RNA", reduction="pca")


##### 4. Clustering of potential cancer stem cells #####
#-------------------------------------------------------

mat <- sample_seurat@assays$RNA@scale.data[csc_markers,potential_csc] # select input matrix = here the scaled data

### Annotations
col_ha <- HeatmapAnnotation(sample=sample_seurat@meta.data[potential_csc,"orig.ident"], 
                            subtype=sample_seurat@meta.data[potential_csc,"H_LP_M_clean_groups_1"],
                            col = list(sample=c("CHC2959T"="#f8766d", "CHC2960T"="#b79f00", "CHC3133T"="#00ba38", "CHC3377T"="#00bfc4", "CHC3610T"="#619bff", "CHC3662T"="#f564e4"),
                                       subtype=c("H"="#E8AAD6", "H+LP"="#BC53CD", "LP"="#510787", "M"="brown")))

pdf(paste0(resdir, "All_CSC_markers_complexHeatmap_min_", umi_threshold, "_UMI_SCT.pdf"))
ht <- draw(Heatmap(as.matrix(mat), show_column_names = F,show_row_names=T, clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D",
           top_annotation=col_ha, row_names_gp = grid::gpar(fontsize = 5)))
dev.off()

