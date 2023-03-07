
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(circlize)

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

resdir = paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/Differential_expression")

include_GLI3 = FALSE # do we include GLI3 cluster in CHC3610T?

###########################
##### 0. DATA LOADING #####
###########################

sample_seurat <- geco.load(file.path(datadir, "Seurat_object_analyses.RData"))


######################################
##### 1. DIFFERENTIAL EXPRESSION #####
######################################

DefaultAssay(sample_seurat) <- "RNA" # the assay to use for differential expression is RNA and not SCT (recommended by satijalab)

### Remove GLI3 cluster from analyses

sample_seurat$H_LP_M <- as.character(sample_seurat$H_LP_M_clean_groups_1)
if(include_GLI3==FALSE){sample_seurat@meta.data[which(sample_seurat@meta.data$SCT_snn_res.0.3==15), "H_LP_M"] <- "M_GLI3_cluster"}
DimPlot(sample_seurat, group.by = "H_LP_M")

Idents(sample_seurat) <- "H_LP_M"

### H vs LP
diff_exp_H_vs_LP <- FindMarkers(sample_seurat, ident.1 = "H", ident.2 = "LP", logfc.threshold=0, min.pct = 0) # findmarkers can compute logFC even when gene not detected in a population because avg logFC = differences of the log(mean(expm1(x)+1))
save(diff_exp_H_vs_LP, file=file.path(resdir, "diff_exp_H_vs_LP.RData"))

### M vs H+LP
diff_exp_M_vs_H_LP <- FindMarkers(sample_seurat, ident.1 = "M", ident.2 = c("H", "H+LP", "LP"), logfc.threshold=0, min.pct = 0)
save(diff_exp_M_vs_H_LP, file=file.path(resdir, "diff_exp_M_vs_H_LP.RData"))

