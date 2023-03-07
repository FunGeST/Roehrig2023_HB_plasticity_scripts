

# strategy based on https://www.biorxiv.org/content/10.1101/2022.02.26.482041v1.full.pdf (The chromatin landscape of Th17cells reveals mechanisms of diversification of regulatory and pro inflammatory states, Thakore et al, 2022)

library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(corrplot)
library(dendextend)
set.seed(1234)

geco.load <- function (filename) # File to be loaded
{
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error - function geco.load : file ", filename, 
            " not found"))
  NULL
}

### parameters

rna_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/"
resdir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/4.Integration/2.Gene_regulatory_network/Final_GRN/H_LP_M"


### data loading

cor_mat <- geco.load(file.path(resdir, "Correlation_expression_TF_vs_target_genes.RData")) # load matrix of snRNAseq correlations between TF and targets
list_GRN <- geco.load(file.path(resdir, "GRN_list.RData")) # load list of GRN

TF_ordered <- geco.load(file.path(resdir, "Heatmap_TF_targets_short_nb_targets_TF_ordered.RData")) # load TF ordered on heatmap
TF_M <- TF_ordered[which(TF_ordered=="TET1"):which(TF_ordered=="PRDM6")]
TF_LP <- TF_ordered[which(TF_ordered=="SOX4"):which(TF_ordered=="LHX1")]
TF_H_LP <- TF_ordered[which(TF_ordered=="MLXIPL"):which(TF_ordered=="MIXL1")]
TF_H <- TF_ordered[which(TF_ordered=="AR"):which(TF_ordered=="FOXP2")]
targets_M_all <- unique(unlist(list_GRN[TF_M])); length(targets_M_all)
targets_LP_all <- unique(unlist(list_GRN[TF_LP])); length(targets_LP_all)
targets_H_LP_all <- unique(unlist(list_GRN[TF_H_LP])); length(targets_H_LP_all)
targets_H_all <- unique(unlist(list_GRN[TF_H])); length(targets_H_all)

rna_sample_all <- geco.load(file.path(rna_dir, "Seurat_object_analyses.RData")) # load Seurat snRNAseq object
rna_exp <- as.matrix(GetAssayData(rna_sample_all, assay="RNA", slot="data")) # take RNA data assay (log normalized) for expression visualization
CHC3610T_cells <- grep("CHC3610T", colnames(rna_exp), value=T) # isolate name of mesenchymal cells to remove for some correlations (ex PC2 vs TF) 
cells_H <- rownames(rna_sample_all@meta.data[which(rna_sample_all$H_LP_M_clean_groups_1=="H"),])
cells_LP <- rownames(rna_sample_all@meta.data[which(rna_sample_all$H_LP_M_clean_groups_1=="LP"),])
cells_M <- rownames(rna_sample_all@meta.data[which(rna_sample_all$H_LP_M_clean_groups_1=="M"),])


### correlation matrix

cor_mat <- cor_mat[which(rownames(cor_mat) %in% unique(unlist(list_GRN))),] # restrict to genes from the GRN
write.table(cor_mat, file.path(resdir, "Correlation_expression_TF_vs_target_genes.csv"), sep=";", col.names=T, row.names=T) # load matrix of snRNAseq correlations between TF and targets


### mean expression per TF or target

df <- data.frame(gene=unique(c(colnames(cor_mat), rownames(cor_mat))), TF_or_target=NA, GRN_module=NA, exp_H=NA, exp_LP=NA, exp_M=NA)
df$TF_or_target <- c(rep("TF", times=ncol(cor_mat)), rep("target", times=nrow(cor_mat))) # describe if gene is TF or target 
df$exp_H <- apply(rna_exp[df$gene, cells_H], 1, mean) # add mean expression in scH, LP and M cells for each gene
df$exp_LP <- apply(rna_exp[df$gene, cells_LP], 1, mean)
df$exp_M <- apply(rna_exp[df$gene, cells_M], 1, mean)
# assign each gene to module
df$H_module <- ifelse(df$gene %in% TF_H | df$gene %in% targets_H_all, "scH", NA)
df$H_LP_module <- ifelse(df$gene %in% TF_H_LP | df$gene %in% targets_H_LP_all, "sc-epi", NA)
df$LP_module <- ifelse(df$gene %in% TF_LP | df$gene %in% targets_LP_all, "scLP", NA)
df$M_module <- ifelse(df$gene %in% TF_M | df$gene %in% targets_M_all, "scM", NA)
df <- df %>% unite("GRN_module", c("H_module", "H_LP_module", "LP_module", "M_module"), sep=", ", na.rm=T)

write.table(df, file.path(resdir, "Mean_expression_in_cell_populations.csv"), sep=";", col.names=T, row.names=F) # load matrix of snRNAseq correlations between TF and targets
