
# strategy based on https://www.biorxiv.org/content/10.1101/2022.02.26.482041v1.full.pdf (The chromatin landscape of Th17cells reveals mechanisms of diversification of regulatory and pro inflammatory states, Thakore et al, 2022)


library(dplyr)
library(tidyr)
library(ArchR)
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

geco.changeRange <- function (v, newmin = 1, newmax = 10) 
{
  oldmin <- min(v, na.rm = TRUE)
  oldmax <- max(v, na.rm = TRUE)
  newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin))
}


######################
##### PARAMETERS #####
######################

rna_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/"
imputation_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/3.Denoising/Imputation_in_house_tumor_cells"
archrdir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/3.snATACseq/Final_analyses/All_samples_merged/"

grn_dir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/4.Integration/2.Gene_regulatory_network/Final_GRN"
resdir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/4.Integration/2.Gene_regulatory_network/Final_GRN/H_LP_M"


###########################
##### 0. DATA LOADING #####
###########################

### Load Multiome GRN
list_GRN <- geco.load(file.path(resdir, "GRN_list.RData")) # list of targets per TF
TF_ordered <- geco.load(file.path(resdir, "Heatmap_TF_targets_short_nb_targets_TF_ordered.RData"))

### TF defined by bulk/Multiome composite scores
df_TF_all <- geco.load(file.path(grn_dir, "TF_characteristics_bulk_snRNAseq_snATACseq_composite_scores.RData"))

### TF clusters
TF_ordered <- setdiff(TF_ordered, TF_ordered[which(TF_ordered=="TET1"):which(TF_ordered=="PRDM6")]) # remove mesenchymal TF
TF_LP <- TF_ordered[which(TF_ordered=="SOX4"):which(TF_ordered=="LHX1")]
TF_H_LP <- TF_ordered[which(TF_ordered=="MLXIPL"):which(TF_ordered=="MIXL1")]
TF_H <- TF_ordered[which(TF_ordered=="AR"):which(TF_ordered=="FOXP2")]
TF_GRN_cluster <- c(rep("LP", length(TF_LP)), rep("H_LP", length(TF_H_LP)), rep("H", length(TF_H))); names(TF_GRN_cluster) <- c(TF_LP, TF_H_LP, TF_H) # associate each TF to its GRN cluster

# Expression tables
rna_sample_all <- geco.load(file.path(rna_dir, "Seurat_object_analyses.RData"))
rna_exp <- rna_sample_all@assays$RNA@data
pca_cells <- rna_sample_all@reductions$pca@cell.embeddings # extract pca values per cell

# Imputated expression (cluster-based)
rna_exp_imputation <- geco.load(file.path(imputation_dir, "exp_mat_clusters_per_cell.RData")) 
imputation_clusters_annot <- geco.load(file.path(imputation_dir, "cluster_annot.RData"))
cells_imputation <- geco.load(file.path(imputation_dir, "cells_in_clusters.RData")) 

# Order cells based on PC1/2 values
cells_order_M <- grep("CHC3610T", colnames(rna_sample_all), value=T)
cells_order_H_LP <- colnames(rna_sample_all)[which(!colnames(rna_sample_all)%in%cells_order_M)]
cells_order_H_LP <- cells_order_H_LP[order(pca_cells[cells_order_H_LP,"PC_2"], decreasing = T)] # order H and LP cells by decreasing PC2 value
cells_order <- c(cells_order_H_LP) # no mesenchymal cells

# TF motif deviations, cisbp v1
TF_motif_dev_v1 <- geco.load(file.path(archrdir, "Motif_TF/MotifMatrix_cisbp.RData"))
TF_motif_dev_v1 <- TF_motif_dev_v1@assays@data[["z"]]
colnames(TF_motif_dev_v1) <- gsub("#", "_", gsub("-1", "", colnames(TF_motif_dev_v1))) # change to seurat nomenclature
rownames(TF_motif_dev_v1) <- paste0(gsub("\\_.*", "", rownames(TF_motif_dev_v1)), "_v1") # change nomenclature of TF names: remove the number after the underscore, add the version of cisbp
# compute cluster imputation for TF motif deviation (from snRNAseq clusters)
TF_motif_dev_imputation_v1 <- as.data.frame(matrix(NA, nrow=nrow(TF_motif_dev_v1), ncol=length(intersect(colnames(rna_exp_imputation), colnames(TF_motif_dev_v1))))) 
rownames(TF_motif_dev_imputation_v1) <- rownames(TF_motif_dev_v1); colnames(TF_motif_dev_imputation_v1) <- intersect(colnames(rna_exp_imputation), colnames(TF_motif_dev_v1))
for(cluster in imputation_clusters_annot$cluster){
  cells_in_cluster <- intersect(names(cells_imputation[which(cells_imputation==cluster)]), colnames(TF_motif_dev_imputation_v1)) 
  TF_motif_dev_imputation_v1[,cells_in_cluster] <- apply(TF_motif_dev_v1[,cells_in_cluster], 1, mean)
}
# TF motif deviations, cisbp v2
TF_motif_dev_v2 <- geco.load(file.path(archrdir, "Motif_TF/MotifMatrix_cisbp_v2.00.RData"))
TF_motif_dev_v2 <- TF_motif_dev_v2@assays@data[["z"]]
colnames(TF_motif_dev_v2) <- gsub("#", "_", gsub("-1", "", colnames(TF_motif_dev_v2))) # change to seurat nomenclature
rownames(TF_motif_dev_v2) <- paste0(gsub("\\_.*", "", rownames(TF_motif_dev_v2)), "_v2") # change nomenclature of TF names: remove the number after the underscore, add the version of cisbp
# compute cluster imputation for TF motif deviation (from snRNAseq clusters)
TF_motif_dev_imputation_v2 <- as.data.frame(matrix(NA, nrow=nrow(TF_motif_dev_v2), ncol=length(intersect(colnames(rna_exp_imputation), colnames(TF_motif_dev_v2))))) 
rownames(TF_motif_dev_imputation_v2) <- rownames(TF_motif_dev_v2); colnames(TF_motif_dev_imputation_v2) <- intersect(colnames(rna_exp_imputation), colnames(TF_motif_dev_v2))
for(cluster in imputation_clusters_annot$cluster){
  cells_in_cluster <- intersect(names(cells_imputation[which(cells_imputation==cluster)]), colnames(TF_motif_dev_imputation_v2)) 
  TF_motif_dev_imputation_v2[,cells_in_cluster] <- apply(TF_motif_dev_v2[,cells_in_cluster], 1, mean)
}
TF_motif_dev <- rbind(TF_motif_dev_v1, TF_motif_dev_v2) # combine cisbp v1 and v2
TF_motif_dev_imputation <- rbind(TF_motif_dev_imputation_v1, TF_motif_dev_imputation_v2) # combine cisbp v1 and v2


#############################################
##### 1. TF ORDERING BY POSITION ON PC2 #####
#############################################

TF_ordered_PC2 <- rep(NA, length(TF_ordered)); names(TF_ordered_PC2) <- TF_ordered # create a vector of the genes in heatmap to contain their order based on global max expression 

for(gene in names(TF_ordered_PC2)){
  snexp_curve <- data.frame(cell_pos=1:length(cells_order), exp=rna_exp[gene, cells_order], stringsAsFactors = F)
  # in order to smooth the data to find the global maximum, we compute the mean of every 100 cells using aggregate function (aggregate computes summary stat of data subsets): 
  snexp_curve_50 <- aggregate(snexp_curve, list(rep(1:(nrow(snexp_curve) %/% 100 + 1), each = 100, len = nrow(snexp_curve))), mean) # link aggregated values to the corresponding cells
  snexp_curve_50 <- snexp_curve_50[,-2] # we remove 2nd column = cell position because it is now the median of the cell positions --> not interesting
  snexp_curve_50 <- snexp_curve_50[-nrow(snexp_curve_50),] # remove last row that contains less than the expected number of cells and can bias
  
  print(ggplot(snexp_curve_50, aes(x=Group.1, y=exp))+
          geom_point()+
          ggtitle(gene)+
          xlab("Groups of 50 cells ordered according to PC values")+
          ylab("mean snRNAseq exp (data slot)")+
          geom_vline(xintercept = which(snexp_curve_50[,"exp"]==max(snexp_curve_50[,"exp"])), col="red", linetype="dashed"))
  
  # assign to the gene the cell group position where the maximum of the expression curve is reached
  TF_ordered_PC2[gene] <- which(snexp_curve_50[,"exp"]==max(snexp_curve_50[,"exp"]))
}

df_TF_ordered_PC2_with_value <- data.frame(TF=names(TF_ordered_PC2), value=TF_ordered_PC2) # order TFs by name when same value
df_TF_ordered_PC2_with_value <- df_TF_ordered_PC2_with_value %>% arrange(value, TF)

TF_ordered_PC2 <- df_TF_ordered_PC2_with_value$TF # order the genes by order of expression peaks along groups of 50 cells following PC correlation
save(TF_ordered_PC2, file=file.path(resdir, "Heatmaps_TF/TF_ordered_by_PC2_position.RData"))


#################################################
##### 2. MULTIOME TF EXPRESSION IN SNRNASEQ #####
#################################################

##### Heatmap #####
#------------------

DefaultAssay(rna_sample_all) <- "RNA" # set RNA assay for scaling and heatmap
rna_sample_all <- ScaleData(rna_sample_all, features=TF_ordered_PC2) # need to scale for Doheatmap
rna_sample_all$id="cells" # we need to define a column to group by on, otherwise DoHeatmap will cluster the cells by identity

# Smooth gene expression with bins of 100 cells
rna_exp_aggr <- as.data.frame(matrix(NA, ncol=length(cells_order), nrow=length(TF_ordered_PC2))); rownames(rna_exp_aggr) <- names(TF_ordered_PC2); colnames(rna_exp_aggr) <- cells_order

for(gene in TF_ordered_PC2){
  snexp_curve <- data.frame(cell_pos=1:length(cells_order), exp=rna_exp[gene, cells_order], stringsAsFactors = F)
  # in order to smooth the data to find the global maximum, we compute the mean of every 100 cells using aggregate function (aggregate computes summary stat of data subsets): 
  snexp_curve_50 <- aggregate(snexp_curve, list(rep(1:(nrow(snexp_curve) %/% 100 + 1), each = 100, len = nrow(snexp_curve))), mean) # link aggregated values to the corresponding cells
  rna_exp_aggr[gene, cells_order] <- rep(snexp_curve_50$exp, each = 100, len = nrow(snexp_curve))
}

pdf(file.path(resdir, "Heatmaps_TF/Heatmap_snRNAseq_TF_ordered_PC1_PC2.pdf"))
### raw RNA data values
DoHeatmap(rna_sample_all, features=TF_ordered_PC2, cells=cells_order, group.by = "id", size = 2) + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.y = element_text(size = 5))

### Use imputed snRNAseq data
rna_exp_imputed <- rna_exp_aggr[TF_ordered_PC2, cells_order]
rna_exp_imputed_scaled <- scale(as.matrix(t(rna_exp_imputed))) # scale the data to better see the expression changes (scale performs on each column)
rna_exp_imputed_scaled[which(rna_exp_imputed_scaled < -3)] <- -3; rna_exp_imputed_scaled[which(rna_exp_imputed_scaled > 3)] <- 3 # set minimum/maximum values to remove outliers

# Cell annotation
ha=ComplexHeatmap::HeatmapAnnotation(sample=rna_sample_all@meta.data[cells_order,"orig.ident"],
                                     subtype=rna_sample_all@meta.data[cells_order,"H_LP_M_clean_groups_1"],
                                     col=list(sample=c("CHC2959T"="#F87B72", "CHC2960T"="#B9A106", "CHC3133T"= "#00BA37", "CHC3377T"="#0EC2C6", "CHC3610T"= "#659FFF", "CHC3662T"="#F564E3"), subtype=c("H"="#E8AAD6", "H+LP"="#BC53CD", "LP"="#510787", "M"="#6F3323")))
# TF annotation
ha2=ComplexHeatmap::rowAnnotation(GRN_group=TF_GRN_cluster[TF_ordered_PC2], 
                                  col=list(GRN_group=c("H"="#E8AAD6", "H_LP"="#BC53CD", "LP"="#510787", "M"="#6F3323")))
par(mar=c(3,4,2,2)) 
ht <- draw(ComplexHeatmap::Heatmap(t(rna_exp_imputed_scaled), cluster_rows=F, cluster_columns=F, show_column_names=F, show_row_names=T, use_raster=F, top_annotation = ha, right_annotation = ha2, row_names_gp = gpar(fontsize = 6)),  
           heatmap_legend_side = "bottom") # here we do not rasterize to have all cells displayed ant not a more global view 
dev.off()


##### Individual graphs on PCA #####
#-----------------------------------

### Create PCA dataframe
df_pca <- as.data.frame(rna_sample_all@reductions$pca@cell.embeddings[,c("PC_1", "PC_2")]) # extract PC cell coordinates
if(identical(rownames(df_pca), rownames(rna_sample_all@meta.data))){ # add cell annotations
  df_pca$cell <- rownames(rna_sample_all@meta.data)
  df_pca$sample <- rna_sample_all@meta.data$orig.ident
  df_pca$subtype <- rna_sample_all@meta.data$H_LP_M_clean_groups_1}
df_pca <- df_pca %>% filter(sample != "CHC3610T") # here we do not show M cels because focus on H LP transition = 11832 cells
df_pca <- df_pca[sample(rownames(df_pca)),] # shuffle randomly the cell position in the table to have random overlap of the points in the grahp
df_pca$groups_by <- "all" # create column for jitter plot with only one group

pdf(file.path(resdir, "Heatmaps_TF/Plot_snRNAseq_TF_on_PCA.pdf"), width=5, height=3)
### Plot subtypes on graph
ggplot(df_pca, aes(x=sample, y=PC_2, col=sample))+
  geom_jitter(size=1.5, width =length(unique(df_pca$sample))/20)+
  theme_classic()+
  scale_color_manual(values=c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")
ggplot(df_pca, aes(x=groups_by, y=PC_2, col=PC_2))+
  geom_jitter(size=1, width =length(unique(df_pca$groups_by))/10)+
  theme_classic()+
  scale_color_gradientn(colours = c("#F3D5EA", "#E8AAD6", "#BC53CD", "#510787", "#1B002A"), values=scales::rescale(c(-50,-25,0,25,50)))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")

### Show TF expression
for(TF in TF_ordered_PC2){
  df_pca$TF_exp <- rna_exp_imputation[TF, df_pca$cell] # take imputed expression for easier visualization
  df_pca <- df_pca[sample(rownames(df_pca)),] # shuffle randomly the cell position in the table to have random overlap of the points in the grahp
  print(ggplot(df_pca, aes(x=sample, y=PC_2, col=TF_exp))+
          geom_jitter(size=1.5, width =length(unique(df_pca$sample))/10)+
          theme_classic()+
          scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
          theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
          ggtitle(paste0(TF, " expression (snRNA-seq)")))
}
dev.off()


########################################################
##### 3. MULTIOME TF TARGET EXPRESSION IN SNRNASEQ #####
########################################################

##### Heatmap #####
#------------------

# create mean expression matrix for targets for each TF
exp_targets <- matrix(NA, nrow=length(TF_ordered_PC2), ncol=length(cells_order)); rownames(exp_targets) <- TF_ordered_PC2; colnames(exp_targets) <- cells_order
for(TF in TF_ordered_PC2){
  targets <- list_GRN[[TF]]
  exp_targets[TF, cells_order] <- apply(rna_exp[intersect(targets, rownames(rna_exp)), cells_order], 2, mean)
}

# smooth gene expression by aggregating in bins of 100 cells
rna_exp_aggr <- as.data.frame(matrix(NA, ncol=length(cells_order), nrow=length(TF_ordered_PC2))); rownames(rna_exp_aggr) <- TF_ordered_PC2; colnames(rna_exp_aggr) <- cells_order
for(gene in TF_ordered_PC2){
  snexp_curve <- data.frame(cell_pos=1:length(cells_order), exp=exp_targets[gene, cells_order], stringsAsFactors = F)
  # in order to smooth the data to find the global maximum, we compute the mean of every 50 cells using aggregate function (aggregate computes summary stat of data subsets): 
  snexp_curve_50 <- aggregate(snexp_curve, list(rep(1:(nrow(snexp_curve) %/% 100 + 1), each = 100, len = nrow(snexp_curve))), mean) # link aggregated values to the corresponding cells
  rna_exp_aggr[gene, cells_order] <- rep(snexp_curve_50$exp, each = 100, len = nrow(snexp_curve))
}


pdf(file.path(resdir, "Heatmaps_TF/Heatmap_snRNAseq_targets_ordered_PC1_PC2.pdf"))

### Use imputed snRNAseq data
rna_exp_imputed <- rna_exp_aggr[TF_ordered_PC2, cells_order]
rna_exp_imputed_scaled <- scale(as.matrix(t(rna_exp_imputed))) # scale the data to better see the expression changes (scale performs on each column)
rna_exp_imputed_scaled[which(rna_exp_imputed_scaled < -3)] <- -3; rna_exp_imputed_scaled[which(rna_exp_imputed_scaled > 3)] <- 3 # set minimum/maximum values to remove outliers

# Cell annotation
ha=ComplexHeatmap::HeatmapAnnotation(sample=rna_sample_all@meta.data[cells_order,"orig.ident"],
                                     subtype=rna_sample_all@meta.data[cells_order,"H_LP_M_clean_groups_1"],
                                     col=list(sample=c("CHC2959T"="#F87B72", "CHC2960T"="#B9A106", "CHC3133T"= "#00BA37", "CHC3377T"="#0EC2C6", "CHC3610T"= "#659FFF", "CHC3662T"="#F564E3"), subtype=c("H"="#E8AAD6", "H+LP"="#BC53CD", "LP"="#510787", "M"="#6F3323")))
# TF annotation
ha2=ComplexHeatmap::rowAnnotation(GRN_group=TF_GRN_cluster[TF_ordered_PC2], 
                                  col=list(GRN_group=c("H"="#E8AAD6", "H_LP"="#BC53CD", "LP"="#510787", "M"="#6F3323")))
par(mar=c(3,4,2,2)) 
ht <- draw(ComplexHeatmap::Heatmap(t(rna_exp_imputed_scaled), cluster_rows=F, cluster_columns=F, show_column_names=F, show_row_names=T, use_raster=F, top_annotation = ha, right_annotation = ha2, row_names_gp = gpar(fontsize = 6)),  
           heatmap_legend_side = "bottom") # here we do not rasterize to have all cells displayed ant not a more global view 
dev.off()


##### Individual graphs on PCA #####
#-----------------------------------

### Create PCA dataframe
df_pca <- as.data.frame(rna_sample_all@reductions$pca@cell.embeddings[,c("PC_1", "PC_2")]) # extract PC cell coordinates
if(identical(rownames(df_pca), rownames(rna_sample_all@meta.data))){ # add cell annotations
  df_pca$cell <- rownames(rna_sample_all@meta.data)
  df_pca$sample <- rna_sample_all@meta.data$orig.ident
  df_pca$subtype <- rna_sample_all@meta.data$H_LP_M_clean_groups_1}
df_pca <- df_pca %>% filter(sample != "CHC3610T") # here we do not show M cels because focus on H LP transition = 11832 cells
df_pca <- df_pca[sample(rownames(df_pca)),] # shuffle randomly the cell position in the table to have random overlap of the points in the grahp

pdf(file.path(resdir, "Heatmaps_TF/Plot_snRNAseq_TF_targets_on_PCA.pdf"))
### Plot subtypes on graph
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=sample))+
  geom_point(size=1)+
  theme_classic()+
  scale_color_manual(values=c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=PC_2))+
  geom_point(size=1)+
  theme_classic()+
  scale_color_gradientn(colours = c("#F3D5EA", "#E8AAD6", "#BC53CD", "#510787", "#1B002A"), values=scales::rescale(c(-50,-25,0,25,50)))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")

### Show TF expression
for(TF in TF_ordered_PC2){
  targets <- list_GRN[[TF]] # identify TF targets
  df_pca$target_exp <- apply(rna_exp_imputation[intersect(targets, rownames(rna_exp_imputation)), df_pca$cell], 2, mean) # take imputed expression for easier visualization
  df_pca <- df_pca[sample(rownames(df_pca)),] # shuffle randomly the cell position in the table to have random overlap of the points in the grahp
  print(ggplot(df_pca, aes(x=PC_1, y=PC_2, col=target_exp))+
          geom_point(size=1)+
          theme_classic()+
          scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
          theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
          ggtitle(paste0(TF, " target expression (snRNA-seq)")))
}
dev.off()


#################################################
##### 4. MULTIOME TF DEVIATION IN SNATACSEQ #####
#################################################

# Define cisbp version for each TF
top_TF_motifs <- c()
for(g in TF_ordered_PC2){ # add the corresponding version to the TF name
  if(df_TF_all[g, "cisbp_version"] == "v1"){top_TF_motifs <- c(top_TF_motifs, paste0(g, "_v1"))}
  if(df_TF_all[g, "cisbp_version"] == "v2"){top_TF_motifs <- c(top_TF_motifs, paste0(g, "_v2"))}
}

# smooth TF motif deciations by aggregating in bins of 100 cells
atac_exp_aggr <- as.data.frame(matrix(NA, ncol=length(intersect(cells_order, colnames(TF_motif_dev))), nrow=length(TF_ordered_PC2))); rownames(atac_exp_aggr) <- TF_ordered_PC2; colnames(atac_exp_aggr) <- intersect(cells_order, colnames(TF_motif_dev))
for(gene in TF_ordered_PC2){
  snexp_curve <- data.frame(cell_pos=1:length(intersect(cells_order, colnames(TF_motif_dev))), exp=TF_motif_dev[top_TF_motifs[which(TF_ordered_PC2==gene)], intersect(cells_order, colnames(TF_motif_dev))], stringsAsFactors = F)
  # in order to smooth the data to find the global maximum, we compute the mean of every 50 cells using aggregate function (aggregate computes summary stat of data subsets): 
  snexp_curve_50 <- aggregate(snexp_curve, list(rep(1:(nrow(snexp_curve) %/% 100 + 1), each = 100, len = nrow(snexp_curve))), mean) # link aggregated values to the corresponding cells
  atac_exp_aggr[gene, intersect(cells_order, colnames(TF_motif_dev))] <- rep(snexp_curve_50$exp, each = 100, len = nrow(snexp_curve))
}

##### Heatmap #####
#------------------

DefaultAssay(rna_sample_all) <- "RNA" # set RNA assay for scaling and heatmap
rna_sample_all <- ScaleData(rna_sample_all, features=TF_ordered_PC2) # need to scale for Doheatmap
rna_sample_all$id="cells" # we need to define a column to group by on, otherwise DoHeatmap will cluster the cells by identity

pdf(file.path(resdir, "Heatmaps_TF/Heatmap_snRNAseq_TF_motif_deviations_PC1_PC2.pdf"))


### Use imputed snATACseq motif data
motifs. <- atac_exp_aggr[TF_ordered_PC2, intersect(cells_order, colnames(atac_exp_aggr))] # restrict to TF under study
motifs. <- scale(as.matrix(t(motifs.))) # scale the data to better see the expression changes (scale performs on each column)
motifs.[which(motifs. < -4)] <- -4; motifs.[which(motifs. > 4)] <- 4 # set minimum/maximum values to remove outliers

# Cell annotation
ha=ComplexHeatmap::HeatmapAnnotation(sample=rna_sample_all@meta.data[intersect(cells_order, colnames(atac_exp_aggr)),"orig.ident"],
                                     subtype=rna_sample_all@meta.data[intersect(cells_order, colnames(atac_exp_aggr)),"H_LP_M_clean_groups_1"],
                                     col=list(sample=c("CHC2959T"="#F87B72", "CHC2960T"="#B9A106", "CHC3133T"= "#00BA37", "CHC3377T"="#0EC2C6", "CHC3610T"= "#659FFF", "CHC3662T"="#F564E3"), subtype=c("H"="#E8AAD6", "H+LP"="#BC53CD", "LP"="#510787", "M"="#6F3323")))
# TF annotation
ha2=ComplexHeatmap::rowAnnotation(GRN_group=TF_GRN_cluster[TF_ordered_PC2], 
                                  col=list(GRN_group=c("H"="#E8AAD6", "H_LP"="#BC53CD", "LP"="#510787", "M"="#6F3323")))
par(mar=c(3,4,2,2)) 
ht <- draw(ComplexHeatmap::Heatmap(t(motifs.), cluster_rows=F, cluster_columns=F, show_column_names=F, show_row_names=T, use_raster=F, top_annotation = ha, right_annotation = ha2, row_names_gp = gpar(fontsize = 6), circlize::colorRamp2(c(-2, 0, 2), c("darkblue", "grey95", "gold"))),  
           heatmap_legend_side = "bottom",) # here we do not rasterize to have all cells displayed ant not a more global view 
dev.off()


##### Individual graphs on PCA #####
#-----------------------------------

### Create PCA dataframe
df_pca <- as.data.frame(rna_sample_all@reductions$pca@cell.embeddings[,c("PC_1", "PC_2")]) # extract PC cell coordinates
if(identical(rownames(df_pca), rownames(rna_sample_all@meta.data))){ # add cell annotations
  df_pca$cell <- rownames(rna_sample_all@meta.data)
  df_pca$sample <- rna_sample_all@meta.data$orig.ident
  df_pca$subtype <- rna_sample_all@meta.data$H_LP_M_clean_groups_1}
df_pca <- df_pca %>% filter(sample != "CHC3610T") # here we do not show M cels because focus on H LP transition = 11832 cells
df_pca <- df_pca[sample(rownames(df_pca)),] # shuffle randomly the cell position in the table to have random overlap of the points in the graph
df_pca <- df_pca[which(df_pca$cell %in% colnames(TF_motif_dev_imputation)),] # restrict to cells in ATAC = 8975 cells 

pdf(file.path(resdir, "Heatmaps_TF/Plot_snATACseq_TF_deviations_on_PCA.pdf"))
### Plot subtypes on graph
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=sample))+
  geom_point(size=1)+
  theme_classic()+
  scale_color_manual(values=c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=PC_2))+
  geom_point(size=1)+
  theme_classic()+
  scale_color_gradientn(colours = c("#F3D5EA", "#E8AAD6", "#BC53CD", "#510787", "#1B002A"), values=scales::rescale(c(-50,-25,0,25,50)))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")

### Show TF expression
for(TF in TF_ordered_PC2){
  df_pca$TF_dev <- as.numeric(TF_motif_dev_imputation[top_TF_motifs[which(TF_ordered_PC2==TF)], df_pca$cell]) # take imputed expression for easier visualization
  df_pca <- df_pca[sample(rownames(df_pca)),] # shuffle randomly the cell position in the table to have random overlap of the points in the grahp
  print(ggplot(df_pca, aes(x=PC_1, y=PC_2, col=TF_dev))+
          geom_point(size=1)+
          theme_classic()+
          scale_color_gradient2( low = "blue", mid = "white", high="red")+
          theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
          ggtitle(paste0(TF, " motif deviation (snATAC-seq)")))
}
dev.off()


###########################################################
##### 5. MULTIOME TF EXPRESSION VS TF MOTIF DEVIATION #####
###########################################################

pdf(file.path(resdir, "Heatmaps_TF/Curves_TF_expression_motif_deviation_along_PC2.pdf"), width=6, height=3)
for(gene in TF_ordered_PC2){
  if(df_TF_all[gene, "cisbp_version"] == "v1"){gene_motif <- c(paste0(gene, "_v1"))} # define TF motif
  if(df_TF_all[gene, "cisbp_version"] == "v2"){gene_motif <- c(paste0(gene, "_v2"))}
  
  TF_exp <- rna_exp[gene, intersect(cells_order, colnames(TF_motif_dev))] # TF expression in common cells RNA/ATAC
  TF_targets <- apply(rna_exp[unique(unlist(list_GRN[gene])), intersect(cells_order, colnames(TF_motif_dev))], 2, mean) # TF expression in common cells RNA/ATAC
  TF_motif <- TF_motif_dev[gene_motif, intersect(cells_order, colnames(TF_motif_dev))] # TF motif deviation in common cells RNA/ATAC
  
  df_rna_atac <- data.frame(cell=names(TF_exp), PC2_value=pca_cells[names(TF_exp),"PC_2"], TF_exp=TF_exp, TF_targets=TF_targets, TF_motif=TF_motif, stringsAsFactors = F) # create dataframe containing both TF expression and motif deviation
  # aggregate into bins of 100 cells along PC2
  df_rna_atac_100 <- aggregate(df_rna_atac, list(rep(1:(nrow(df_rna_atac) %/% 200 + 1), each = 200, len = nrow(df_rna_atac))), mean) # link aggregated values to the corresponding cells
  df_rna_atac_100$cell <- NULL # remove column of cells which only contains NA
  # set same scale for TF exp and motif deviation (we just want to see the maximum location)
  df_rna_atac_100$TF_exp <- geco.changeRange(df_rna_atac_100$TF_exp, newmin=0, newmax=1)
  df_rna_atac_100$TF_targets <- geco.changeRange(df_rna_atac_100$TF_targets, newmin=0, newmax=1)
  df_rna_atac_100$TF_motif <- geco.changeRange(df_rna_atac_100$TF_motif, newmin=0, newmax=1)
  
  print(ggplot(df_rna_atac_100) + stat_smooth(aes(x=PC2_value, y=TF_exp), method = "gam", se = F, col="red", size=2)+
          stat_smooth(aes(x=PC2_value, y=TF_targets), method = "gam", se = F, col="orange", size=2)+theme_classic()+
          stat_smooth(aes(x=PC2_value, y=TF_motif), method = "gam", se = F, col="blue", size=2)+
      ggtitle(gene))
}
dev.off()





