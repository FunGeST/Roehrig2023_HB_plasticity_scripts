
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

# retrieve CNV annotation:
cnv_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/"

### Cancer stem cell markers (found in literature,liver specific)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4814640/#:~:text=Liver%20cancer%20stem%20cell%20markers%20have%20been%20defined%20that%20can,potential%20to%20serve%20as%20CSCs.
# https://www.frontiersin.org/articles/10.3389/fgene.2020.00112/full
# https://www.sciencedirect.com/science/article/pii/S1873506122000502
# https://journals.sagepub.com/doi/10.1177/1758835918816287

csc_markers_liver <- c("PROM1", "THY1", "CD44", "ANPEP", "EPCAM", "CD24", "ALDH3A1", "KRT19", "ICAM1", "CD34", "SOX9", "ABCG2", "SOX12", "CD47") # originally there was ALDH1A1 but it was expressed in too many cells
csc_markers_liver_main <- c("PROM1", "THY1", "CD44", "EPCAM", "CD24") # liver CSC markers that are the most recurrent

### Define UMI counts threshold to keep cells with sufficient expression
umi_threshold = 5


###########################
##### 0. DATA LOADING #####
###########################

### Load seurat object for merged tumor samples (tumor cells only)
sample_seurat <- geco.load(file.path(datadir, "Seurat_object_analyses.RData"))
DefaultAssay(sample_seurat) <- "RNA" # switch to "RNA" assay for visualization
sample_seurat <- ScaleData(sample_seurat, features=csc_markers_liver) # scale for the csc markers for heatmap

### restrict to CSC markers present in seurat object
setdiff(csc_markers_liver, rownames(sample_seurat)) # check which markers are not present (sometimes it is a different gene name used)

csc_markers_liver <- intersect(csc_markers_liver, rownames(sample_seurat@assays$RNA@counts)) # 14 genes


##############################################
##### 1. CLUSTERING ON LIVER CSC MARKERS #####
##############################################

##### 1. Restrict exp data to CSC markers #####
#----------------------------------------------

counts_csc <- sample_seurat@assays$SCT@counts[csc_markers_liver,] # choose SCT corrected counts to avoid depth sequencing bias for cell selection
exp_csc <- sample_seurat@assays$RNA@data[csc_markers_liver,] # but to display expression we stick to RNA assay


##### 2. Keep cells with a certain amount of UMI counts for > 1 CSC marker #####
#-------------------------------------------------------------------------------

is_marker_detected <- apply(counts_csc, 2, function(x) length(x[x>=umi_threshold])>=1)
table(is_marker_detected) # for umi_threshold1=3: ~6000 potential cancer stem cells, for umi_threshold1=5: ~2400 potential cancer stem cells


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

mat <- sample_seurat@assays$RNA@scale.data[csc_markers_liver,potential_csc] # select input matrix = here the scaled data

### Annotations
col_ha <- HeatmapAnnotation(sample=sample_seurat@meta.data[potential_csc,"orig.ident"], 
                            subtype=sample_seurat@meta.data[potential_csc,"H_LP_M_clean_groups_1"],
                            col = list(sample=c("CHC2959T"="#f8766d", "CHC2960T"="#b79f00", "CHC3133T"="#00ba38", "CHC3377T"="#00bfc4", "CHC3610T"="#619bff", "CHC3662T"="#f564e4"),
                                       subtype=c("H"="#E8AAD6", "H+LP"="#BC53CD", "LP"="#510787", "M"="brown")))

pdf(paste0(resdir, "Liver_CSC_markers_complexHeatmap_min_", umi_threshold, "_UMI_SCT.pdf"))
ht <- draw(Heatmap(as.matrix(mat), show_column_names = F,show_row_names=T, clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D",
                   top_annotation=col_ha, col = circlize::colorRamp2(c(min(quantile(as.matrix(mat), probs=c(0.01,0.99))), 0, max(quantile(as.matrix(mat), probs=c(0.01,0.99)))), c("blue", "white", "red"))))
dev.off()


############################################
##### 2. POTENTIAL CSC FOR EACH MARKER #####
############################################

# count matrix
counts_csc <- sample_seurat@assays$SCT@counts[csc_markers_liver,] # choose SCT corrected counts to avoid depth sequencing bias for cell selection


##### 1. Identify the potential csc based on marker counts #####
#---------------------------------------------------------------

df_potential_csc <- as.data.frame(matrix(FALSE, ncol=ncol(sample_seurat), nrow=length(csc_markers_liver))) # create dataframe that will tell for each marker if the cell has strong expression or not
rownames(df_potential_csc) <- csc_markers_liver; colnames(df_potential_csc) <- colnames(sample_seurat)

pdf(paste0(resdir, "Liver_CSC_markers_distribution_samples_subtypes_min_", umi_threshold, "_UMI_SCT.pdf"))
for(marker in csc_markers_liver){
  # keep cells with > x UMI counts for this marker (= cells positive for this marker)
  cells_for_this_marker <- counts_csc[marker,]  
  cells_for_this_marker <- names(cells_for_this_marker[which(cells_for_this_marker>=umi_threshold)])
  df_potential_csc[marker, cells_for_this_marker] <- TRUE

  # Add to seurat object
  sample_seurat$cells_for_csc_marker <- FALSE 
  sample_seurat@meta.data[cells_for_this_marker,"cells_for_csc_marker"] <- TRUE
  
  ### Compare cells for this marker to samples
  tt <- as.data.frame(t(as.matrix(table(sample_seurat$cells_for_csc_marker, sample_seurat$orig.ident)) *100) /  as.numeric(table(sample_seurat$orig.ident))) # we divide by the total cells per sample to have percent
  colnames(tt) <- c("sample", "csc_for_this_marker","Freq")
  
  print(ggplot(tt, aes(x=sample, y=Freq, fill=csc_for_this_marker))+
          geom_bar(stat="identity")+
          scale_fill_manual(values=c("grey", "orange"))+
          coord_cartesian(ylim = c(0, 10))+ # zoom to 10% of frequency and not 100% because otherwise it is very hard to see the potential csc on the bar charts
          ylab("Freq (%)")+
          theme_classic()+
          theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
          ggtitle(paste0("Distribution of cells with > 5 UMI counts for ", marker)))
  
  ### Compare cells for this marker to subtypes
  tt2 <- as.data.frame(t(as.matrix(table(sample_seurat$cells_for_csc_marker, sample_seurat$H_LP_M_clean_groups_1)) *100) /  as.numeric(table(sample_seurat$H_LP_M_clean_groups_1))) # we divide by the total cells per sample to have percent
  colnames(tt2) <- c("subtype", "csc_for_this_marker","Freq")
  
  print(ggplot(tt2, aes(x=subtype, y=Freq, fill=csc_for_this_marker))+
          geom_bar(stat="identity")+
          scale_fill_manual(values=c("grey", "orange"))+
          coord_cartesian(ylim = c(0, 10))+ # zoom to 10% of frequency and not 100% because otherwise it is very hard to see the potential csc on the bar charts
          ylab("Freq (%)")+
          theme_classic()+
          ggtitle(paste0("Distribution of cells with > 5 UMI counts for ", marker)))
}
dev.off()

# remove markers that are not sufficiently detected 
sum_potential_csc_per_marker <- apply(df_potential_csc, 1, sum)
markers_to_keep <- names(sum_potential_csc_per_marker)[which(sum_potential_csc_per_marker!=0)]
df_potential_csc <- df_potential_csc[markers_to_keep,] # 12 markers

# compute nb of markers highly expressed for each cell

nb_markers_per_cell <- apply(df_potential_csc, 2, sum)
df_potential_csc <- rbind(df_potential_csc, nb_markers_per_cell); rownames(df_potential_csc)[nrow(df_potential_csc)] <- "nb_markers_per_cell"

save(df_potential_csc, file=paste0(resdir, "Liver_CSC_markers_potential_csc_min_", umi_threshold, "_UMI_SCT.RData"))


##### 2. Define possible groups of potential csc #####
#-----------------------------------------------------

nb_common_cells_between_markers <- as.data.frame(matrix(NA, nrow=nrow(df_potential_csc)-1, ncol=nrow(df_potential_csc)-1)) # remove 1 row because corresponds to the row "nb markers per cell" in df_potential_csc
rownames(nb_common_cells_between_markers) <- colnames(nb_common_cells_between_markers) <- rownames(df_potential_csc)[1:(nrow(df_potential_csc)-1)]
for(i in 1:nrow(nb_common_cells_between_markers)){
  for(j in 1:ncol(nb_common_cells_between_markers)){
    nb_common_cells_between_markers[i,j] <- length(which(df_potential_csc[i,]==1 & df_potential_csc[j,]==1))
  }
}

save(nb_common_cells_between_markers, file=paste0(resdir, "Liver_CSC_markers_potential_csc_shared_markers_min_", umi_threshold, "_UMI_SCT.RData"))


###################################
##### 3. POTENTIAL CSC GLOBAL #####
###################################

##### 1. Identify potential CSC #####
#------------------------------------

### keep cells which are potential csc for at least one of the markers csc_markers_liver_main = PROM1, EPCAM, THY1, CD24, CD44

df_potential_csc_main <- df_potential_csc[csc_markers_liver_main,]

potential_csc_global <- apply(df_potential_csc_main, 2, sum) # if the cells are potential csc for at least one marker the sum will be at least 1
potential_csc_global. <- names(potential_csc_global)[which(potential_csc_global!=0)]; length(potential_csc_global.) # 427 cells, remove cells that are csc for no marker among PROM1, EPCAM, THY1, CD24, CD44

save(potential_csc_global., file=file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/potential_cancer_stem_cells.RData"))
df_potential_csc_main <- df_potential_csc_main[,potential_csc_global.]
save(df_potential_csc_main, file=paste0(resdir, "Liver_main_CSC_markers_potential_csc_min_", umi_threshold, "_UMI_SCT.RData"))

### identify to which marker the potential csc is associated

potential_csc_global_associated_markers <- rep(NA, times=length(potential_csc_global.)); names(potential_csc_global_associated_markers) <- potential_csc_global. # create vector that will contain the info
for(i in 1:length(potential_csc_global_associated_markers)){
  potential_csc_global_associated_markers[i] <- rownames(df_potential_csc_main)[as.logical(df_potential_csc_main[,i])] # for each cell take the name of the corresponding marker (works only when 1 marker exactly)
} 

# for cells with 2 markers:
potential_csc_global_associated_markers[names(potential_csc_global)[which(potential_csc_global==2)]] <- "2 markers"


### plot on UMAP/PCA

sample_seurat$potential_csc <- "not CSC" # create column in seurat object
sample_seurat@meta.data[names(potential_csc_global_associated_markers),"potential_csc"] <- paste("CSC,", potential_csc_global_associated_markers, sep=" ")

pdf(file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/UMAP_potential_cancer_stem_cells.pdf"))
DimPlot(sample_seurat, group.by = "potential_csc", order=c( unique(sample_seurat$potential_csc)[-1], "not CSC"), cols = c("grey", "red", "orange", "yellow", "green", "grey20", "blue"), reduction="umap", pt.size=0.9)
DimPlot(sample_seurat, group.by = "potential_csc", order=c( unique(sample_seurat$potential_csc)[-1], "not CSC"), cols = c("grey", "red", "orange", "yellow", "green", "grey20", "blue"), reduction="pca", pt.size=0.9)
dev.off()


##### 2. Compare potential CSC to annotations #####
#--------------------------------------------------

##### 1. Samples 

tt1 <- as.data.frame(table(sample_seurat$potential_csc, sample_seurat$orig.ident))
tt1$CSC <- ifelse(tt1$Var1=="not CSC", "not CSC", "CSC") # to separe globally not CSC from CSC for barplot

pdf(file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/potential_cancer_stem_cells_vs_samples.pdf"))
ggplot(tt1, aes(x=CSC, y=Freq, fill=Var2))+
  geom_bar(stat="identity", position="fill")+
  theme_classic()
ggplot(tt1, aes(x=Var1, y=Freq, fill=Var2))+
  geom_bar(stat="identity", position="fill")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(tt1, aes(x=Var2, y=Freq, fill=Var1))+
  geom_bar(stat="identity", position="fill")+
  theme_classic()+
  scale_fill_manual(values=c("yellow", "orange", "green", "blue","grey20", "red", "grey"))+
  coord_cartesian(ylim = c(1,0.9)) # zoom to 10% of frequency and not 100% because otherwise it is very hard to see the potential csc on the bar charts
  dev.off()

save(tt1, file=file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/potential_cancer_stem_cells_vs_samples.RData"))


##### 2. Subtypes

tt2 <- as.data.frame(table(sample_seurat$potential_csc, sample_seurat$H_LP_M_clean_groups_1))
tt2$CSC <- ifelse(tt2$Var1=="not CSC", "not CSC", "CSC") # to separe globally not CSC from CSC for barplot

pdf(file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/potential_cancer_stem_cells_vs_subtypes.pdf"))
ggplot(tt2, aes(x=CSC, y=Freq, fill=Var2))+
  geom_bar(stat="identity", position="fill")+
  theme_classic()+
  scale_fill_manual(values=c("#E8AAD6", "#BC53CD", "#510787", "#6F3323"))
ggplot(tt2, aes(x=Var1, y=Freq, fill=Var2))+
  geom_bar(stat="identity", position="fill")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c("#E8AAD6", "#BC53CD", "#510787", "#6F3323"))
ggplot(tt2, aes(x=Var2, y=Freq, fill=Var1))+
  geom_bar(stat="identity", position="fill")+
  theme_classic()+
  scale_fill_manual(values=c("yellow", "orange", "green", "blue","grey20", "red", "grey"))+
  coord_cartesian(ylim = c(1,0.9)) # zoom to 10% of frequency and not 100% because otherwise it is very hard to see the potential csc on the bar charts
dev.off()

save(tt2, file=file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/potential_cancer_stem_cells_vs_samples.RData"))


##### 3. CNV clusters

pdf(file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/potential_cancer_stem_cells_vs_cnv_clusters_per_sample.pdf"))
for(name_sample in unique(sample_seurat$orig.ident)){ # the CNV clusters are made per sample and not on the merged dataset
  cnv_clusters <- geco.load(file.path(cnv_dir, name_sample, "Final_clusters/final_CNV_clusters.RData"))
  sample_seurat$cnv_clusters_per_sample <- NA
  sample_seurat@meta.data[cnv_clusters$cells,"cnv_clusters_per_sample"] <- cnv_clusters$final_CNV_cluster # add CNV cluster info to the corresponding cells
  
  tt3 <- as.data.frame(table(sample_seurat$potential_csc, sample_seurat$cnv_clusters_per_sample))
  tt3$CSC <- ifelse(tt3$Var1=="not CSC", "not CSC", "CSC") # to separe globally not CSC from CSC for barplot
  print(ggplot(tt3, aes(x=Var2, y=Freq, fill=Var1))+
    geom_bar(stat="identity", position="fill")+
    theme_classic()+
    scale_fill_manual(values=c("yellow", "orange", "green", "blue","grey20", "red", "grey"))+
    ggtitle(name_sample)+
    coord_cartesian(ylim = c(1,0.9)))
  
  save(tt3, file=paste0(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/potential_cancer_stem_cells_vs_cnv_clusters_", name_sample, ".RData"))
}
dev.off()


##### 3. Fisher test between CSC marker pairs #####
#--------------------------------------------------

# create table that will contain the fisher test p values 
fisher_test_marker_pairs <- as.data.frame(matrix(NA, nrow=length(csc_markers_liver_main), ncol=length(csc_markers_liver_main)))
rownames(fisher_test_marker_pairs) <- colnames(fisher_test_marker_pairs) <-csc_markers_liver_main

df_potential_csc. <- df_potential_csc[csc_markers_liver_main,] # restrict to the main csc markers but keep all tumor cells = different from df_potential_csc_main which has only potential CSC

for(i in 1:(nrow(df_potential_csc.)-1)){ # we will compare each marker pair wise based on all tumor cells
  for(j in (i+1):nrow(df_potential_csc.)){
    tt <- table(df_potential_csc.[i,], df_potential_csc.[j,])
    ft <- fisher.test(tt)
    fisher_test_marker_pairs[i,j] <- ft$p.value  
  }
}

write.table(fisher_test_marker_pairs, file=file.path(resdir,"Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/potential_cancer_stem_cells_fisher_test_markers.txt"), sep="\t", row.names=T, col.names=T)







#-------------------------------------------------------------------------------------------------------
### OLD ###

##### 6. Enrichment of csc markers in the subtypes #####
#-------------------------------------------------------

# create categories of subtypes LP vs all, M vs all

sample_seurat$LP_vs_all <- ifelse(sample_seurat$H_LP_M_clean_groups_1 == "LP", "LP", "other")
sample_seurat$M_vs_all <- ifelse(sample_seurat$H_LP_M_clean_groups_1 == "M", "M", "other")

### CD44 in M tumor cells
sample_seurat$csc_marker <- exp_csc["CD44",colnames(sample_seurat)]
print(ggplot(sample_seurat@meta.data, aes(x=M_vs_all, y=csc_marker))+
  geom_boxplot()+
  theme_classic()+
  ggtitle("CD44 expression")+
  stat_compare_means())
ggsave(file.path(resdir, "CD44_exp_M_vs_other_RNA_data.pdf"))


### other markers in LP tumor cells

for(g in csc_markers_liver){
  sample_seurat$csc_marker <- exp_csc[g,colnames(sample_seurat)]
  print(ggplot(sample_seurat@meta.data, aes(x=LP_vs_all, y=csc_marker))+
          geom_boxplot()+
          theme_classic()+
          ggtitle(paste0(g, " expression"))+
          stat_compare_means())
  ggsave(paste0(resdir, g, "_exp_LP_vs_other_RNA_data.pdf"))
  
}



