
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

datadir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/4.Pseudo_CNV/1.Results/All_samples_merged/"
resdir <- file.path(datadir, "UMAP_visualization"); if(!dir.exists(resdir)){dir.create(resdir)}
seurat_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/No_regression/Merged/All_samples"

preliminary_results <- c("yes", "no")[1] # check "yes" if the inferCNV algo has not finished yet and only the preliminary results are available (but it should not change the groups obtained, only the CNV values that should be clearer after filtering)
if(preliminary_results=="yes"){
  preliminary_suffix <- "preliminary."
}else{preliminary_suffix <- ""}

samples <- c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")

j=5 # nb de clusters


###########################
##### 0. DATA LOADING #####
###########################

dendogram <- read.newick(file=paste0(datadir, "/infercnv.", preliminary_suffix, "observations_dendrogram.txt"))
object <- readRDS(paste0(datadir, preliminary_suffix, "infercnv_obj"))
object_gene_order <- object@gene_order # extract the gene coordinates
object_matrix <- object@expr.data # extract the values displayed on the inferCNV heatmap 

sample <- geco.load(paste0(seurat_dir, "/Seurat_object_analyses.RData"))

ref_cells <- geco.load("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/4.Pseudo_CNV/0.Input_data/normal_cells_for_CNV_CHC2959N_clusters_0_2.RData"); ref_cells <- gsub("CHC2959N_", "", ref_cells)
ref_cells <- paste0("CHC2959N_", ref_cells)


#######################################
##### 1. CELL ANNOTATION FOR UMAP #####
#######################################

dendogram_sample <- as.dendrogram(dendogram)
clusters_infercnv <- cutree(dendogram_sample, k=j)

#save(clusters_infercnv, file=file.path(resdir, "inferCNV_clusters.RData"))

### Display on UMAP ###

sample$infercnv_clusters <- NA

sample@meta.data[names(clusters_infercnv), "infercnv_clusters"] <- clusters_infercnv
sample@meta.data[ref_cells, "infercnv_clusters"] <- "ref"
sample$tumor_normal <- NA
sample@meta.data[which(sample@meta.data$infercnv_clusters%in%c(1,2,3,5)),"tumor_normal"] <- "tumor"
sample@meta.data[which(sample@meta.data$infercnv_clusters%in%c(4,"ref")),"tumor_normal"] <- "normal/immune"

pdf(paste0(resdir, "/UMAP_tumor_vs_normal_immune_inferCNV_no_cluster_by_sample.pdf"))
DimPlot(sample, group.by="infercnv_clusters")+ggtitle(paste0("inferCNV clusters (k=", j, ")"))
DimPlot(sample, group.by="tumor_normal")
dev.off()

tumor_cells_inferCNV <- rownames(sample@meta.data[which(sample@meta.data$tumor_normal=="tumor"),])
normal_cells_inferCNV <- rownames(sample@meta.data[which(sample@meta.data$tumor_normal=="normal/immune"),])

save(tumor_cells_inferCNV, file=file.path(resdir, "tumor_cells_inferCNV.RData"))
save(normal_cells_inferCNV, file=file.path(resdir, "normal_cells_inferCNV.RData"))


########################################################
##### 2. IMPROVE INFERCNV CLUSTERING VISUALIZATION #####
########################################################

##### 1. Reorder inferCNV clusters for nicer visualization #####
#---------------------------------------------------------------
clusters_current.order <- cutree(dendogram_sample, k=7)
clusters_new.order <- c(clusters_current.order[c(which(clusters_current.order==1), which(clusters_current.order==2), which(clusters_current.order==3), which(clusters_current.order==5), which(clusters_current.order==6), which(clusters_current.order==7), which(clusters_current.order==4))])
clusters.reordered <- dendextend::rotate(as.dendrogram(dendogram_sample), order = names(clusters_new.order))
clusters_new.order <- cutree(clusters.reordered, k=7)


##### 2. Project inferCNV clustering with complexheatmap #####
#-------------------------------------------------------------
### Create gene annotation (chromosome) ###
if(identical(rownames(object_gene_order), rownames(object_matrix))){
  row_ha =rowAnnotation(chromosomes = object_gene_order$chr)
} # check that the genes in matrix are ordered by position

### Create cell annotation (sample of origin) 
col_ha = columnAnnotation(sample = sample@meta.data[names(clusters_new.order),"orig.ident"],
                          col = list(sample = c("CHC2959N"="grey20", "CHC3377N"="grey50", "CHC2959T"="#F8766D",
                                                "CHC2960T"="#B79F00", "CHC3133T"="#00BA38", "CHC3377T"="#00BFC4",
                                                "CHC3610T"="#619CFF", "CHC3662T"="#F564E3"))) 

### Print heatmap
# takes ~10 min

pdf(paste0(resdir, "/inferCNV_clustering_reordered.pdf"))
ht <- draw(Heatmap(as.matrix(object_matrix)[,names(clusters_new.order)], show_column_names = F,show_row_names=F, cluster_rows=F,  column_dend_height = unit(4, "cm"),
                   cluster_columns = clusters.reordered, column_title = "InferCNV clustering reordered",  use_raster = T, left_annotation = row_ha,
                   top_annotation = col_ha, col=circlize::colorRamp2(c(0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3), c("#00008E", "#4848AB", "#B6B6DD", "white", "#DDB6B6", "#AB4848", "#630000"))))
dev.off()

#dendro <- column_dend(ht) # save the whole dendrogram
#save(dendro, file=paste0(resdir, "/dendrogram_", clustering_distance, "_", clustering_method, ".RData"))



#############################################################
##### 3. BARPLOT OF SAMPLE DISTRIBUTION PER CNV CLUSTER #####
#############################################################

k = 7 # define nb of inferCNV clusters to look at
clusters_infercnv2 <- cutree(dendogram, k=k)
# Rename clusters for the paper
clusters_infercnv2[which(clusters_infercnv2==1)] <- "Alt_6"; clusters_infercnv2[which(clusters_infercnv2==2)] <- "Alt_5"; clusters_infercnv2[which(clusters_infercnv2==3)] <- "Alt_4"
clusters_infercnv2[which(clusters_infercnv2==4)] <- "No_alt"
clusters_infercnv2[which(clusters_infercnv2==5)] <- "Alt_3"; clusters_infercnv2[which(clusters_infercnv2==6)] <- "Alt_2"; clusters_infercnv2[which(clusters_infercnv2==7)] <- "Alt_1"

# Add to Seurat object
sample$infercnv_clusters2 <- NA
sample@meta.data[names(clusters_infercnv2), "infercnv_clusters2"] <- clusters_infercnv2
sample@meta.data[ref_cells, "infercnv_clusters2"] <- "ref"

tt <- table(sample$orig.ident, sample$infercnv_clusters2) # confront inferCNV clusters to samples
tt <- reshape2::melt(tt); tt <- tt %>% filter(Var2 != "ref") # remove reference cells because of no use here

ggplot(tt, aes(x=Var2, y=value,  fill=Var1))+
  geom_bar(position="fill", stat="identity")+
  theme_classic()
ggsave(file.path(resdir, "sample_distribution_per_infercnv_cluster.png"))


