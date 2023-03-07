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

sample_names <- c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")

datadir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/No_regression/Merged/All_samples"
infercnv_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/4.Pseudo_CNV/1.Results/All_samples_merged/UMAP_visualization"
resdir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune"; if(!dir.exists(resdir)){dir.create(resdir)}
atac_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/3.snATACseq/2.Basic_analyses/2.Individual_samples"


###########################
##### 0. DATA LOADING ##### 
###########################

sample_all <- geco.load(file.path(datadir, "Seurat_object_analyses.RData"))
tumor_cells_infercnv <- geco.load(file.path(infercnv_dir, "tumor_cells_inferCNV.RData"))
normal_cells_infercnv <- geco.load(file.path(infercnv_dir, "normal_cells_inferCNV.RData"))

gs.db <- readLines("F:/Amelie-Datas/Data-hepatopartage/msigdb.v6.1.symbols.gmt")

inflam_markers <- unlist(strsplit(gs.db[[grep("HALLMARK_INFLAMMATORY_RESPONSE", gs.db)]], "\t"))[3:length(unlist(strsplit(gs.db[[grep("HALLMARK_INFLAMMATORY_RESPONSE", gs.db)]], "\t")))]
inflam_markers <- inflam_markers[which(inflam_markers%in%rownames(sample_all))]


############################################################
##### 1. DISPLAY IMPORTANT FEATURES ON UMAP ALL CELLS  #####
############################################################

##### 1. InferCNV results #####
#------------------------------

sample_all$tumor_normal_infercnv <- NA
sample_all@meta.data[which(rownames(sample_all@meta.data)%in%tumor_cells_infercnv),"tumor_normal_infercnv"] <- "tumor"
sample_all@meta.data[which(rownames(sample_all@meta.data)%in%normal_cells_infercnv),"tumor_normal_infercnv"] <- "normal/immune"

pdf(file.path(resdir, "UMAP_all_cells_tumor_normal_infercnv.pdf"))
DimPlot(sample_all, group.by = "tumor_normal_infercnv")
dev.off()


##### 2. Normal/tumor clusters with Seurat (marker visualization) #####
#----------------------------------------------------------------------

### Identification of immune clusters thanks to PCA: PC2 and PC3

pca_coord <- as.data.frame(sample_all@reductions$pca@cell.embeddings)

sample_all$mean_inflam <- apply(GetAssayData(object=sample_all, slot="data", assay="RNA")[inflam_markers,], 2, mean)
FeaturePlot(sample_all, features="mean_inflam", reduction="pca", dims=c(2,3)) + scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))

FeaturePlot(sample_all, features="CD247", reduction="pca", dims=c(2,3))
FeaturePlot(sample_all, features="CDH5", reduction="pca", dims=c(2,3))
FeaturePlot(sample_all, features="CD163", reduction="pca", dims=c(2,3))
FeaturePlot(sample_all, features="DCN", reduction="pca", dims=c(2,3))
FeaturePlot(sample_all, features="COL3A1", reduction="pca", dims=c(2,3))

pdf(file.path(resdir, "Identification_immune_cells_PCA_2_3.pdf"))
DimPlot(sample_all, reduction="pca", dims=c(2,3), group.by = "orig.ident")+
  geom_vline(xintercept = -2, linetype="dashed")+
  geom_hline(yintercept = 25, linetype="dashed")

immune_cells_pc2_3 <- rownames(pca_coord)[which(pca_coord$PC_2< -2 & pca_coord$PC_3<25)]
sample_all$immune_cells_pca <- NA
sample_all@meta.data[which(rownames(sample_all@meta.data)%in%immune_cells_pc2_3),"immune_cells_pca"] <- "immune"

DimPlot(sample_all, reduction="pca", dims=c(2,3), group.by = "immune_cells_pca") 
DimPlot(sample_all, group.by = "immune_cells_pca")
dev.off()

sample_all <- FindClusters(sample_all, resolution =4) # take low resolution to have few clusters to better identify tumor/non tumor
DimPlot(sample_all, group.by = "SCT_snn_res.4", label=T) + NoLegend()
#sample_all <- FindClusters(sample_all, resolution =5) # take low resolution to have few clusters to better identify tumor/non tumor

sample_all$tumor_normal_seurat <- "tumor"
sample_all@meta.data[which(sample_all@meta.data$SCT_snn_res.4%in%c(6,10,15,21,22,24,28,29,33:35,37,42:45,52,53,55,56,58:66,69)),"tumor_normal_seurat"] <- "normal/immune"
# keep GLI3 cluster for the moment (suspicious tumor cells)

pdf(file.path(resdir, "UMAP_all_cells_tumor_normal_seurat.pdf"))
DimPlot(sample_all, group.by = "tumor_normal_seurat")
dev.off()


##### 3. Doublets according to ATAC #####
#----------------------------------------

sample_all$doublet_enrich_ATAC <- NA 
sample_all$doublets_ATAC <- "no_doublet"

for(s in sample_names){
  atac_metadata <- as.data.frame(geco.load(file.path(atac_dir, s, "metadata.RData")))
  atac_metadata$cell_barcode <- gsub(paste0(s, "#"), paste0(s, "_"), gsub("-1", "", rownames(atac_metadata)))
  cells_in_RNA_not_ATAC <- rownames(sample_all@meta.data[which(sample_all@meta.data$orig.ident==s&!rownames(sample_all@meta.data)%in%atac_metadata$cell_barcode),])
  sample_all@meta.data[which(rownames(sample_all@meta.data)%in%cells_in_RNA_not_ATAC),"doublets_ATAC"] <- NA
  
  atac_metadata <- atac_metadata[which(atac_metadata$cell_barcode%in%colnames(sample_all)),] # common cells ATAC/RNA
  
  sample_all@meta.data[atac_metadata$cell_barcode,"doublet_enrich_ATAC"] <- atac_metadata$DoubletEnrichment 
  
  # choose doublet enrichment threshold to have ~15% of doublets (to be sure not to mistake clusters)
  if(s=="CHC2959N"){sample_all@meta.data[which(sample_all@meta.data$orig.ident==s&sample_all@meta.data$doublet_enrich_ATAC>15),"doublets_ATAC"] <- "doublet"}
  if(s=="CHC2959T"){sample_all@meta.data[which(sample_all@meta.data$orig.ident==s&sample_all@meta.data$doublet_enrich_ATAC>12),"doublets_ATAC"] <- "doublet"}
  if(s=="CHC2960T"){sample_all@meta.data[which(sample_all@meta.data$orig.ident==s&sample_all@meta.data$doublet_enrich_ATAC>4),"doublets_ATAC"] <- "doublet"}
  if(s=="CHC3133T"){sample_all@meta.data[which(sample_all@meta.data$orig.ident==s&sample_all@meta.data$doublet_enrich_ATAC>1),"doublets_ATAC"] <- "doublet"}
  if(s=="CHC3377N"){sample_all@meta.data[which(sample_all@meta.data$orig.ident==s&sample_all@meta.data$doublet_enrich_ATAC>1),"doublets_ATAC"] <- "doublet"}
  if(s=="CHC3377T"){sample_all@meta.data[which(sample_all@meta.data$orig.ident==s&sample_all@meta.data$doublet_enrich_ATAC>30),"doublets_ATAC"] <- "doublet"}
  if(s=="CHC3610T"){sample_all@meta.data[which(sample_all@meta.data$orig.ident==s&sample_all@meta.data$doublet_enrich_ATAC>6),"doublets_ATAC"] <- "doublet"}
  if(s=="CHC3662T"){sample_all@meta.data[which(sample_all@meta.data$orig.ident==s&sample_all@meta.data$doublet_enrich_ATAC>5),"doublets_ATAC"] <- "doublet"}
  
}

sample_all <- FindClusters(sample_all, resolution =5) # take very high resolution to try to idolate doublet clusters by expression
doublet_clusters <- as.data.frame(table(sample_all$doublets_ATAC, sample_all$SCT_snn_res.5))

pdf(file.path(resdir, "UMAP_all_cells_doublets_ATAC.pdf"))
DimPlot(sample_all, group.by = "doublets_ATAC")
DimPlot(sample_all, group.by = "SCT_snn_res.5") + NoLegend()
ggplot(data=doublet_clusters, aes(x=Var2, y=Freq, fill=Var1))+
  geom_bar(stat="identity", position="fill")+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  geom_hline(yintercept = 0.7, linetype="dashed")+
  ggtitle("ATAC doublet distribution in Seurat clusters (res=5")


### Identify possible and probable doublets clusters using high clustering resolution and ATAC doublet enrichment results
clusters_doublets_probable <- c(43,49,61,62,72,73,74) # at least 50% of doublets in cluster
clusters_doublets_possible <- c(13,64,65,67,70) # at least 30% of doublets in cluster

sample_all$doublet_clusters_potential <- "no doublet cluster"
sample_all@meta.data[which(sample_all@meta.data$SCT_snn_res.5%in%clusters_doublets_probable),"doublet_clusters_potential"] <- "probable doublet cluster"
sample_all@meta.data[which(sample_all@meta.data$SCT_snn_res.5%in%clusters_doublets_possible),"doublet_clusters_potential"] <- "possible doublet cluster"
DimPlot(sample_all, group.by = "doublet_clusters_potential")
dev.off()


#####################################################################
##### 2. SELECTION OF CELLS FOR IMMUNE/NORMAL AND TUMOR STUDIES #####
#####################################################################

##### 1. Identification of tumor cells with large thresholds = inferCNV tumor cells #####
#----------------------------------------------------------------------------------------

sample_all$selection_proposition_large <- sample_all$tumor_normal_infercnv

##### 2. Identification of tumor cells with strict thresholds = inferCNV tumor cells + visual non immune clusters and not doublets/etc on Seurat UMAP #####
#----------------------------------------------------------------------------------------------------------------------------------------------------------

sample_all$selection_proposition_strict <- "undefined"

sample_all@meta.data[which(sample_all@meta.data$tumor_normal_infercnv=="tumor" & 
                             sample_all@meta.data$tumor_normal_seurat=="tumor" &
                             sample_all@meta.data$doublet_clusters_potential != "probable doublet cluster" &
                             ! sample_all@meta.data$orig.ident %in% c("CHC2959N", "CHC3377N")),"selection_proposition_strict"] <- "tumor (clean)"
sample_all@meta.data[which(sample_all@meta.data$tumor_normal_infercnv=="normal/immune" &
                             sample_all@meta.data$immune_cells_pca=="immune"),"selection_proposition_strict"] <- "normal/immune (clean)"
sample_all@meta.data[which(sample_all@meta.data$tumor_normal_infercnv=="normal/immune" &
                             sample_all@meta.data$orig.ident%in%c("CHC2959N", "CHC3377N")),"selection_proposition_strict"] <- "normal/immune (clean)"


pdf(file.path(resdir, "UMAP_all_cells_selection_proposition.pdf"))
DimPlot(sample_all, group.by = "selection_proposition_large")
DimPlot(sample_all, group.by = "selection_proposition_strict")
dev.off()

sample_all$clusters_Seurat_res_5 <- sample_all$SCT_snn_res.5
annotation_all_cells <- sample_all@meta.data[,c("orig.ident", "tumor_normal_infercnv", "tumor_normal_seurat", "doublet_enrich_ATAC", "doublets_ATAC", "doublet_clusters_potential", "clusters_Seurat_res_5", "selection_proposition_large", "selection_proposition_strict")]

### Display tumor cells kept for each sample individually

pdf(file.path(resdir, "UMAP_tumor_cells_per_sample.pdf"))
for(s in c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3377T", "CHC3610T", "CHC3662T")){
  sample_all$tumor_cells_per_sample_large <- NA
  sample_all@meta.data[which(sample_all@meta.data$orig.ident==s&sample_all@meta.data$selection_proposition_large=="tumor"),"tumor_cells_per_sample_large"] <- "tumor cells"
  col_sample = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#609CFF", "#F564E3")[which(c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3377T", "CHC3610T", "CHC3662T")==s)]
  print(DimPlot(sample_all, group.by = "tumor_cells_per_sample_large", cols=col_sample))
  sample_all$tumor_cells_per_sample_strict <- NA
  sample_all@meta.data[which(sample_all@meta.data$orig.ident==s&sample_all@meta.data$selection_proposition_strict=="tumor (clean)"),"tumor_cells_per_sample_strict"] <- "tumor cells"
  col_sample = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#609CFF", "#F564E3")[which(c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3377T", "CHC3610T", "CHC3662T")==s)]
  print(DimPlot(sample_all, group.by = "tumor_cells_per_sample_strict", cols=col_sample))
}
dev.off()

### Display distribution of tumor/non-tumor cells in each sample
tt <- table(sample_all$orig.ident, sample_all$selection_proposition_strict)
tt <- reshape2::melt(tt); tt$Var2 <- factor(tt$Var2, levels=c("undefined", "tumor (clean)", "normal/immune (clean)")); tt$Var1 <- factor(tt$Var1, levels=c("CHC2959N", "CHC3377N", "CHC3377T", "CHC2960T", "CHC2959T", "CHC3662T", "CHC3133T", "CHC3610T"))

ggplot(tt, aes(x=Var1, y=value, fill=Var2))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=c("grey", "firebrick3", "forestgreen"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
ggsave(file.path(resdir, "nb_tumor_non_tumor_cells_by_sample.png"))


###################
##### 3. SAVE #####
###################

save(annotation_all_cells, file=file.path(resdir, "annotation_all_cells.RData"))
sample_all@meta.data$seurat_clusters <- sample_all@meta.data$tumor_normal <- NULL
save(sample_all, file=file.path(resdir, "Seurat_object_analyses_all_cells.RData"))


