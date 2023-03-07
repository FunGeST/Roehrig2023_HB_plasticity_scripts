
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(circlize)
library(ggpubr)

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
resdir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/Gene_signatures"


###########################
##### 0. DATA LOADING #####
###########################

sample_seurat <- geco.load(file.path(datadir, "Seurat_object_analyses.RData"))

### Load Song signatures
song_sig <- read.table(file.path(resdir, "Song_2022_Nature_Comm.csv"), sep=";", header=T)
song_sig_HB1 <- unique(song_sig$Hepatoblast_I); song_sig_HB1 <- intersect(song_sig_HB1, rownames(sample_seurat))
song_sig_HB2 <- unique(song_sig$Hepatoblast_II); song_sig_HB2 <- intersect(song_sig_HB2, rownames(sample_seurat))
song_sig_erythroid <- unique(song_sig$Erythroid_like); song_sig_erythroid <- intersect(song_sig_erythroid, rownames(sample_seurat))
song_sig_fibroblast <- unique(song_sig$Fibroblast_like); song_sig_fibroblast <- intersect(song_sig_fibroblast, rownames(sample_seurat))
song_sig_neuroendocrine <- unique(song_sig$Neuroendocrine); song_sig_neuroendocrine <- intersect(song_sig_neuroendocrine, rownames(sample_seurat))

### Huang signatures
huang_sig_HB1 <- c("AFP", "GPC3", "MZT2B")
huang_sig_HB1_stemness <- c("AFP", "DUSP9", "GPC3", "DLK1")
huang_sig_HB1_prolif <- c("CDK1", "TOP2A", "CCNB1", "CCNB2")
huang_sig_HB1_cholangio <- c("KRT7", "KRT8", "KRT19", "ANXA4")
huang_sig_HB1_hepato <- c("CYP2E1", "CYP2C9", "ADH1B", "GLUL")
huang_sig_HB2 <- c("VIM", "COL1A1", "COL1A2", "WIF1")
huang_sig_HB3 <- c("ALB", "HP" ,"ADH1A", "ADH1B", "UGT2B10")
huang_sig_HB4 <- c()

### Create PCA dataframe to project
df_pca <- as.data.frame(sample_seurat@reductions$pca@cell.embeddings[,c("PC_1", "PC_2")]) # extract PC cell coordinates
if(identical(rownames(df_pca), rownames(sample_seurat@meta.data))){ # add cell annotations
  df_pca$cell <- rownames(sample_seurat@meta.data)
  df_pca$sample <- factor(sample_seurat@meta.data$orig.ident, levels=unique(sample_seurat@meta.data$orig.ident))
  df_pca$subtype <- sample_seurat@meta.data$H_LP_M_clean_groups_1}

### Compute mean expression for each gene set
df_pca$mean_song_HB1 <- apply(sample_seurat@assays$RNA@data[intersect(song_sig_HB1, rownames(sample_seurat)),],2,mean)
df_pca$mean_song_HB2 <- apply(sample_seurat@assays$RNA@data[intersect(song_sig_HB2, rownames(sample_seurat)),],2,mean)
df_pca$mean_song_erythroid <- apply(sample_seurat@assays$RNA@data[intersect(song_sig_erythroid, rownames(sample_seurat)),],2,mean)
df_pca$mean_song_fibroblast <- apply(sample_seurat@assays$RNA@data[intersect(song_sig_fibroblast, rownames(sample_seurat)),],2,mean)
df_pca$mean_song_neuroendocrine <- apply(sample_seurat@assays$RNA@data[intersect(song_sig_neuroendocrine, rownames(sample_seurat)),],2,mean)

df_pca$mean_huang_HB1 <- apply(sample_seurat@assays$RNA@data[intersect(huang_sig_HB1, rownames(sample_seurat)),],2,mean)
df_pca$mean_huang_HB1_stemness <- apply(sample_seurat@assays$RNA@data[intersect(huang_sig_HB1_stemness, rownames(sample_seurat)),],2,mean)
df_pca$mean_huang_HB1_prolif <- apply(sample_seurat@assays$RNA@data[intersect(huang_sig_HB1_prolif, rownames(sample_seurat)),],2,mean)
df_pca$mean_huang_HB1_cholangio <- apply(sample_seurat@assays$RNA@data[intersect(huang_sig_HB1_cholangio, rownames(sample_seurat)),],2,mean)
df_pca$mean_huang_HB1_hepato <- apply(sample_seurat@assays$RNA@data[intersect(huang_sig_HB1_hepato, rownames(sample_seurat)),],2,mean)
df_pca$mean_huang_HB2 <- apply(sample_seurat@assays$RNA@data[intersect(huang_sig_HB2, rownames(sample_seurat)),],2,mean)
df_pca$mean_huang_HB3 <- apply(sample_seurat@assays$RNA@data[intersect(huang_sig_HB3, rownames(sample_seurat)),],2,mean)

### shuffle cells for visualization
df_pca <- df_pca[sample(rownames(df_pca)),] # shuffle randomly the cell position in the table to have random overlap of the points in the grahp


################################
##### 1. PROJECTION ON PCA #####
################################

pdf(file.path(resdir, "Song_signatures_PCA.pdf"))
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_song_HB1))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Song Hepatoblast I")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_song_HB2))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Song Hepatoblast II")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_song_erythroid))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Song Erythroid-like")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_song_fibroblast))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Song Fibroblast-like")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_song_neuroendocrine))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Song Neuroendocrine")
dev.off()

pdf(file.path(resdir, "Huang_signatures_PCA.pdf"))
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_huang_HB1))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Huang HB1")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_huang_HB1_prolif))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Huang HB1 proliferation")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_huang_HB1_stemness))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Huang HB1 stemness")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_huang_HB1_cholangio))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Huang HB1 cholangiocyte")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_huang_HB1_hepato))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Huang HB1 hepatocytes")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_huang_HB2))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Huang HB2")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_huang_HB3))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Huang HB3")
dev.off()


#######################
##### 2. BOXPLOTS #####
#######################

df_pca$subtype <- factor(df_pca$subtype, levels=c("M", "LP", "H+LP", "H"))
pdf(file.path(resdir, "Other_HB_scRNAseq_studies_signatures_boxplots.pdf"))
### Song
ggplot(df_pca, aes(x=subtype, y=mean_song_HB1, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Song Hepatoblast I")
ggplot(df_pca, aes(x=subtype, y=mean_song_HB2, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Song Hepatoblast II")
ggplot(df_pca, aes(x=subtype, y=mean_song_erythroid, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Song Erythroid-like")
ggplot(df_pca, aes(x=subtype, y=mean_song_fibroblast, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Song Fibroblast-like")
ggplot(df_pca, aes(x=subtype, y=mean_song_neuroendocrine, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Song Neuroendocrine")

### Huang
ggplot(df_pca, aes(x=subtype, y=mean_huang_HB1, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Huang HB1")
ggplot(df_pca, aes(x=subtype, y=mean_huang_HB1_stemness, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Huang HB1 stemness")
ggplot(df_pca, aes(x=subtype, y=mean_huang_HB1_prolif, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Huang HB1 proliferation")
ggplot(df_pca, aes(x=subtype, y=mean_huang_HB1_cholangio, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Huang HB1 cholangiocytes")
ggplot(df_pca, aes(x=subtype, y=mean_huang_HB1_hepato, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Huang HB1 hepaocytes")
ggplot(df_pca, aes(x=subtype, y=mean_huang_HB2, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Huang HB2")
ggplot(df_pca, aes(x=subtype, y=mean_huang_HB3, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Huang HB3")
dev.off()
