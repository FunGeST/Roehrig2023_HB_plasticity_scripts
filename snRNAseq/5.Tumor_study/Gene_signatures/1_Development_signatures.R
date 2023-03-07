
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
resdir = file.path(datadir, "Gene_signatures"); if(!dir.exists(resdir)){dir.create(resdir)}

### Cairo genes (Supp Table 13)
cairo_genes <- read.table("F:/Amelie-Datas/Data-hepatopartage/Genes_early_late_development_cairo.csv", sep=";", header=T)
cairo_early <- unique(cairo_genes[which(cairo_genes$Early.Late.Ratio>1),"Human.Gene.Identifier"])
cairo_late <- unique(cairo_genes[which(cairo_genes$Early.Late.Ratio<1),"Human.Gene.Identifier"])

### Lotto genes (Fig S4a)
lotto_endoderm <- c("APELA", "TRH", "LIN28A", "POU5F1", "CER1", "CLDN7", "PODXL", "DNMT3B", "KRT8", "CLDN4", "CDH6", "PYY", "AMOT", "KRT18", "SOX4")
lotto_migrating_hepatoblasts <- c("PROX1", "NR5A2", "MKI67", "MGA", "NEDD4", "DSP", "RELN", "VCAN", "ROBO2", "FLRT2", "MYH10", "KMT2A", "SETD2", "TFRC", "FZD3", "MET", "ONECUT1", "ONECUT2", "TET1", "HHEX", "TET2", "FBXO30", "EPCAM", "TBX3", "TPM4", "HMGA2", "MEIS1", "SOX11", "LMNB1")
lotto_hepatoblasts <- c("TRF", "APOA1"," APOM", "ALB", "APOA2", "F2", "TST", "APOE", "CITED1", "FGA", "DLK1", "F10", "VCAM1", "FST", "TTR", "ICAM1", "FOXA3")
lotto_hepatomesenchyme <- c("YAF2", "AKIRIN2", "JUND", "PTOV1", "COL3A1", "HOXB4", "MFAP4", "HOXA5", "HOXB2", "BCL11A", "TMEM11", "HOXC4", "PTN", "CDH11", "MMP2", "FOXF1", "COL1A1", "NR2F1", "TCF21", "TWIST2", "MFAP2", "IGF1")

### Wesley genes (https://www.biorxiv.org/content/10.1101/2022.03.08.482299v1.full#F1, Fig1e)
wesley_hepatoblast1 <- c("SPINK1", "BRI3", "NPW", "FST", "MAP2K2")
wesley_hepatoblast2 <- c("MT1H", "IGF2", "AMN", "APOM", "DLK1", "RSPO3", "GSTP1", "BAAT")
wesley_fetal_hepatocyte1 <- c("AFP", "CYP3A7", "F2", "AHSG", "APOA1", "LIPC", "AGT")
wesley_fetal_hepatocyte2 <- c("ALB", "ANGPTL3", "CFHR1", "LEPR", "APOB", "GC", "VTN", "AKR1C1")
wesley_adult_hepatocyte <- c("CYP2E1", "SAA1", "HP", "C9", "ADH1B", "F9", "CES1", "TAT", "AZGP1", "APCS", "NNMT")


###########################
##### 0. DATA LOADING #####
###########################

sample_seurat <- geco.load(file.path(datadir, "Seurat_object_analyses.RData"))
pca_cells <- sample_seurat@reductions$pca@cell.embeddings # contains the PCA value per cell

### Create PCA dataframe to project
df_pca <- as.data.frame(sample_seurat@reductions$pca@cell.embeddings[,c("PC_1", "PC_2")]) # extract PC cell coordinates
if(identical(rownames(df_pca), rownames(sample_seurat@meta.data))){ # add cell annotations
  df_pca$cell <- rownames(sample_seurat@meta.data)
  df_pca$sample <- factor(sample_seurat@meta.data$orig.ident, levels=unique(sample_seurat@meta.data$orig.ident))
  df_pca$subtype <- sample_seurat@meta.data$H_LP_M_clean_groups_1}

### Compute mean expression for each gene set
df_pca$mean_cairo_early <- apply(sample_seurat@assays$RNA@data[intersect(cairo_early, rownames(sample_seurat)),],2,mean)
df_pca$mean_cairo_late <- apply(sample_seurat@assays$RNA@data[intersect(cairo_late, rownames(sample_seurat)),],2,mean)

df_pca$mean_lotto_endoderm <- apply(sample_seurat@assays$RNA@data[intersect(lotto_endoderm, rownames(sample_seurat)),],2,mean)
df_pca$mean_lotto_migrating_hepatoblasts <- apply(sample_seurat@assays$RNA@data[intersect(lotto_migrating_hepatoblasts, rownames(sample_seurat)),],2,mean)
df_pca$mean_lotto_hepatoblasts <- apply(sample_seurat@assays$RNA@data[intersect(lotto_hepatoblasts, rownames(sample_seurat)),],2,mean)
df_pca$mean_lotto_hepatomesenchyme <- apply(sample_seurat@assays$RNA@data[intersect(lotto_hepatomesenchyme, rownames(sample_seurat)),],2,mean)

df_pca$mean_wesley_hepatoblast1 <- apply(sample_seurat@assays$RNA@data[intersect(wesley_hepatoblast1, rownames(sample_seurat)),],2,mean)
df_pca$mean_wesley_hepatoblast2 <- apply(sample_seurat@assays$RNA@data[intersect(wesley_hepatoblast2, rownames(sample_seurat)),],2,mean)
df_pca$mean_wesley_fetal_hepatocyte1 <- apply(sample_seurat@assays$RNA@data[intersect(wesley_fetal_hepatocyte1, rownames(sample_seurat)),],2,mean)
df_pca$mean_wesley_fetal_hepatocyte2 <- apply(sample_seurat@assays$RNA@data[intersect(wesley_fetal_hepatocyte2, rownames(sample_seurat)),],2,mean)
df_pca$mean_wesley_adult_hepatocyte <- apply(sample_seurat@assays$RNA@data[intersect(wesley_adult_hepatocyte, rownames(sample_seurat)),],2,mean)

### shuffle cells for visualization
df_pca <- df_pca[sample(rownames(df_pca)),] # shuffle randomly the cell position in the table to have random overlap of the points in the grahp


################################
##### 1. PROJECTION ON PCA #####
################################

pdf(file.path(resdir, "Developmental_signatures_PCA.pdf"))
### Cairo
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_cairo_early))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Cairo early")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_cairo_late))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Cairo late")

### Lotto
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_lotto_endoderm))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Lotto endoderm")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_lotto_migrating_hepatoblasts))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Lotto migrating hepatoblasts")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_lotto_hepatoblasts))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Lotto hepatoblasts")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_lotto_hepatomesenchyme))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Lotto hepatomesenchyme")

### Wesley
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_wesley_hepatoblast1))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Wesley hepatoblast1")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_wesley_hepatoblast2))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Wesley hepatoblast2")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_wesley_fetal_hepatocyte1))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Wesley fetal hepatocyte1")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_wesley_fetal_hepatocyte2))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Wesley fetal hepatocyte2")
ggplot(df_pca, aes(x=PC_1, y=PC_2, col=mean_wesley_adult_hepatocyte))+
  geom_point(size=0.6)+
  theme_classic()+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  ggtitle("Wesley adult hepatocyte")
dev.off()


#######################
##### 2. BOXPLOTS #####
#######################

df_pca$subtype <- factor(df_pca$subtype, levels=c("M", "LP", "H+LP", "H"))
pdf(file.path(resdir, "Developmental_signatures_boxplots.pdf"))
### Cairo
ggplot(df_pca, aes(x=subtype, y=mean_cairo_early, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Cairo early")
ggplot(df_pca, aes(x=subtype, y=mean_cairo_late, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Cairo late")

### Lotto
ggplot(df_pca, aes(x=subtype, y=mean_lotto_endoderm, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Lotto endoderm")
ggplot(df_pca, aes(x=subtype, y=mean_lotto_migrating_hepatoblasts, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Lotto migrating hepatoblasts")
ggplot(df_pca, aes(x=subtype, y=mean_lotto_hepatoblasts, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Lotto hepatoblasts")
ggplot(df_pca, aes(x=subtype, y=mean_lotto_hepatomesenchyme, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Lotto hepatomesenchyme")

### Wesley
ggplot(df_pca, aes(x=subtype, y=mean_wesley_hepatoblast1, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Wesley hepatoblast1")
ggplot(df_pca, aes(x=subtype, y=mean_wesley_hepatoblast2, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Wesley hepatoblast2")
ggplot(df_pca, aes(x=subtype, y=mean_wesley_fetal_hepatocyte1, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Wesley fetal hepatocyte1")
ggplot(df_pca, aes(x=subtype, y=mean_wesley_fetal_hepatocyte2, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Wesley fetal hepatocyte2")
ggplot(df_pca, aes(x=subtype, y=mean_wesley_adult_hepatocyte, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  ggtitle("Wesley adult hepatocyte")
dev.off()
