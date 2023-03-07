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

resdir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict"; if(!dir.exists(resdir)){dir.create(resdir)}

### Define 2 sets of thresholds to identify H and LP subtypes on PC2
t1_1 = -15; t2_1 = 15 
t1_2 = -20; t2_2 = 0; t3_2 = 17; t4_2 = 35


###########################
##### 0. DATA LOADING ##### 
###########################

sample <- geco.load(file.path(resdir, "Seurat_object_analyses.RData")) # object with all tumor cells


###########################################################################
##### 1. IDENTIFICATION OF H/LP/M GLOBAL CLUSTERS: CLASSIFICATION N°1 #####
###########################################################################

##### Define visual thresholds on PCA based on signatures #####
#--------------------------------------------------------------
t1 = t1_1; t2 = t2_1 

pdf(file.path(resdir, "H_LP_M_definition_thresholds_1.pdf"))
FeaturePlot(sample, features="mean_H", dims=c(1,2), reduction="pca")+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  geom_hline(yintercept = t1, linetype="dashed")+
  geom_hline(yintercept = t2, linetype="dashed")
FeaturePlot(sample, features="mean_LP", dims=c(1,2), reduction="pca")+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  geom_hline(yintercept = t1, linetype="dashed")+
  geom_hline(yintercept = t2, linetype="dashed")
DimPlot(sample, group.by = "orig.ident", reduction="pca", dims=c(1,2))+
  geom_hline(yintercept = t1, linetype="dashed")+
  geom_hline(yintercept = t2, linetype="dashed")
dev.off()


##### Define clusters #####
#--------------------------
# separate the mesenchymal CHC3610T for easier identification

sample$H_LP_M_clean_groups_1 <- NA
sample@meta.data[which(sample@meta.data$orig.ident!="CHC3610T" 
                       & sample@reductions$pca@cell.embeddings[,"PC_2"] < t1),"H_LP_M_clean_groups_1"] <- "H"
sample@meta.data[which(sample@meta.data$orig.ident!="CHC3610T"
                       & sample@reductions$pca@cell.embeddings[,"PC_2"] > t1 
                       & sample@reductions$pca@cell.embeddings[,"PC_2"] < t2),"H_LP_M_clean_groups_1"] <- "H+LP"
sample@meta.data[which(sample@meta.data$orig.ident!="CHC3610T"
                       & sample@reductions$pca@cell.embeddings[,"PC_2"] > t2),"H_LP_M_clean_groups_1"] <- "LP"
sample@meta.data[which(sample@meta.data$orig.ident=="CHC3610T"),"H_LP_M_clean_groups_1"] <- "M"
sample$H_LP_M_clean_groups_1 <- factor(sample$H_LP_M_clean_groups_1, levels=c("H", "H+LP", "LP", "M"))

H_LP_M_per_sample <- as.data.frame(table(sample$orig.ident, sample$H_LP_M_clean_groups_1))
colnames(H_LP_M_per_sample) <- c("sample", "H_LP_M_group", "Freq")

pdf(file.path(resdir, "H_LP_M_clean_groups_1.pdf"))
DimPlot(sample, group.by = "H_LP_M_clean_groups_1", reduction="pca", dims=c(1,2), cols=c("#E8AAD6", "#BC53CD", "#510787", "#6F3323"))
DimPlot(sample, group.by = "H_LP_M_clean_groups_1", reduction="umap", cols=c("#E8AAD6", "#BC53CD", "#510787", "#6F3323"))
dev.off()

pdf(file.path(resdir, "H_LP_M_clean_groups_1_vs_samples.pdf"))
ggplot(data=H_LP_M_per_sample, aes(x=sample, y=Freq, fill=H_LP_M_group ))+
  geom_bar(stat="identity", position="fill") +
  xlab("")+
  ylab("%")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c("#E8AAD6","#CD82DA", "#510787", "#6F3323"))
ggplot(data=H_LP_M_per_sample, aes(x=H_LP_M_group, y=Freq, fill=sample ))+
  geom_bar(stat="identity", position="fill") +
  xlab("")+
  ylab("%")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


###########################################################################
##### 2. IDENTIFICATION OF H/LP/M GLOBAL CLUSTERS: CLASSIFICATION N°2 #####
###########################################################################

##### Define visual thresholds on PCA based on signatures #####
#--------------------------------------------------------------

t1 = t1_2; t2 = t2_2; t3 = t3_2; t4 = t4_2

pdf(file.path(resdir, "H_LP_M_definition_thresholds_2.pdf"))
FeaturePlot(sample, features="mean_H", dims=c(1,2), reduction="pca")+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  geom_hline(yintercept = t1, linetype="dashed")+
  geom_hline(yintercept = t2, linetype="dashed")+
  geom_hline(yintercept = t3, linetype="dashed")+
  geom_hline(yintercept = t4, linetype="dashed")
FeaturePlot(sample, features="mean_LP", dims=c(1,2), reduction="pca")+
  scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))+
  geom_hline(yintercept = t1, linetype="dashed")+
  geom_hline(yintercept = t2, linetype="dashed")+
  geom_hline(yintercept = t3, linetype="dashed")+
  geom_hline(yintercept = t4, linetype="dashed")
DimPlot(sample, group.by = "orig.ident", reduction="pca", dims=c(1,2))+
  geom_hline(yintercept = t1, linetype="dashed")+
  geom_hline(yintercept = t2, linetype="dashed")+
  geom_hline(yintercept = t3, linetype="dashed")+
  geom_hline(yintercept = t4, linetype="dashed")
dev.off()


##### Define clusters #####
#--------------------------

# separate the mesenchymal CHC3610T for easier identification

sample$H_LP_M_clean_groups_2 <- NA
sample@meta.data[which(sample@meta.data$orig.ident!="CHC3610T" 
                       & sample@reductions$pca@cell.embeddings[,"PC_2"] < t1),"H_LP_M_clean_groups_2"] <- "H ++"
sample@meta.data[which(sample@meta.data$orig.ident!="CHC3610T"
                       & sample@reductions$pca@cell.embeddings[,"PC_2"] > t1 
                       & sample@reductions$pca@cell.embeddings[,"PC_2"] < t2),"H_LP_M_clean_groups_2"] <- "H +"
sample@meta.data[which(sample@meta.data$orig.ident!="CHC3610T"
                       & sample@reductions$pca@cell.embeddings[,"PC_2"] > t2 
                       & sample@reductions$pca@cell.embeddings[,"PC_2"] < t3),"H_LP_M_clean_groups_2"] <- "H+LP low"
sample@meta.data[which(sample@meta.data$orig.ident!="CHC3610T"
                       & sample@reductions$pca@cell.embeddings[,"PC_2"] > t3
                       & sample@reductions$pca@cell.embeddings[,"PC_2"] < t4),"H_LP_M_clean_groups_2"] <- "LP +"
sample@meta.data[which(sample@meta.data$orig.ident!="CHC3610T"
                       & sample@reductions$pca@cell.embeddings[,"PC_2"] > t4),"H_LP_M_clean_groups_2"] <- "LP ++"
sample@meta.data[which(sample@meta.data$orig.ident=="CHC3610T"),"H_LP_M_clean_groups_2"] <- "M"
sample$H_LP_M_clean_groups_2 <- factor(sample$H_LP_M_clean_groups_2, levels=c("H ++", "H +", "H+LP low", "LP +", "LP ++", "M"))

H_LP_M_per_sample <- as.data.frame(table(sample$orig.ident, sample$H_LP_M_clean_groups_2))
colnames(H_LP_M_per_sample) <- c("sample", "H_LP_M_group", "Freq")

pdf(file.path(resdir, "H_LP_M_clean_groups_2.pdf"))
DimPlot(sample, group.by = "H_LP_M_clean_groups_2", reduction="pca", dims=c(1,2), cols=c("#E8AAD6", "#E088C3", "#BC53CD", "#8C3FC1", "#510787", "#6F3323"))
DimPlot(sample, group.by = "H_LP_M_clean_groups_2", reduction="umap", cols=c("#E8AAD6", "#E088C3", "#BC53CD", "#8C3FC1", "#510787", "#6F3323"))
dev.off()

pdf(file.path(resdir, "H_LP_M_clean_groups_2_vs_samples.pdf"))
ggplot(data=H_LP_M_per_sample, aes(x=sample, y=Freq, fill=H_LP_M_group ))+
  geom_bar(stat="identity", position="fill") +
  xlab("")+
  ylab("%")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c("#E8AAD6", "#E088C3", "#CD82DA", "#8C3FC1", "#510787", "#6F3323"))
ggplot(data=H_LP_M_per_sample, aes(x=H_LP_M_group, y=Freq, fill=sample ))+
  geom_bar(stat="identity", position="fill") +
  xlab("")+
  ylab("%")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


#####################################
##### 3. SAVE CELLS PER SUBTYPE #####
#####################################

dir.create(file.path(resdir, "Cells_per_subtype"))

### Clean groups 1

cells_H <- rownames(sample@meta.data)[which(sample@meta.data$H_LP_M_clean_groups_1=="H")]
cells_H_LP <- rownames(sample@meta.data)[which(sample@meta.data$H_LP_M_clean_groups_1=="H+LP")]
cells_LP <- rownames(sample@meta.data)[which(sample@meta.data$H_LP_M_clean_groups_1=="LP")]
cells_M_no_GLI3 <- rownames(sample@meta.data)[which(sample@meta.data$H_LP_M_clean_groups_1=="M"&sample@meta.data$SCT_snn_res.0.3!=15)]
cells_GLI3 <- rownames(sample@meta.data)[which(sample@meta.data$H_LP_M_clean_groups_1=="M"&sample@meta.data$SCT_snn_res.0.3==15)]

save(cells_H, file=file.path(resdir, "Cells_per_subtype/cells_H.RData"))
save(cells_H_LP, file=file.path(resdir, "Cells_per_subtype/cells_H_LP.RData"))
save(cells_LP, file=file.path(resdir, "Cells_per_subtype/cells_LP.RData"))
save(cells_M_no_GLI3, file=file.path(resdir, "Cells_per_subtype/cells_M_no_GLI3.RData"))
save(cells_GLI3, file=file.path(resdir, "Cells_per_subtype/cells_GLI3.RData"))

### Clean groups 2

cells_H_very_high <- rownames(sample@meta.data)[which(sample@meta.data$H_LP_M_clean_groups_2=="H ++")]
cells_H_high <- rownames(sample@meta.data)[which(sample@meta.data$H_LP_M_clean_groups_2=="H +")]
cells_H_LP_low <- rownames(sample@meta.data)[which(sample@meta.data$H_LP_M_clean_groups_2=="H+LP low")]
cells_LP_high <- rownames(sample@meta.data)[which(sample@meta.data$H_LP_M_clean_groups_2=="LP +")]
cells_LP_very_high <- rownames(sample@meta.data)[which(sample@meta.data$H_LP_M_clean_groups_2=="LP ++")]
cells_M_no_GLI3 <- rownames(sample@meta.data)[which(sample@meta.data$H_LP_M_clean_groups_2=="M"&sample@meta.data$SCT_snn_res.0.3!=15)]
cells_GLI3 <- rownames(sample@meta.data)[which(sample@meta.data$H_LP_M_clean_groups_2=="M"&sample@meta.data$SCT_snn_res.0.3==15)]

save(cells_H_very_high, file=file.path(resdir, "Cells_per_subtype/cells_H_very_high.RData"))
save(cells_H_high, file=file.path(resdir, "Cells_per_subtype/cells_H_high.RData"))
save(cells_H_LP_low, file=file.path(resdir, "Cells_per_subtype/cells_H_LP_low.RData"))
save(cells_LP_high, file=file.path(resdir, "Cells_per_subtype/cells_LP_high.RData"))
save(cells_LP_very_high, file=file.path(resdir, "Cells_per_subtype/cells_LP_very_high.RData"))


###################
##### 3. SAVE #####
###################

save(sample, file=file.path(resdir,"Seurat_object_analyses.RData"))




