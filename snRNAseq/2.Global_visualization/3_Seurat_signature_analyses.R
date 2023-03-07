
library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)

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

sample_name="CHC3133T"
keep_MT_genes=c("yes", "no")[1]
regression=c("No_regression", "MT_regression", "Cell_cycle_regression")[1]
norm=c("NormalizeData", "SCTransform")[2]

if(keep_MT_genes=="yes"){genes_MT_or_not <- "With_MT_genes"}else{genes_MT_or_not <- "Without_MT_genes"}

datadir <- paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/", norm, "/", regression, "/", sample_name, "/", genes_MT_or_not)
resdir <- file.path(datadir, "Genes_or_signatures_on_UMAP"); if(!dir.exists(resdir)){dir.create(resdir)}


###########################
##### 0. DATA LOADING #####
###########################

sample <- geco.load(file.path(datadir, "Seurat_object_analyses.RData"))
DefaultAssay(sample) <- "RNA"
sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor=10000) # il faut faire la moyenne sur les données NormalizeData et pas SCT (meilleure pratique selon Satijalab)


##### External files #####
#-------------------------

gs.db <- readLines("F:/Amelie-Datas/Data-hepatopartage/msigdb.v6.1.symbols.gmt")

diff_markers <- unlist(strsplit(gs.db[[grep("HSIAO_LIVER_SPECIFIC_GENES", gs.db)]], "\t"))[3:length(unlist(strsplit(gs.db[[grep("HSIAO_LIVER_SPECIFIC_GENES", gs.db)]], "\t")))]
inflam_markers <- unlist(strsplit(gs.db[[grep("HALLMARK_INFLAMMATORY_RESPONSE", gs.db)]], "\t"))[3:length(unlist(strsplit(gs.db[[grep("HALLMARK_INFLAMMATORY_RESPONSE", gs.db)]], "\t")))]
prolif_markers <- unlist(strsplit(gs.db[[grep("REACTOME_CELL_CYCLE", gs.db)[1]]], "\t"))[3:length(unlist(strsplit(gs.db[[grep("REACTOME_CELL_CYCLE", gs.db)[1]]], "\t")))]
diff_markers <- diff_markers[which(diff_markers%in%rownames(sample))]
inflam_markers <- inflam_markers[which(inflam_markers%in%rownames(sample))]
prolif_markers <- prolif_markers[which(prolif_markers%in%rownames(sample))]

diff_LP <- geco.load("F:/Amelie-Datas/Pediatric tumors/2.Expression/RNAseq/Without new NL samples/Supervised/cluster_E_vs_HB/top_resullts_limma.RData")
diff_H <- geco.load("F:/Amelie-Datas/Pediatric tumors/2.Expression/RNAseq/Without new NL samples/Supervised/cluster_F_vs_HB/top_resullts_limma.RData")
diff_M <- geco.load("F:/Amelie-Datas/Pediatric tumors/2.Expression/RNAseq/Without new NL samples/Supervised/cluster_M_vs_HB/top_resullts_limma.RData")
LP_markers <- rownames(diff_LP)[which(diff_LP$logFC>3)]
H_markers <- rownames(diff_H)[which(diff_H$logFC>3)]
M_markers <- rownames(diff_M)[which(diff_M$logFC>3)]
LP_markers <- LP_markers[which(LP_markers%in%rownames(sample))]
H_markers <- H_markers[which(H_markers%in%rownames(sample))]
M_markers <- M_markers[which(M_markers%in%rownames(sample))]

hepatic_markers <- c("ALB", "BAAT", "ADH1A", "FGA"); hepatic_markers <- hepatic_markers[which(hepatic_markers%in%rownames(sample))]
progenitor_markers <- c("MYCN", "KRT19", "AFP", "PROM1"); progenitor_markers <- progenitor_markers[which(progenitor_markers%in%rownames(sample))]
mesenchymal_markers <- c("VIM", "THY1"); mesenchymal_markers <- mesenchymal_markers[which(mesenchymal_markers%in%rownames(sample))]
immune_markers <- c("CD3D", "CD8A", "CD160", "CD19", "HLA-DOA", "CD1A", "CD80", "CXCR1", "CDH5"); immune_markers <- immune_markers[which(immune_markers%in%rownames(sample))]
proliferation_markers <- c("CCNA2", "CDC20", "BUB1", "AURKA"); proliferation_markers <- proliferation_markers[which(proliferation_markers%in%rownames(sample))]

module_H <- c("OLIG3", "ISX", "TBX10", "MLXIPL", "NR1I2", "CREB3L3", "CEBPA", "HNF1A", "HNF4A"); module_H <- module_H[which(module_H%in%rownames(sample))]
module_LP <- c("LIN28B", "MIXL1", "LHX1", "SALL3", "DMRT2", "EVX1"); module_LP <- module_LP[which(module_LP%in%rownames(sample))]
module_inflam <- c("ESR1", "TCF21", "IKZF3", "SATB1", "ZMAT1", "NR3C2"); module_inflam <- module_inflam[which(module_inflam%in%rownames(sample))]
module_M <- c("TBX5", "ZIC1", "ZIC4", "SOX11", "MSX2", "DLX3", "WT1", "HAND2", "TWIST1", "FOXF2", "TRPS1", "GLI3"); module_M <- module_M[which(module_M%in%rownames(sample))]

endoderm <- c("APELA", "TRH", "LIN28A", "POU5F1", "CER1", "CLDN7", "PODXL", "DNMT3B", "KRT8", "CLDN4", "CDH6", "PYY", "AMOT", "KRT18", "SOX4"); endoderm <- endoderm[which(endoderm%in%rownames(sample))]; length(endoderm)
migrating_hepatoblasts <- c("PROX1", "NR5A2", "MKI67", "MGA", "NEDD4", "DSP", "RELN", "VCAN", "ROBO2", "FLRT2", "MYH10", "KMT2A", "SETD2", "TFRC", "FZD3", "MET", "ONECUT1", "ONECUT2", "TET1", "HHEX", "TET2", "FBXO30", "EPCAM", "TBX3", "TPM4", "HMGA2", "MEIS1", "SOX11", "LMNB1")
migrating_hepatoblasts <- migrating_hepatoblasts[which(migrating_hepatoblasts%in%rownames(sample))]; length(migrating_hepatoblasts)
hepatoblasts <- c("TRF", "APOA1"," APOM", "ALB", "APOA2", "F2", "TST", "APOE", "CITED1", "FGA", "DLK1", "F10", "VCAM1", "FST", "TTR", "ICAM1", "FOXA3"); hepatoblasts <- hepatoblasts[which(hepatoblasts%in%rownames(sample))]; length(hepatoblasts)
hepatomesenchyme <- c("YAF2", "AKIRIN2", "JUND", "PTOV1", "COL3A1", "HOXB4", "MFAP4", "HOXA5", "HOXB2", "BCL11A", "TMEM11", "HOXC4", "PTN", "CDH11", "MMP2", "FOXF1", "COL1A1", "NR2F1", "TCF21", "TWIST2", "MFAP2", "IGF1")
hepatomesenchyme <- hepatomesenchyme[which(hepatomesenchyme%in%rownames(sample))]; length(hepatomesenchyme)

genes_dvp <- read.table("F:/Amelie-Datas/Data-hepatopartage/Genes_early_late_development_cairo.csv", sep=";", header=T)
genes_early <- genes_dvp[which(genes_dvp$Early.Late.Ratio>1),"Human.Gene.Identifier"]; genes_early <- genes_early[which(genes_early%in%rownames(sample))]
genes_late <- genes_dvp[which(genes_dvp$Early.Late.Ratio<1),"Human.Gene.Identifier"]; genes_late <- genes_late[which(genes_late%in%rownames(sample))]


#########################
##### 1. SIGNATURES #####
#########################

##### 1. Inflammation/Differentiation/Proliferation signatures #####
#|------------------------------------------------------------------

pdf(file.path(resdir, "Signature_diff_inflam_prolif.pdf"), width=10, height=10)
mean_diff <- apply(GetAssayData(object=sample, slot="data")[diff_markers,], 2, mean); sample$mean_diff <- mean_diff
mean_inflam <- apply(GetAssayData(object=sample, slot="data")[inflam_markers,], 2, mean); sample$mean_inflam <- mean_inflam
mean_prolif <- apply(GetAssayData(object=sample, slot="data")[prolif_markers,], 2, mean); sample$mean_prolif <- mean_prolif
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
p1 <- FeaturePlot(sample, features = c("mean_diff", "mean_prolif", "mean_inflam"), order=T, combine = FALSE)
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
# pour avoir même échelle
scale_max <- max(c(max(mean_diff), max(mean_inflam), max(mean_prolif)))
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'), limits=c(0,scale_max))
p1 <- FeaturePlot(sample, features = c("mean_diff", "mean_prolif", "mean_inflam"), order=T, combine = FALSE)
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
dev.off()


##### 2. Tumor signatures #####
#------------------------------

pdf(file.path(resdir, "Signature_H_LP_M_bulk_diff_exp.pdf"), width=10, height=10)
mean_LP <- apply(GetAssayData(object=sample, slot="data")[LP_markers,], 2, mean); sample$mean_LP <- mean_LP
mean_H <- apply(GetAssayData(object=sample, slot="data")[H_markers,], 2, mean); sample$mean_H <- mean_H
mean_M <- apply(GetAssayData(object=sample, slot="data")[M_markers,], 2, mean); sample$mean_M <- mean_M
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
p1 <- FeaturePlot(sample, features = c("mean_H", "mean_LP", "mean_M"),order=T, combine = FALSE ) 
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
# pour avoir même échelle
scale_max <- max(c(max(mean_LP), max(mean_H), max(mean_M)))
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'), limits=c(0,scale_max))
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
dev.off()


##### 3. Bulk TF module signatures #####
#---------------------------------------

pdf(file.path(resdir, "Signature_bulk_TF_modules.pdf"), width=10, height=10)
mean_module_LP <- apply(GetAssayData(object=sample, slot="data")[module_LP,], 2, mean); sample$mean_module_LP <- mean_module_LP
mean_module_H <- apply(GetAssayData(object=sample, slot="data")[module_H,], 2, mean); sample$mean_module_H <- mean_module_H
mean_module_M <- apply(GetAssayData(object=sample, slot="data")[module_M,], 2, mean); sample$mean_module_M <- mean_module_M
mean_module_inflam <- apply(GetAssayData(object=sample, slot="data")[module_inflam,], 2, mean); sample$mean_module_inflam <- mean_module_inflam
if(length(module_LP)!=0&length(module_H)!=0&length(module_M)!=0&length(module_inflam)!=0){
  fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
  p1 <- FeaturePlot(sample, features = c("mean_module_H", "mean_module_LP", "mean_module_M", "mean_module_inflam"),order=T, combine = FALSE ) 
  p2 <- lapply(p1, function (x) x + fix.sc_affy)
  CombinePlots(p2)
  # pour avoir même échelle
  scale_max <- max(c(max(mean_module_LP), max(mean_module_H), max(mean_module_M), max(mean_module_inflam)))
  fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'), limits=c(0,scale_max))
  p1 <- FeaturePlot(sample, features = c("mean_module_H", "mean_module_LP", "mean_module_M", "mean_module_inflam"),order=T, combine = FALSE ) 
  p2 <- lapply(p1, function (x) x + fix.sc_affy)
  CombinePlots(p2)
}else if(length(module_LP)==0){
  fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
  p1 <- FeaturePlot(sample, features = c("mean_module_H", "mean_module_M", "mean_module_inflam"),order=T, combine = FALSE ) 
  p2 <- lapply(p1, function (x) x + fix.sc_affy)
  CombinePlots(p2)
  # pour avoir même échelle
  scale_max <- max(c(max(mean_module_H), max(mean_module_M), max(mean_module_inflam)))
  fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'), limits=c(0,scale_max))
  p1 <- FeaturePlot(sample, features = c("mean_module_H", "mean_module_M", "mean_module_inflam"),order=T, combine = FALSE ) 
  p2 <- lapply(p1, function (x) x + fix.sc_affy)
  CombinePlots(p2)
}else if(length(module_M)==0){
  fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
  p1 <- FeaturePlot(sample, features = c("mean_module_H", "mean_module_LP", "mean_module_inflam"),order=T, combine = FALSE ) 
  p2 <- lapply(p1, function (x) x + fix.sc_affy)
  CombinePlots(p2)
  # pour avoir même échelle
  scale_max <- max(c(max(mean_module_H), max(mean_module_LP), max(mean_module_inflam)))
  fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'), limits=c(0,scale_max))
  p1 <- FeaturePlot(sample, features = c("mean_module_H", "mean_module_LP", "mean_module_inflam"),order=T, combine = FALSE ) 
  p2 <- lapply(p1, function (x) x + fix.sc_affy)
  CombinePlots(p2)
}else if(length(module_inflam)==0){
  fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
  p1 <- FeaturePlot(sample, features = c("mean_module_H", "mean_module_LP", "mean_module_M"),order=T, combine = FALSE ) 
  p2 <- lapply(p1, function (x) x + fix.sc_affy)
  CombinePlots(p2)
  # pour avoir même échelle
  scale_max <- max(c(max(mean_module_H), max(mean_module_LP), max(mean_module_M)))
  fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'), limits=c(0,scale_max))
  p1 <- FeaturePlot(sample, features = c("mean_module_H", "mean_module_LP", "mean_module_M"),order=T, combine = FALSE ) 
  p2 <- lapply(p1, function (x) x + fix.sc_affy)
  CombinePlots(p2)
}else if(length(module_H)==0){
  fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
  p1 <- FeaturePlot(sample, features = c("mean_module_LP", "mean_module_M", "mean_module_inflam"),order=T, combine = FALSE ) 
  p2 <- lapply(p1, function (x) x + fix.sc_affy)
  CombinePlots(p2)
  # pour avoir même échelle
  scale_max <- max(c(max(mean_module_LP), max(mean_module_M), max(mean_module_inflam)))
  fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'), limits=c(0,scale_max))
  p1 <- FeaturePlot(sample, features = c("mean_module_LP", "mean_module_M", "mean_module_inflam"),order=T, combine = FALSE ) 
  p2 <- lapply(p1, function (x) x + fix.sc_affy)
  CombinePlots(p2)
}
dev.off()


##### 4. Lotto developmental signatures #####
#--------------------------------------------

pdf(file.path(resdir, "Signature_Lotto_development.pdf"), width=10, height=10)
mean_endoderm <- apply(GetAssayData(object=sample, slot="data")[endoderm,], 2, mean); sample$mean_endoderm <- mean_endoderm
mean_migrating_hepatoblasts <- apply(GetAssayData(object=sample, slot="data")[migrating_hepatoblasts,], 2, mean); sample$mean_migrating_hepatoblasts <- mean_migrating_hepatoblasts
mean_hepatoblasts <- apply(GetAssayData(object=sample, slot="data")[hepatoblasts,], 2, mean); sample$mean_hepatoblasts <- mean_hepatoblasts
mean_hepatomesenchyme <- apply(GetAssayData(object=sample, slot="data")[hepatomesenchyme,], 2, mean); sample$mean_hepatomesenchyme <- mean_hepatomesenchyme
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
p1 <- FeaturePlot(sample, features = c("mean_endoderm", "mean_migrating_hepatoblasts", "mean_hepatoblasts", "mean_hepatomesenchyme"),order=T, combine = FALSE ) 
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
# pour avoir même échelle
scale_max <- max(c(max(mean_endoderm), max(mean_migrating_hepatoblasts), max(mean_hepatoblasts), max(mean_hepatomesenchyme)))
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'), limits=c(0,scale_max))
p1 <- FeaturePlot(sample, features = c("mean_endoderm", "mean_migrating_hepatoblasts", "mean_hepatoblasts", "mean_hepatomesenchyme"),order=T, combine = FALSE ) 
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
dev.off()


##### 5. Cairo early/late signatures #####
#-----------------------------------------

pdf(file.path(resdir, "Signature_Cairo_early_late.pdf"), width=10, height=5)
mean_early <- apply(GetAssayData(object=sample, slot="data")[genes_early,], 2, mean); sample$mean_early <- mean_early
mean_late <- apply(GetAssayData(object=sample, slot="data")[genes_late,], 2, mean); sample$mean_late <- mean_late
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
p1 <- FeaturePlot(sample, features = c("mean_early", "mean_late"),order=T, combine = FALSE ) 
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
# pour avoir même échelle
scale_max <- max(c(max(mean_early), max(mean_late)))
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'), limits=c(0,scale_max))
p1 <- FeaturePlot(sample, features = c("mean_early", "mean_late"),order=T, combine = FALSE ) 
p2 <- lapply(p1, function (x) x + fix.sc_affy)
CombinePlots(p2)
dev.off()


mean_signatures <- sample@meta.data[,c("mean_diff", "mean_prolif", "mean_inflam", "mean_H", "mean_LP", "mean_M")]
save(mean_signatures, file=file.path(resdir, "Signatures_mean_NormData.RData"))


#################################
##### 2. INDIVIDUAL MARKERS #####
#################################

##### 1. Typical markers #####
#-----------------------------

pdf(file.path(resdir, "Markers_hepatic.pdf"), width=8, height=8)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in hepatic_markers){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_progenitor.pdf"), width=8, height=8)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in progenitor_markers){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_mesenchymal.pdf"), width=8, height=8)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in mesenchymal_markers){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_immune.pdf"), width=8, height=8)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in immune_markers){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_proliferation.pdf"), width=8, height=8)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in proliferation_markers){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()


##### 2. Immune markers #####
#----------------------------

pdf(file.path(resdir, "Markers_B_cells.pdf"), width=6, height=6)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in c("HLA-DRA", "IGHM")){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_cholangiocytes.pdf"), width=6, height=6)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in c("CD24", "SPP1")){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_endothelial.pdf"), width=6, height=6)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in c("CDH5", "PECAM1", "ADAMTS9")){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_epithelial.pdf"), width=6, height=6)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in c("FLT1")){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_fibroblasts.pdf"), width=6, height=6)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in c("ZEB2", "DCN")){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_hepatic_stellate.pdf"), width=6, height=6)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in c("COL6A3", "COL1A1", "COL3A1")){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_hepatocytes.pdf"), width=6, height=6)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in c("HPGD", "CYP3A4", "APOC3", "HP")){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_kupffer.pdf"), width=6, height=6)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in c("CD163", "CD68", "VCAM1")){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_macrophages.pdf"), width=6, height=6)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in c("MRC1", "CD163", "MS4A4A")){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_T_cells.pdf"), width=6, height=6)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in c("CD247", "TEHMIS", "CD96")){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()
pdf(file.path(resdir, "Markers_T_NK_cells.pdf"), width=6, height=6)
fix.sc_affy <- scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20'))
for(gene in c("CD3D", "CD3E")){
  p1 <- FeaturePlot(sample, features=gene, order=T); p2 <- p1 + fix.sc_affy
  print(p2)
}
dev.off()



