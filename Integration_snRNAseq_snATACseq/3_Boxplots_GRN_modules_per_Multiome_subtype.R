
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


geco.enrichmentTest <- function (gene.sets, mylist, possibleIds, sep = ";", silent = F) 
{
  possibleIds <- unique(possibleIds)
  mylist <- unique(mylist)
  gene.sets <- lapply(gene.sets, unique)
  nids <- length(possibleIds)
  gene.sets <- lapply(gene.sets, function(x) intersect(x, possibleIds))
  nref <- sapply(gene.sets, length)
  if (all(nref == 0)) stop("Error: no intersection between gene sets and possible IDs.")
  if (any(nref == 0)) print("Warning: some of the gene sets have no intersection with possibleIds")
  if (!all(mylist %in% possibleIds)) stop("Error: some genes in mylist are not in possibleIds")
  if (!silent) cat(paste("NB : enrichment tests are based on", nids, "distinct ids.\n"))
  gene.sets <- gene.sets[nref > 0]
  n <- length(mylist)
  fun <- function(x) {
    y <- intersect(x, mylist)
    nx <- length(x)
    ny <- length(y)
    pval <- phyper(ny - 1, nx, nids - nx, n, lower.tail = F)
    c(nx, ny, pval,paste(y, collapse = sep))
  }
  tmp <- as.data.frame(t(sapply(gene.sets, fun)))
  rownames(tmp) <- names(gene.sets)
  for (i in 1:3) tmp[,i] <- as.numeric(as.character(tmp[,i]))
  tmp <- data.frame(tmp[,1:3],p.adjust(tmp[,3],method="BH"),tmp[,4])
  names(tmp) <- c("Nb_of_genes","Nb_of_deregulated_genes","p-value","q-value","Deregulated_genes")
  tmp
}


######################
##### PARAMETERS #####
######################

rna_dir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/"
resdir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/4.Integration/2.Gene_regulatory_network/Final_GRN/H_LP_M"


###########################
##### 0. DATA LOADING #####
###########################

### Load Multiome GRN
TF_ordered <- geco.load(file.path(resdir, "Heatmap_TF_targets_short_nb_targets_TF_ordered.RData"))
targets_ordered <- geco.load(file.path(resdir, "Heatmap_TF_targets_short_nb_targets_targets_ordered.RData"))

### Load Multiome GRN
list_GRN <- geco.load(file.path(resdir, "GRN_list.RData")) # list of targets per TF

### TF clusters
TF_M <- TF_ordered[which(TF_ordered=="TET1"):which(TF_ordered=="PRDM6")]
TF_LP <- TF_ordered[which(TF_ordered=="SOX4"):which(TF_ordered=="LHX1")]
TF_H_LP <- TF_ordered[which(TF_ordered=="MLXIPL"):which(TF_ordered=="MIXL1")]
TF_H <- TF_ordered[which(TF_ordered=="AR"):which(TF_ordered=="FOXP2")]

### Targets clusters
targets_M <- targets_ordered[which(targets_ordered=="ALKAL2"):which(targets_ordered=="GFRA1")]
targets_LP <- targets_ordered[which(targets_ordered=="SLC27A1"):which(targets_ordered=="IL17RB")]
targets_H_LP <- targets_ordered[which(targets_ordered=="PNPLA3"):which(targets_ordered=="AGT")]
targets_H <- targets_ordered[which(targets_ordered=="CD6"):which(targets_ordered=="RTN4RL1")]

targets_M_all <- unique(unlist(list_GRN[TF_M])); length(targets_M_all)
targets_LP_all <- unique(unlist(list_GRN[TF_LP])); length(targets_LP_all)
targets_H_LP_all <- unique(unlist(list_GRN[TF_H_LP])); length(targets_H_LP_all)
targets_H_all <- unique(unlist(list_GRN[TF_H])); length(targets_H_all)

### Expression tables
rna_sample_all <- geco.load(file.path(rna_dir, "Seurat_object_analyses.RData"))
rna_exp <- rna_sample_all@assays$RNA@data

# Load data for pathway enrichment
MSIG.ls <- geco.load("F:/Amelie-Datas/Data-hepatopartage/MSIG_v6.1.RData")


########################################################
##### 1. BOXPLOTS OF TF/TARGET MODULES IN SUBTYPES #####
########################################################

### Generate table of mean expression per cell
df <- data.frame(cell=colnames(rna_exp), subtype=NA, mean_TF_M=NA, mean_TF_LP=NA, mean_TF_H_LP=NA, mean_TF_H=NA, mean_targets_M=NA, mean_targets_LP=NA, mean_targets_H_LP=NA, mean_targets_H=NA, stringsAsFactors = F)
rownames(df) <- df$cell
df$subtype <- rna_sample_all@meta.data[df$cell, "H_LP_M_clean_groups_1"]
df$mean_TF_M <- apply(rna_exp[TF_M, df$cell], 2, mean)
df$mean_TF_LP <- apply(rna_exp[TF_LP, df$cell], 2, mean)
df$mean_TF_H_LP <- apply(rna_exp[TF_H_LP, df$cell], 2, mean)
df$mean_TF_H <- apply(rna_exp[TF_H, df$cell], 2, mean)
df$mean_targets_M <- apply(rna_exp[targets_M_all, df$cell], 2, mean)
df$mean_targets_LP <- apply(rna_exp[targets_LP_all, df$cell], 2, mean)
df$mean_targets_H_LP <- apply(rna_exp[targets_H_LP_all, df$cell], 2, mean)
df$mean_targets_H <- apply(rna_exp[targets_H_all, df$cell], 2, mean)

### Plot
df$subtype <- factor(df$subtype, levels=c("M", "LP", "H+LP", "H")) #order the subtypes
pdf(file.path(resdir, "Boxplots_TF_targets_modules_in_multiome_subtypes.pdf"), width=5, height=5)
ggplot(df, aes(x=subtype, y=mean_TF_M, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))
ggplot(df, aes(x=subtype, y=mean_TF_LP, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))
ggplot(df, aes(x=subtype, y=mean_TF_H_LP, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))
ggplot(df, aes(x=subtype, y=mean_TF_H, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))
ggplot(df, aes(x=subtype, y=mean_targets_M, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))
ggplot(df, aes(x=subtype, y=mean_targets_LP, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))
ggplot(df, aes(x=subtype, y=mean_targets_H_LP, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))
ggplot(df, aes(x=subtype, y=mean_targets_H, col=subtype, fill=subtype))+
  geom_boxplot(alpha=0.5, size=1)+
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")+
  scale_color_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))+
  scale_fill_manual(values=c("#6F3322", "#510787", "#BC53CD", "#E8AAD6"))
dev.off()


########################################
##### 2. HYPERGEOMETRIC ENRICHMENT #####
########################################

enrich_targets_H <- geco.enrichmentTest(gene.sets=MSIG.ls, mylist= targets_H_all, possibleIds= rownames(rna_sample_all))
enrich_targets_H_LP <- geco.enrichmentTest(gene.sets=MSIG.ls, mylist= targets_H_LP_all, possibleIds= rownames(rna_sample_all))
enrich_targets_LP <- geco.enrichmentTest(gene.sets=MSIG.ls, mylist= targets_LP_all, possibleIds= rownames(rna_sample_all))
enrich_targets_M <- geco.enrichmentTest(gene.sets=MSIG.ls, mylist= targets_M_all, possibleIds= rownames(rna_sample_all))

colnames(enrich_targets_H)[which(colnames(enrich_targets_H)=="q-value")] <- "qvalue"
colnames(enrich_targets_H_LP)[which(colnames(enrich_targets_H_LP)=="q-value")] <- "qvalue"
colnames(enrich_targets_LP)[which(colnames(enrich_targets_LP)=="q-value")] <- "qvalue"
colnames(enrich_targets_M)[which(colnames(enrich_targets_M)=="q-value")] <- "qvalue"

enrich_targets_H <- enrich_targets_H %>% arrange(qvalue) %>% filter(qvalue < 0.05)
enrich_targets_H_LP <- enrich_targets_H_LP %>% arrange(qvalue) %>% filter(qvalue < 0.05)
enrich_targets_LP <- enrich_targets_LP %>% arrange(qvalue) %>% filter(qvalue < 0.05)
enrich_targets_M <- enrich_targets_M %>% arrange(qvalue) %>% filter(qvalue < 0.05)

save(enrich_targets_H, file=file.path(resdir, "hypergeom_enrichment_H_targets.RData"))
save(enrich_targets_H_LP, file=file.path(resdir, "hypergeom_enrichment_H_LP_targets.RData"))
save(enrich_targets_LP, file=file.path(resdir, "hypergeom_enrichment_LP_targets.RData"))
save(enrich_targets_M, file=file.path(resdir, "hypergeom_enrichment_M_targets.RData"))

  