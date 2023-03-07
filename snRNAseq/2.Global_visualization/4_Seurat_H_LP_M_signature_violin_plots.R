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

datadir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/No_regression/Merged_tumors"
resdir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/No_regression/Comparison_signatures"

keep_MT = "With_MT_genes"

bulk_dir = "F:/Amelie-Datas/Pediatric tumors/2.Expression/RNAseq/Without new NL samples/Supervised"


###########################
##### 0. DATA LOADING #####
###########################

sample <- geco.load(file.path(datadir, "Seurat_object_analyses.RData"))
DefaultAssay(sample) <- "RNA"
sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor=10000) 


##### Gene signatures #####
#--------------------------

diff_H <- geco.load(paste0(bulk_dir, "/cluster_F_vs_HB/top_resullts_limma.RData"))
diff_LP <- geco.load(paste0(bulk_dir, "/cluster_E_vs_HB/top_resullts_limma.RData"))
diff_M <- geco.load(paste0(bulk_dir, "/cluster_M_vs_HB/top_resullts_limma.RData"))

LP_markers <- rownames(diff_LP)[which(diff_LP$logFC>3)] #133 genes
H_markers <- rownames(diff_H)[which(diff_H$logFC>3)] # 71 genes
M_markers <- rownames(diff_M)[which(diff_M$logFC>3)] # 304 genes
LP_markers <- LP_markers[which(LP_markers%in%rownames(sample))]
H_markers <- H_markers[which(H_markers%in%rownames(sample))]
M_markers <- M_markers[which(M_markers%in%rownames(sample))]


###########################################################
##### 1. VIOLIN PLOTS OF H/LP/M SIGNATURES PER SAMPLE #####
###########################################################

mean_LP <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[LP_markers,], 2, mean)
mean_H <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[H_markers,], 2, mean)
mean_M <- apply(GetAssayData(object=sample, slot="data", assay="RNA")[M_markers,], 2, mean)

exp_table <- factoall(data.frame(mean_H=mean_H, mean_LP=mean_LP, mean_M=mean_M))
rownames(exp_table) <- colnames(sample)
exp_table$sample <- sample$orig.ident 

pdf(file.path(resdir, "violin_plots_H_signature.pdf"))
ggplot(exp_table, aes(x=sample, y=mean_H, fill=sample)) + 
  geom_violin()+
  theme(legend.position="none")
ggplot(exp_table, aes(x=sample, y=mean_H, fill=sample)) + 
  geom_violin()+
  geom_jitter(size=0.5, position=position_jitter(0.2), col="black")+
  theme(legend.position="none")
dev.off()

pdf(file.path(resdir, "violin_plots_LP_signature.pdf"))
ggplot(exp_table, aes(x=sample, y=mean_LP, fill=sample)) + 
  geom_violin()+
  theme(legend.position="none")
ggplot(exp_table, aes(x=sample, y=mean_LP, fill=sample)) + 
  geom_violin()+
  geom_jitter(size=0.5, position=position_jitter(0.2), col="black")+
  theme(legend.position="none")
dev.off()

pdf(file.path(resdir, "violin_plots_M_signature.pdf"))
ggplot(exp_table, aes(x=sample, y=mean_M, fill=sample)) + 
  geom_violin()+
  theme(legend.position="none")
ggplot(exp_table, aes(x=sample, y=mean_M, fill=sample)) + 
  geom_violin()+
  geom_jitter(size=0.5, position=position_jitter(0.2), col="black")+
  theme(legend.position="none")
dev.off()


  