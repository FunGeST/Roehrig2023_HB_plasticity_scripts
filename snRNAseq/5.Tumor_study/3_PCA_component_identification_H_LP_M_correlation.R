

library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)


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

geco.changeRange <- function (v, newmin = 1, newmax = 10) 
{
  oldmin <- min(v, na.rm = TRUE)
  oldmax <- max(v, na.rm = TRUE)
  newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin))
}

######################
##### PARAMETERS #####
######################

datadir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/"
bulk_dir = "F:/Amelie-Datas/Pediatric tumors/2.Expression/RNAseq/Without new NL samples/Supervised"

resdir = file.path(datadir, "PCA_components_study")


###########################
##### 0. DATA LOADING #####
###########################

sample <- geco.load(file.path(datadir, "Seurat_object_analyses.RData"))

pca_coord <- sample@reductions$pca@cell.embeddings
pca_genes <- sample@reductions$pca@feature.loadings

save(pca_coord, file=file.path(resdir, "PCA_coordinates.RData"))
save(pca_genes, file=file.path(resdir, "PCA_genes.RData"))

### bulk H/LP/M markers

diff_H <- geco.load(paste0(bulk_dir, "/cluster_F_vs_HB/top_resullts_limma.RData"))
diff_LP <- geco.load(paste0(bulk_dir, "/cluster_E_vs_HB/top_resullts_limma.RData"))
diff_M <- geco.load(paste0(bulk_dir, "/cluster_M_vs_HB/top_resullts_limma.RData"))

LP_markers <- rownames(diff_LP)[which(diff_LP$logFC>3)] 
H_markers <- rownames(diff_H)[which(diff_H$logFC>3)] 
M_markers <- rownames(diff_M)[which(diff_M$logFC>3)] 
LP_markers <- LP_markers[which(LP_markers%in%rownames(sample))] # 133 genes
H_markers <- H_markers[which(H_markers%in%rownames(sample))] # 71 genes
M_markers <- M_markers[which(M_markers%in%rownames(sample))] # 300 genes

# TF identified by differential expression snRNAseq

TF_H <- geco.load("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/7.Gene_regulatory_networks/Tumor_cells/Merge/Analyses/TF_differential_expression/TF_H_snRNAseq_X_bulk.RData")
TF_LP <- geco.load("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/7.Gene_regulatory_networks/Tumor_cells/Merge/Analyses/TF_differential_expression/TF_LP_snRNAseq_X_bulk.RData")


##############################################
##### 1. FIND SAMPLE-SPECIFIC COMPONENTS #####
##############################################

if(identical(rownames(pca_coord), colnames(sample))){pc_vs_samples <- factoall(data.frame(cells=rownames(pca_coord), orig.ident=sample$orig.ident, PC=NA))}

pdf(file.path(resdir, "Sample_specific_components.pdf"))
for(pc in colnames(pca_coord)){
  pc_vs_samples$PC <- pca_coord[,pc]
  test <- kruskal.test(pc_vs_samples$PC, pc_vs_samples$orig.ident)
  
  print(ggplot(pc_vs_samples, aes(x=orig.ident, y=PC))+
    geom_violin()+
    ggtitle(paste0(pc, ", pval=", test$p.value)))
}
dev.off()


#### A REFAIRE EN PLUSIEURS ETAPES POUR COMPARER 
#sample_specif_PC <- c(3,4,5,7,8,9) # don't take into account PC1 because it is the mesenchymal component (n=1)

#nbcomp <- 1:15
#nbcomp <- setdiff(nbcomp, sample_specif_PC)

# run new UMAP without sample specific PCA components

#sample <- RunUMAP(sample, dims=nbcomp, verbose=F, umap.method = "uwot")
#DimPlot(sample, group.by = "orig.ident")

########################################################################
##### 2. IDENTIFY ADDITIONAL INTERESTING COMPONENTS TO PC1 AND PC2 #####
########################################################################

pdf(file.path(resdir, "PC1_PC2_additional_PC.pdf"))
for(pc in colnames(pca_coord)[3:ncol(pca_coord)]){
  print(FeaturePlot(sample, reduction="pca", dims=c(1,2), features=pc))
}
dev.off()


pdf(file.path(resdir, "PC_on_UMAP.pdf"))
for(pc in colnames(pca_coord)){
  print(FeaturePlot(sample, reduction="umap", features=pc))
}
dev.off()


######################################################################################
##### 3. CHECK IF H/LP/M CAN RECAPITULATE THE PCA OR IF THERE ARE OTHER SUBTYPES #####
######################################################################################

##### 1. H/LP/M scores #####
#---------------------------

# compute H/LP/M scores
DefaultAssay(sample) <- "RNA" # better to do means on RNA assay 
sample$mean_H <- apply(sample@assays$RNA@data[H_markers,], 2, mean); sample$mean_H <- geco.changeRange(sample$mean_H, newmin=0, newmax=1)
sample$mean_LP <- apply(sample@assays$RNA@data[LP_markers,], 2, mean); sample$mean_LP <- geco.changeRange(sample$mean_LP, newmin=0, newmax=1)
sample$mean_M <- apply(sample@assays$RNA@data[M_markers,], 2, mean); sample$mean_M <- geco.changeRange(sample$mean_M, newmin=0, newmax=1)

# compare H/LP/M scores to see if they identify all cells or if there are potential new subtypes

H_LP_M_scores <- factoall(data.frame(orig.ident=sample$orig.ident, subtype=sample$H_LP_M_clean_groups_1, cells=colnames(sample), H_score=sample$mean_H, LP_score=sample$mean_LP, M_score=sample$mean_M))

pdf(file.path(resdir, "H_vs_LP_vs_M_scores.pdf"))
ggplot(H_LP_M_scores, aes(x=H_score, y=LP_score, color=orig.ident))+
  geom_point(size=0.8)

ggplot(H_LP_M_scores, aes(x=H_score, y=LP_score, color=M_score))+
  geom_point(size=0.8)

# remove mesenchymal sample to see
ggplot(H_LP_M_scores[which(H_LP_M_scores$orig.ident!="CHC3610T"),], aes(x=H_score, y=LP_score, color=orig.ident))+
  geom_point(size=0.8)

ggplot(H_LP_M_scores[which(H_LP_M_scores$orig.ident!="CHC3610T"),], aes(x=H_score, y=LP_score, color=subtype))+
  geom_point(size=0.8)
dev.off()


##### 2. Intermediate PC: are some PCs linked to mixed H+LP subtype? #####
#-------------------------------------------------------------------------

if(identical(rownames(pca_coord), colnames(sample))){pc_vs_subtype <- factoall(data.frame(cells=rownames(pca_coord), subtype=sample$H_LP_M_clean_groups_1, PC=NA))}

pdf(file.path(resdir, "Link_PC_H_LP_M_subtypes.pdf"))
for(pc in colnames(pca_coord)){
  pc_vs_subtype$PC <- pca_coord[,pc]
  test <- kruskal.test(pc_vs_subtype$PC, pc_vs_subtype$subtype)
  
  print(ggplot(pc_vs_subtype, aes(x=subtype, y=PC))+
          geom_violin()+
          ggtitle(paste0(pc, ", pval=", test$p.value)))
}
dev.off()


# test with scores

H_LP_M_scores$sum_scores <- H_LP_M_scores$H_score + H_LP_M_scores$LP_score + H_LP_M_scores$M_score
H_LP_M_scores_save <- H_LP_M_scores
H_LP_M_scores$H_score <- H_LP_M_scores$H_score / H_LP_M_scores$sum_scores
H_LP_M_scores$LP_score <- H_LP_M_scores$LP_score / H_LP_M_scores$sum_scores
H_LP_M_scores$M_score <- H_LP_M_scores$M_score / H_LP_M_scores$sum_scores

pdf(file.path(resdir, "H_vs_LP_vs_M_scores_total_is_1.pdf"))
ggplot(H_LP_M_scores, aes(x=H_score, y=LP_score, color=M_score))+
  geom_point(size=0.8)

ggplot(H_LP_M_scores, aes(x=H_score, y=LP_score, color=orig.ident))+
  geom_point(size=0.8)

ggplot(H_LP_M_scores, aes(x=H_score, y=LP_score, color=subtype))+
  geom_point(size=0.8)
dev.off()

# test with PC genes

PC1_genes_pos <- as.data.frame(pca_genes) %>% arrange(desc(PC_1)) %>% row.names %>% head(n=100)
PC2_genes_pos <- as.data.frame(pca_genes) %>% arrange(desc(PC_2)) %>% row.names %>% head(n=100)
PC2_genes_neg <- as.data.frame(pca_genes) %>% arrange(PC_2) %>% row.names %>% head(n=100)

DefaultAssay(sample) <- "RNA" # better to do means on RNA assay 
sample$mean_PC1_pos <- apply(sample@assays$RNA@data[PC1_genes_pos,], 2, mean)
sample$mean_PC2_pos <- apply(sample@assays$RNA@data[PC2_genes_pos,], 2, mean)
sample$mean_PC2_neg <- apply(sample@assays$RNA@data[PC2_genes_neg,], 2, mean)

PC1_PC2_scores <- factoall(data.frame(orig.ident=sample$orig.ident, subtype=sample$H_LP_M_clean_groups_1, cells=colnames(sample), PC1_pos_score=sample$mean_PC1_pos, PC2_pos_score=sample$mean_PC2_pos, PC2_neg_score=sample$mean_PC2_neg))

ggplot(PC1_PC2_scores, aes(x=PC2_pos_score, y=PC2_neg_score, color=PC1_pos_score))+
  geom_point(size=0.8)

ggplot(PC1_PC2_scores, aes(x=PC2_pos_score, y=PC2_neg_score, color=orig.ident))+
  geom_point(size=0.8)

ggplot(PC1_PC2_scores, aes(x=PC2_pos_score, y=PC2_neg_score, color=subtype))+
  geom_point(size=0.8)



##############################################
##### 4. GENE CORRELATION TO PC1 AND PC2 #####
##############################################

# Idea: instead of taking the genes contributing to PC1 and PC2 from PCA analysis we take the genes that display correlation to these components (less sample-specific?)

DefaultAssay(sample) <- "RNA" # take RNA assay for correlation 

# create correlation table between gene expression and PC contributions
cor_gene_PCA <- as.data.frame(matrix(NA, nrow=nrow(sample@assays$RNA@counts), ncol=5))
colnames(cor_gene_PCA) <- c("gene", "PC1_cor", "PC1_pval", "PC2_cor", "PC2_pval") # we only focus on PC2 because this is where sample-specific effects cause problem
cor_gene_PCA$gene <- rownames(sample@assays$RNA@counts)

# before doing PC2 we need to remove M tumor cells because many H TF will not be expressed but the correlation will be badly impacted since M tumor cells are in the middle of PC2

sample_noM <- subset(sample, subset=orig.ident != "CHC3610T")
# for PC1 keep all samples because it is H+LP vs M

for(g in cor_gene_PCA$gene){
  a <- cor.test(sample@assays$RNA@data[g,], pca_coord[colnames(sample),"PC_1"])
  cor_gene_PCA[which(cor_gene_PCA$gene==g),"PC1_cor"] <- a$estimate
  cor_gene_PCA[which(cor_gene_PCA$gene==g),"PC1_pval"] <- a$p.value
  b <- cor.test(sample_noM@assays$RNA@data[g,], pca_coord[colnames(sample_noM),"PC_2"])
  cor_gene_PCA[which(cor_gene_PCA$gene==g),"PC2_cor"] <- b$estimate
  cor_gene_PCA[which(cor_gene_PCA$gene==g),"PC2_pval"] <- b$p.value
} 

# adjust p values to have q values
cor_gene_PCA$PC1_qval <- p.adjust(cor_gene_PCA$PC1_pval, method="BH")
cor_gene_PCA$PC2_qval <- p.adjust(cor_gene_PCA$PC2_pval, method="BH")
cor_gene_PCA <- cor_gene_PCA[,c("gene", "PC1_cor", "PC1_pval", "PC1_qval", "PC2_cor", "PC2_pval", "PC2_qval")]

save(cor_gene_PCA, file=file.path(resdir, "gene_vs_PC1_PC2_correlation_pval.RData"))

# take 100 genes most correlated
PC2_genes_pos <- cor_gene_PCA %>% arrange(desc(PC2_cor)) %>% pull(gene) %>% head(n=100)
PC2_genes_neg <- cor_gene_PCA %>% arrange(PC2_cor) %>% pull(gene) %>% head(n=100)

DefaultAssay(sample) <- "RNA" # better to do means on RNA assay 
sample$mean_PC2_pos <- apply(sample@assays$RNA@data[PC2_genes_pos,], 2, mean)
sample$mean_PC2_neg <- apply(sample@assays$RNA@data[PC2_genes_neg,], 2, mean)

PC2_scores <- factoall(data.frame(orig.ident=sample$orig.ident, subtype=sample$H_LP_M_clean_groups_1, cells=colnames(sample), PC2_pos_score=sample$mean_PC2_pos, PC2_neg_score=sample$mean_PC2_neg))

pdf(file.path(resdir, "Correlated_genes_PC2_top_100.pdf"))
ggplot(PC2_scores, aes(x=PC2_pos_score, y=PC2_neg_score, color=orig.ident))+
  geom_point(size=0.8)

ggplot(PC2_scores, aes(x=PC2_pos_score, y=PC2_neg_score, color=subtype))+
  geom_point(size=0.8)
dev.off()

### intersection with TF H and LP identified by differential expression

a <- intersect(TF_H, PC2_genes_neg)
b <- intersect(TF_LP, PC2_genes_pos)

save(a, file=file.path(resdir, "TF_H_X_PC2_neg_correlated_genes.RData"))
save(b, file=file.path(resdir, "TF_LP_X_PC2_pos_correlated_genes.RData"))

