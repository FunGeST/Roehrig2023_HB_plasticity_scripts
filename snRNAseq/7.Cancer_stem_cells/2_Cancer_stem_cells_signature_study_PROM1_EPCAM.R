
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

geco.changeRange <- function (v, newmin = 1, newmax = 10) 
{
  oldmin <- min(v, na.rm = TRUE)
  oldmax <- max(v, na.rm = TRUE)
  newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin))
}


######################
##### PARAMETERS #####
######################

datadir = paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/")

resdir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/10.Cancer_stem_cells/"; if(!dir.exists(resdir)){dir.create(resdir)}
csc_markers_liver_main <- c("PROM1", "THY1", "CD44", "EPCAM", "CD24") # liver CSC markers that are the most recurrent

name_samples <- c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3377T", "CHC3610T", "CHC3662T")


###########################
##### 0. DATA LOADING #####
###########################

### Load seurat object for merged tumor samples (tumor cells only)
sample_seurat <- geco.load(file.path(datadir, "Seurat_object_analyses.RData"))
DefaultAssay(sample_seurat) <- "RNA" # always set RNA assay for differential analyses

### Define cells to remove (cell cycle, MT% clusters)
cells_to_keep <- c()
for(n in name_samples){
  sample <- geco.load(file.path(datadir, n, "Seurat_object_analyses.RData")) # load object for the sample
  clusters_to_remove <- unique(grep(c("MT%|prolif"), sample$H_LP_M_seurat_clusters, value=T)) # remove cycle and MT% clusters
  cells_to_keep <- c(cells_to_keep, rownames(sample@meta.data[which(!sample@meta.data$H_LP_M_seurat_clusters%in%clusters_to_remove),]))
}
sample_seurat$cell_ID <- rownames(sample_seurat@meta.data)
sample_seurat <- subset(sample_seurat, subset=cell_ID %in% cells_to_keep) # remove the cells from the cell cycle and MT clusters (defined on the individual bases) in the merged object = 13255 cells
sample_seurat$cell_ID <- NULL

### Load potential cancer stem cells
potential_csc <- geco.load(file.path(resdir, "Liver_main_CSC_markers_potential_csc_min_5_UMI_SCT.RData")) # 427 potential csc, 5 main csc markers, 1 if the cell is csc for this marker (>=5 SCT-corrected UMI counts)

### Load file for enrichement test
MSIG.ls <- geco.load("F:/Amelie-Datas/Data-hepatopartage/MSIG_v6.1.RData")

### Load bulk files

load("F:/Amelie-Datas/Data-hepatopartage/Expression_matrix_new_2020/exp.Rdata")
load("F:/Amelie-Datas/Data-hepatopartage/Expression_matrix_new_2020/RNA_order.Rdata")
# identify clusters
NT_samples <- RNA_order[1:4]
H_samples <- RNA_order[5:48]
LP_samples <- RNA_order[49:92]
M_samples <- RNA_order[93:104]
  

#########################################################################
##### 1. DIFFERENTIAL EXPRESSION: LP PROM1+ or EPCAM+/PROM1- EPCAM- #####
#########################################################################

##### 1.  Define the groups to compare #####
#-------------------------------------------

csc_PROM1_EPCAM_plus <- names(potential_csc[,which(potential_csc["PROM1",]==1 | potential_csc["EPCAM",]==1)]); length(csc_PROM1_EPCAM_plus) # 113 CSC for PROM1 and/or EPCAM
# WARNING we restrict to only LP cells (there are some H+LP cells in the PROM1/EPCAM CSC)
csc_PROM1_EPCAM_plus <- intersect(rownames(sample_seurat@meta.data[which(sample_seurat@meta.data$H_LP_M_clean_groups_1=="LP"),]), csc_PROM1_EPCAM_plus); length(csc_PROM1_EPCAM_plus) # 96 CSC for PROM1 and/or EPCAM that are also LP
LP_PROM1_EPCAM_minus <- setdiff(rownames(sample_seurat@meta.data[which(sample_seurat@meta.data$H_LP_M_clean_groups_1=="LP"),]), csc_PROM1_EPCAM_plus); length(LP_PROM1_EPCAM_minus) # 2341 LP cells which have <5 SCT-corrected UMI counts for PROM1 and/or EPCAM

# assign to seurat object
sample_seurat$PROM1_EPCAM_csc <- NA
sample_seurat@meta.data[csc_PROM1_EPCAM_plus,"PROM1_EPCAM_csc"] <- "CSC PROM1+ and/or EPCAM+"
sample_seurat@meta.data[LP_PROM1_EPCAM_minus,"PROM1_EPCAM_csc"] <- "LP PROM1- EPCAM-"

# check subtype of those csc
table(sample_seurat@meta.data[intersect(csc_PROM1_EPCAM_plus, rownames(sample_seurat@meta.data)),"H_LP_M_clean_groups_1"]) # 96 LP cells


##### 2. Differential expression with FindMarkers #####
#------------------------------------------------------

Idents(sample_seurat) <- "PROM1_EPCAM_csc" # set the column for differential analysis
diff_PROM1_EPCAM_csc <- FindMarkers(sample_seurat, ident.1 = "CSC PROM1+ and/or EPCAM+", ident.2 = "LP PROM1- EPCAM-")
diff_PROM1_EPCAM_csc <- diff_PROM1_EPCAM_csc %>% filter(p_val_adj < 0.05) # keep siginificant results: 165 genes
save(diff_PROM1_EPCAM_csc, file=file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/diff_exp_PROM1_EPCAM_csc.RData"))


##### 3. Display top markers on the UMAP/PCA
#-------------------------------------------

top_10_markers_pos <- diff_PROM1_EPCAM_csc %>% arrange(desc(avg_logFC)) %>% row.names() %>% head(10) 
top_10_markers_neg <- diff_PROM1_EPCAM_csc %>% arrange(avg_logFC) %>% row.names() %>% head(10) 

pdf(file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/diff_exp_PROM1_EPCAM_csc_top10_markers_pos.pdf"))
for(g in top_10_markers_pos){
  print(FeaturePlot(sample_seurat, features=g, reduction="umap") + scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20')))
  print(FeaturePlot(sample_seurat, features=g, reduction="pca") + scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20')))
}
dev.off()

pdf(file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/diff_exp_PROM1_EPCAM_csc_top10_markers_neg.pdf"))
for(g in top_10_markers_neg){
  print(FeaturePlot(sample_seurat, features=g, reduction="umap") + scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20')))
  print(FeaturePlot(sample_seurat, features=g, reduction="pca") + scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20')))
}
dev.off()


##### 4. Pathway analysis #####
#------------------------------

### PROM1 + and/or EPCAM + CSC
enrich.test_PROM1_EPCAM_plus <- geco.enrichmentTest(gene.sets=MSIG.ls, mylist= rownames(diff_PROM1_EPCAM_csc)[which(diff_PROM1_EPCAM_csc$avg_logFC>0)], possibleIds= rownames(sample_seurat@assays$RNA@counts))
save(enrich.test_PROM1_EPCAM_plus, file=file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/diff_exp_PROM1_EPCAM_csc_pos_enrichment_test.RData"))

### PROM1 - EPCAM- LP
enrich.test_LP_PROM1_EPCAM_minus <- geco.enrichmentTest(gene.sets=MSIG.ls, mylist= rownames(diff_PROM1_EPCAM_csc)[which(diff_PROM1_EPCAM_csc$avg_logFC<0)], possibleIds= rownames(sample_seurat@assays$RNA@counts))
save(enrich.test_LP_PROM1_EPCAM_minus, file=file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/diff_exp_PROM1_EPCAM_csc_neg_enrichment_test.RData"))


############################################################
##### 2. IDENTIFY LP WITH/WITHOUT CSC: PROM1 AND EPCAM #####
############################################################

##### 1. Mean bulk expression per sample #####
#---------------------------------------------

bulk_mean_exp_csc_PROM1_EPCAM <- factoall(data.frame(sample=RNA_order, cluster=NA, mean_exp=NA)) # create table for bulk samples
# add cluster info
bulk_mean_exp_csc_PROM1_EPCAM$cluster <- ifelse(bulk_mean_exp_csc_PROM1_EPCAM$sample %in% NT_samples, "NT", 
                                        ifelse(bulk_mean_exp_csc_PROM1_EPCAM$sample %in% H_samples, "H", 
                                               ifelse(bulk_mean_exp_csc_PROM1_EPCAM$sample %in% LP_samples, "LP",
                                                      ifelse(bulk_mean_exp_csc_PROM1_EPCAM$sample %in% M_samples, "M", ""))))
# compute mean expression on CSC gene signature in bulk
bulk_mean_exp_csc_PROM1_EPCAM$mean_exp <- apply(exp[c("PROM1", "EPCAM"),RNA_order],2,mean) 
save(bulk_mean_exp_csc_PROM1_EPCAM, file=file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/PROM1_EPCAM_bulk_mean_exp.RData"))


##### 2. Distribution of mean bulk expression #####
#--------------------------------------------------

### Identify the bulk mean exp for the multiome tumors to display on the density curve

pos_2959T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC2959T"),"mean_exp"]
pos_2960T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC2960T"),"mean_exp"]
pos_3133T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC3133T"),"mean_exp"]
pos_3377T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC3377T"),"mean_exp"]
pos_3610T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC3610T"),"mean_exp"]
pos_3662T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC3662T"),"mean_exp"]

### Plot

pdf(file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/PROM1_EPCAM_bulk_mean_exp.pdf"))

# Density curves
ggplot(bulk_mean_exp_csc_PROM1_EPCAM, aes(x=mean_exp, color=cluster)) + 
  geom_density(size=1.5)+
  scale_color_manual(values=c("#E8AAD6", "#510787", "#6F3323", "grey")) +
  theme_classic()+
  geom_vline(xintercept = 9, linetype="dashed", col="red") # select a threshold to separate LP with/without csc
ggplot(bulk_mean_exp_csc_PROM1_EPCAM, aes(x=mean_exp, color=cluster)) + 
  geom_density(size=1.5)+
  scale_color_manual(values=c("#E8AAD6", "#510787", "#6F3323", "grey")) +
  theme_classic()+
  geom_vline(xintercept = pos_2959T, linetype="dashed")+
  geom_text(aes(x=pos_2959T, y=0.05, label="CHC2959T"), angle=90, color="black")+
  geom_vline(xintercept = pos_2960T, linetype="dashed")+
  geom_text(aes(x=pos_2960T, y=0.05, label="CHC2960T"), angle=90, color="black")+
  geom_vline(xintercept = pos_3133T, linetype="dashed")+
  geom_text(aes(x=pos_3133T, y=0.05, label="CHC3133T"), angle=90, color="black")+
  geom_vline(xintercept = pos_3377T, linetype="dashed")+
  geom_text(aes(x=pos_3377T, y=0.05, label="CHC3377T"), angle=90, color="black")+
  geom_vline(xintercept = pos_3610T, linetype="dashed")+
  geom_text(aes(x=pos_3610T, y=0.05, label="CHC3610T"), angle=90, color="black")+
  geom_vline(xintercept = pos_3662T, linetype="dashed")+
  geom_text(aes(x=pos_3662T, y=0.05, label="CHC3662T"), angle=90, color="black")

# Stripchart
ggplot(bulk_mean_exp_csc_PROM1_EPCAM, aes(x=cluster, y=mean_exp, color=cluster)) + 
  geom_boxplot()+
  geom_jitter(size=3, position=position_jitter(0.2))+
  scale_color_manual(values=c("#E8AAD6", "#510787", "#6F3323", "grey")) +
  theme_classic()
ggplot(bulk_mean_exp_csc_PROM1_EPCAM, aes(x=cluster, y=mean_exp, color=cluster)) + 
  geom_violin(trim=F)+
  geom_jitter(size=3, position=position_jitter(0.2))+
  scale_color_manual(values=c("#E8AAD6", "#510787", "#6F3323", "grey")) +
  theme_classic()+
  geom_hline(yintercept = 9, linetype="dashed", col="red")

dev.off()



#############################################################################
##### 2. IDENTIFY LP WITH/WITHOUT CSC: CSC SIGNATURE FROM PROM1+/EPCAM+ #####
#############################################################################

##### 1. Extract CSC signature #####
#-----------------------------------

# Signature = genes diff exp >0 in CSC for the marker

csc_signature <- diff_PROM1_EPCAM_csc %>% filter(avg_logFC > 0) %>% row.names(); length(csc_signature) # 103 genes
setdiff(csc_signature, rownames(exp)) # check if all genes in the CSC signature are found in bulk
csc_signature <- intersect(csc_signature, rownames(exp)); length(csc_signature) # 102 genes


##### 2. Mean bulk expression per sample #####
#---------------------------------------------

bulk_mean_exp_csc_sig <- factoall(data.frame(sample=RNA_order, cluster=NA, mean_exp=NA)) # create table for bulk samples
# add cluster info
bulk_mean_exp_csc_sig$cluster <- ifelse(bulk_mean_exp_csc_sig$sample %in% NT_samples, "NT", 
                                        ifelse(bulk_mean_exp_csc_sig$sample %in% H_samples, "H", 
                                               ifelse(bulk_mean_exp_csc_sig$sample %in% LP_samples, "LP",
                                                      ifelse(bulk_mean_exp_csc_sig$sample %in% M_samples, "M", ""))))
# compute mean expression on CSC gene signature in bulk
bulk_mean_exp_csc_sig$mean_exp <- apply(exp[csc_signature,RNA_order],2,mean) # all genes of csc_signature must be in the rownames of exp
save(bulk_mean_exp_csc_sig, file=file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/PROM1_EPCAM_csc_signature_bulk_mean_exp.RData"))


##### 3. Heatmap of the CSC signature genes in bulk and snRNAseq #####
#---------------------------------------------------------------------

### Set the expression of each gene between 0 and 1
exp_csc <- as.matrix(t(exp[csc_signature,RNA_order]))
for(j in 1:ncol(exp_csc)){
  exp_csc[,j] <- geco.changeRange( exp_csc[,j], newmin = 0, newmax = 1)
}

hmColors <- colorRampPalette(c("royalblue","white","indianred"))(256) # define the scale colors

pdf(file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/PROM1_EPCAM_csc_signature_heatmap_mean_exp.pdf"))
par(mar=c(5,5,2,2))
image(exp_csc[,ncol(exp_csc):1], xaxt= "n", yaxt= "n", col = hmColors)
axis( 2, at=seq(0,1,length.out=ncol(exp_csc[,ncol(exp_csc):1]) ), labels=colnames(exp_csc[,ncol(exp_csc):1]), las= 2, cex.axis=0.3)
axis( 1, at=seq(0,1,length.out=nrow(exp_csc) ), labels=rownames(exp_csc), las=2, cex.axis=0.3)

### Scale the expected genes in snRNAseq
sample_seurat <- ScaleData(sample_seurat, features=csc_signature)
DoHeatmap(sample_seurat, group.by="H_LP_M_clean_groups_1",  features=csc_signature)+ scale_fill_gradientn(colors = c("blue", "white", "red"))+ theme(text = element_text(size = 5))
DoHeatmap(sample_seurat, group.by="orig.ident", features=csc_signature)+ scale_fill_gradientn(colors = c("blue", "white", "red"))+ theme(text = element_text(size = 5))
dev.off()


##### 4. Distribution of mean bulk expression #####
#--------------------------------------------------

### Identify the bulk mean exp for the multiome tumors to display on the density curve

pos_2959T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC2959T"),"mean_exp"]
pos_2960T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC2960T"),"mean_exp"]
pos_3133T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC3133T"),"mean_exp"]
pos_3377T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC3377T"),"mean_exp"]
pos_3610T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC3610T"),"mean_exp"]
pos_3662T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC3662T"),"mean_exp"]

### Plot

pdf(file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/PROM1_EPCAM_csc_signature_bulk_mean_exp.pdf"))

# Density curves
ggplot(bulk_mean_exp_csc_sig, aes(x=mean_exp, color=cluster)) + 
  geom_density(size=1.5)+
  scale_color_manual(values=c("#E8AAD6", "#510787", "#6F3323", "grey")) +
  theme_classic()
ggplot(bulk_mean_exp_csc_sig, aes(x=mean_exp, color=cluster)) + 
  geom_density(size=1.5)+
  scale_color_manual(values=c("#E8AAD6", "#510787", "#6F3323", "grey")) +
  theme_classic()+
  geom_vline(xintercept = pos_2959T, linetype="dashed")+
  geom_text(aes(x=pos_2959T, y=0.05, label="CHC2959T"), angle=90, color="black")+
  geom_vline(xintercept = pos_2960T, linetype="dashed")+
  geom_text(aes(x=pos_2960T, y=0.05, label="CHC2960T"), angle=90, color="black")+
  geom_vline(xintercept = pos_3133T, linetype="dashed")+
  geom_text(aes(x=pos_3133T, y=0.05, label="CHC3133T"), angle=90, color="black")+
  geom_vline(xintercept = pos_3377T, linetype="dashed")+
  geom_text(aes(x=pos_3377T, y=0.05, label="CHC3377T"), angle=90, color="black")+
  geom_vline(xintercept = pos_3610T, linetype="dashed")+
  geom_text(aes(x=pos_3610T, y=0.05, label="CHC3610T"), angle=90, color="black")+
  geom_vline(xintercept = pos_3662T, linetype="dashed")+
  geom_text(aes(x=pos_3662T, y=0.05, label="CHC3662T"), angle=90, color="black")

# Stripchart
ggplot(bulk_mean_exp_csc_sig, aes(x=cluster, y=mean_exp, color=cluster)) + 
  geom_boxplot()+
  geom_jitter(size=3, position=position_jitter(0.2))+
  scale_color_manual(values=c("#E8AAD6", "#510787", "#6F3323", "grey")) +
  theme_classic()
ggplot(bulk_mean_exp_csc_sig, aes(x=cluster, y=mean_exp, color=cluster)) + 
  geom_violin(trim=F)+
  geom_jitter(size=3, position=position_jitter(0.2))+
  scale_color_manual(values=c("#E8AAD6", "#510787", "#6F3323", "grey")) +
  theme_classic()

dev.off()
