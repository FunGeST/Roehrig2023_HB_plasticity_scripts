
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



######################
##### PARAMETERS #####
######################

datadir = paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/")

resdir = "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/10.Cancer_stem_cells/"; if(!dir.exists(resdir)){dir.create(resdir)}
csc_markers_liver_main <- c("PROM1", "THY1", "CD44", "EPCAM", "CD24") # liver CSC markers that are the most recurrent


###########################
##### 0. DATA LOADING #####
###########################

### Load seurat object for merged tumor samples (tumor cells only)
sample_seurat <- geco.load(file.path(datadir, "Seurat_object_analyses.RData"))
DefaultAssay(sample_seurat) <- "RNA" # always set RNA assay for differential analyses

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
  

########################################################
##### 1. DIFFERENTIAL EXPRESSION: LP PROM1+/PROM1- #####
########################################################

##### 1.  Define the groups to compare #####
#-------------------------------------------

csc_PROM1_plus <- names(potential_csc[,which(potential_csc["PROM1",]==1)]); length(csc_PROM1_plus) # 71 CSC for PROM1
LP_PROM1_minus <- setdiff(rownames(sample_seurat@meta.data[which(sample_seurat@meta.data$H_LP_M_clean_groups_1=="LP"),]), csc_PROM1_plus); length(LP_PROM1_minus) # 3083 LP cells which have <5 SCT-corrected UMI counts for PROM1

# assign to seurat object
sample_seurat$PROM1_csc <- NA
sample_seurat@meta.data[csc_PROM1_plus,"PROM1_csc"] <- "CSC PROM1+"
sample_seurat@meta.data[LP_PROM1_minus,"PROM1_csc"] <- "LP PROM1-"

# check subtype of those csc
table(sample_seurat@meta.data[intersect(csc_PROM1_plus, rownames(sample_seurat@meta.data)),"H_LP_M_clean_groups_1"])


##### 2. Differential expression with FindMarkers #####
#------------------------------------------------------

Idents(sample_seurat) <- "PROM1_csc" # set the column for differential analysis
diff_PROM1_csc <- FindMarkers(sample_seurat, ident.1 = "CSC PROM1+", ident.2 = "LP PROM1-")
diff_PROM1_csc <- diff_PROM1_csc %>% filter(p_val_adj < 0.05) # keep siginificant results: 99 genes
save(diff_PROM1_csc, file=file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/diff_exp_PROM1_csc.RData"))


##### 3. Display top markers on the UMAP/PCA
#-------------------------------------------

top_10_markers_pos <- diff_PROM1_csc %>% arrange(desc(avg_logFC)) %>% row.names() %>% head(10) 
top_10_markers_neg <- diff_PROM1_csc %>% arrange(avg_logFC) %>% row.names() %>% head(10) 

pdf(file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/diff_exp_PROM1_csc_top10_markers_pos.pdf"))
for(g in top_10_markers_pos){
  print(FeaturePlot(sample_seurat, features=g, reduction="umap") + scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20')))
  print(FeaturePlot(sample_seurat, features=g, reduction="pca") + scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20')))
}
dev.off()

pdf(file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/diff_exp_PROM1_csc_top10_markers_neg.pdf"))
for(g in top_10_markers_neg){
  print(FeaturePlot(sample_seurat, features=g, reduction="umap") + scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20')))
  print(FeaturePlot(sample_seurat, features=g, reduction="pca") + scale_color_gradientn( colours = c('#FBF6B8', '#F5D8AE', '#E88976', '#9B4577', '#1B0A20')))
}
dev.off()


##### 4. Pathway analysis #####
#------------------------------

# seprate with or without PROM1

### PROM1 + CSC
enrich.test_PROM1_plus <- geco.enrichmentTest(gene.sets=MSIG.ls, mylist= rownames(diff_PROM1_csc)[which(diff_PROM1_csc$avg_logFC>0)], possibleIds= rownames(sample_seurat@assays$RNA@counts))
save(enrich.test_PROM1_plus, file=file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/diff_exp_PROM1_csc_pos_enrichment_test.RData"))

### PROM1 - LP
enrich.test_LP_PROM1_minus <- geco.enrichmentTest(gene.sets=MSIG.ls, mylist= rownames(diff_PROM1_csc)[which(diff_PROM1_csc$avg_logFC<0)], possibleIds= rownames(sample_seurat@assays$RNA@counts))
save(enrich.test_LP_PROM1_minus, file=file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/diff_exp_PROM1_csc_neg_enrichment_test.RData"))


######################################################
##### 2. PROJECT SNRNASEQ CSC SIGNATURES IN BULK #####
######################################################

##### 1. Extract CSC signature #####
#-----------------------------------

# Signature = genes diff exp >0 in CSC for the marker

csc_signature <- diff_PROM1_csc %>% filter(avg_logFC > 0) %>% row.names() # 47 genes
setdiff(csc_signature, rownames(exp)) # check if all genes in the CSC signature are found in bulk
# LRATD1 is not found in bulk but another name is FAM84A and this symbol is found in the bulk exp table
csc_signature <- gsub("LRATD1", "FAM84A", csc_signature) # replace this gene by its other name
setdiff(csc_signature, rownames(exp)) # must be 0 now


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
save(bulk_mean_exp_csc_sig, file=file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/PROM1_csc_signature_bulk_mean_exp.RData"))

##### 3. Distribution of mean bulk expression #####
#--------------------------------------------------

# Identify the bulk mean exp for the multiome tumors to display on the density curve

pos_2959T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC2959T"),"mean_exp"]
pos_2960T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC2960T"),"mean_exp"]
pos_3133T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC3133T"),"mean_exp"]
pos_3377T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC3377T"),"mean_exp"]
pos_3610T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC3610T"),"mean_exp"]
pos_3662T <- bulk_mean_exp_csc_sig[which(bulk_mean_exp_csc_sig$sample=="CHC3662T"),"mean_exp"]


# Plot density curve, group by bulk cluster

pdf(file.path(resdir, "Potential_CSC_PROM1_THY1_EPCAM_CD24_CD44/Differential_expression/PROM1_csc_signature_bulk_mean_exp.pdf"))
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
dev.off()
