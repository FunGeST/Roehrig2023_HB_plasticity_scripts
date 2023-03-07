
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ArchR)
library(Seurat)
library(ggplot2)
set.seed(1234)

geco.load <- function (filename) # File to be loaded
{
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error - function geco.load : file ", filename, 
            " not found"))
  NULL
}


addArchRGenome("hg38") # load the genome for each R session
# addArchRThreads(threads = 16) # takes automatically half of the threads available

######################
##### PARAMETERS #####
######################

datadir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/"
rna_dir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/2.Basic_analyses/SCTransform/No_regression/Merged/All_samples/"

resdir = file.path(datadir, "Peaks_to_genes"); if(!dir.exists(resdir)){dir.create(resdir)}

# parameters to choose
keepChr = paste0("chr", 1:22) # keep same chromosomes as ArchR, remove chr MT and Y but also X like Signac


###########################
##### 0. DATA LOADING #####
###########################

# Load ArchR project
project_all <- readRDS(file.path(datadir, "Save-ArchR-Project.rds"))

# Load RNA Seurat object
seurat_sample_all <- geco.load(file.path(rna_dir, "Seurat_object_analyses.RData")) # here we load the merge object with ALL cells (not only tumor)
seurat_sample_all <- RenameCells(seurat_sample_all, new.names = paste0(gsub("_", "#", colnames(seurat_sample_all)), "-1")) # replace by ArchR nomenclature

# Restrict ArchR project to common ATAC/RNA cells
project_RNA_ATAC <- project_all[which(project_all$cellNames %in% colnames(seurat_sample_all)),] # 15832 common cells = 90% of ArchR cells and 75% of Seurat cells

# Load ArchR gene annotation
gene_annot_archr <- getGeneAnnotation(project_all)$genes #extract from Archr gene annotation file, 25017 genes
gene_annot_archr <- gene_annot_archr[which(!is.na(gene_annot_archr$symbol))] # remove NA genes, 24961 genes 

# Load gene annotation
#gene_annot <- read.table("/mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/filtered_feature_bc_matrix/features.tsv.gz", sep="\t") 
#ensg_to_filter <- c("ENSG00000285053", "ENSG00000261186", "ENSG00000234323", "ENSG00000261480", "ENSG00000286070", "ENSG00000269226", "ENSG00000271858") # keep the ENSG id that are given by HGNC and remove the other ENSG for the same gene 
#gene_annot <- gene_annot[which(!gene_annot$V1 %in%ensg_to_filter),]
#genes_to_keep <- gene_annot[which(gene_annot$V3 == "Gene Expression" & gene_annot$V4 %in% keepChr),"V2"]


############################################
##### 1. ADD RNA DATA TO ARCHR PROJECT #####
############################################

# According to doc https://github.com/GreenleafLab/ArchR/issues/554 with multiome we don't need to run any gene integration function 

##### 1. Create SummarizedExperiment from Seurat object #####
#------------------------------------------------------------

### create Granges with genomic coordinates
#gene_coord <- data.frame(chr=NA, start=NA , end=NA , gene=intersect(genes_to_keep, rownames(seurat_sample_all)))
#gene_coord$chr <- gene_annot[match(gene_coord$gene, gene_annot$V2), "V4"]
#gene_coord$start <- gene_annot[match(gene_coord$gene, gene_annot$V2), "V5"]
#gene_coord$end <- gene_annot[match(gene_coord$gene, gene_annot$V2), "V6"]
#gr_gene_coord <- makeGRangesFromDataFrame(gene_coord, keep.extra.columns=TRUE, ignore.strand=TRUE,
#                                          seqnames.field="chr", start.field="start", end.field="end",
#                                          starts.in.df.are.0based=FALSE)

# keep only genes for which we have ArchR genomic annotation
seurat_sample_all_raw <- GetAssayData(seurat_sample_all[intersect(gene_annot_archr$symbol, rownames(seurat_sample_all)),], assay="RNA", slot="counts") # we first need to keep only raw counts (asked by ArchR https://www.archrproject.com/reference/addGeneExpressionMatrix.html)
if(identical(as.character(gene_annot_archr[which(gene_annot_archr$symbol%in%rownames(seurat_sample_all))]$symbol), intersect(gene_annot_archr$symbol, rownames(seurat_sample_all)))){
  seurat_sample_all_SE <- SummarizedExperiment::SummarizedExperiment(seurat_sample_all_raw, rowRanges = gene_annot_archr[which(gene_annot_archr$symbol %in%  rownames(seurat_sample_all))]) # convert Seurat object to single cell experiment (required by addGeneExpressionMatrix)
} # 19615 genes 



##### 2. Add Gene expression matrix in ArchR project #####
#---------------------------------------------------------

project_RNA_ATAC <- addGeneExpressionMatrix(input = project_RNA_ATAC, seRNA = seurat_sample_all_SE, force = TRUE)

exp_mat_in_ArchR <- getMatrixFromProject(ArchRProj = project_RNA_ATAC, useMatrix = "GeneExpressionMatrix") #19594 genes

# access to the gene set
geneSet <- getFeatures(project_RNA_ATAC, useMatrix="GeneExpressionMatrix") # 19594 genes
rownames(exp_mat_in_ArchR) <- geneSet

save(exp_mat_in_ArchR, file=file.path(datadir, "GeneExpressionMatrix_common_cells_RNA_ATAC.RData"))

#test on one gene
plotEmbedding(
  ArchRProj = project_RNA_ATAC,
  embedding = "UMAP",
  colorBy = "GeneExpressionMatrix",
  name = "FRAS1")


####################################
##### 2. INFER PEAK/GENE LINKS #####
####################################

##### 1. Peak2Gene function #####
#--------------------------------

project_RNA_ATAC <- addPeak2GeneLinks(ArchRProj = project_RNA_ATAC, useMatrix="GeneExpressionMatrix", reducedDims = "IterativeLSI")


##### 2. Extract peak/gene links #####
#-------------------------------------

peaks2gene <- getPeak2GeneLinks(ArchRProj = project_RNA_ATAC, corCutOff = -10, FDRCutOff = 10, varCutOffATAC = 0, varCutOffRNA = 0, resolution = 1, returnLoops = FALSE) #823378 rows, 86138 unique peaks, 14045 unique genes
# try to take all information without any thresholds (default: corCutOff=0.45, FDRCutOff=0.0001, varCutOffATAC=0.25, varCutOffRNA=0.25)
# resolution = 1, this links the center of the peak to the single-base TSS of the gene (onlyw orks with returnLoops =T)
# Decreasing the resolution even further also decreases the total number of peak-to-gene links identified.

peaks2gene_gr <- getPeak2GeneLinks(ArchRProj = project_RNA_ATAC, corCutOff = -10, FDRCutOff = 10, varCutOffATAC = 0, varCutOffRNA = 0, resolution = 1, returnLoops = TRUE) #823378 rows, 86138 unique peaks, 14045 unique genes

### Annotate peaks and genes
peaks2gene_df <- as.data.frame(peaks2gene)
peaks2gene_df$peak <- paste0(seqnames(metadata(peaks2gene)[[1]][peaks2gene$idxATAC]), ":", start(metadata(peaks2gene)[[1]][peaks2gene$idxATAC]), "-", end(metadata(peaks2gene)[[1]][peaks2gene$idxATAC]))
peaks2gene_df$gene <- geneSet[peaks2gene$idxRNA]


################
##### SAVE #####
################

save(peaks2gene, file=file.path(resdir, "peaks2genes_all_correlations.RData"))
save(peaks2gene_df, file=file.path(resdir, "peaks2genes_all_correlations_df.RData"))
save(peaks2gene_gr, file=file.path(resdir, "peaks2genes_all_correlations_gr.RData"))


#saveArchRProject(ArchRProj = project_RNA_ATAC, outputDirectory = resdir, load = TRUE)




