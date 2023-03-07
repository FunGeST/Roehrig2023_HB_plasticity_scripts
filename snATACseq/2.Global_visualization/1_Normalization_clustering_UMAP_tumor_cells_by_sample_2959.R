
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ArchR)
library(pheatmap)
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

resdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged"
cnv_dir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/"

name_sample <- c("CHC2959T", "CHC2960T")


###########################
##### 0. DATA LOADING #####
###########################

### ATAC data 
project_all <- readRDS(file.path(resdir, "Save-ArchR-Project.rds"))
metadata <- getCellColData(project_all) # extract metadata

# Subset project to tumor cells of the sample of interest
cells_to_keep <- rownames(metadata)[which(!is.na(metadata$RNA_subtypes) & metadata$Sample%in%name_sample)]
project_s <- project_all[cells_to_keep,]

### CNV clusters)
CNV_clusters <- geco.load(file.path(cnv_dir, "CHC2959T_CHC2960T/Final_clusters/clusters_euclidean_ward.D.RData"))
names(CNV_clusters) <- paste0(gsub("_", "#", names(CNV_clusters)), "-1") # take ArchR cell nomenclature
if(length(setdiff(rownames(project_s), names(CNV_clusters)))==0){CNV_clusters <- CNV_clusters[rownames(project_s)] }# restrict to cells in ArchR
CNV_clusters[which(CNV_clusters==1)] <- "C2_2959" # rename clsuters, beware here we changed and put 2960T as clone1 and 2959t as clone2
CNV_clusters[which(CNV_clusters==2)] <- "C1_2959"


############################################
##### 1. DIMENSIONALITY REDUCTION: LSI #####
############################################

project_s <- addIterativeLSI(
  ArchRProj = project_s,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = 0.2, 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  threads=1,
  force=T
)


#######################################
##### 2. HARMONY BATCH CORRECTION #####
#######################################

project_s <- addHarmony(
  ArchRProj = project_s,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force=T
)


#################################
##### 3. UMAP VISUALIZATION #####
#################################

project_s <- addUMAP(
  ArchRProj = project_s, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force=T
)
# there can be 2 warning messages regarding invalid uid --> seemingly not problematic

### Extract UMAP coordinates
UMAP_df <- getEmbedding(ArchRProj = project_s, embedding = "UMAP", returnDF = TRUE); colnames(UMAP_df) <- c("UMAP_1", "UMAP_2")
UMAP_df$cells <- rownames(UMAP_df)
UMAP_df$RNA_subtypes <- metadata[UMAP_df$cell, "RNA_subtypes"]
if(identical(UMAP_df$cells, names(CNV_clusters))){UMAP_df$CNV_clusters <- CNV_clusters} # Add CNV cluster info

### Plot
ggplot(UMAP_df, aes(x=UMAP_1, y=UMAP_2, col=CNV_clusters))+
  geom_point(size=0.4)+
  theme_classic()+
  scale_color_manual(values=c("#B79F00", "#F6493C"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")
ggsave(file.path(resdir, paste0("UMAP_tumor_cells_clean/UMAP_CHC2959T_CHC2960T_CNV_clusters.png")), width=4, height=4)
ggplot(UMAP_df, aes(x=UMAP_1, y=UMAP_2, col=CNV_clusters))+
  geom_point(size=0.4)+
  theme_classic()+
  scale_color_manual(values=c("#B79F00", "#F6493C"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
ggsave(file.path(resdir, paste0("UMAP_tumor_cells_clean/UMAP_CHC2959T_CHC2960T_CNV_clusters_with_legend.png")), width=4, height=4)

save(UMAP_df, file=file.path(resdir, paste0("UMAP_tumor_cells_clean/UMAP_CHC2959T_CHC2960T_CNV_clusters.RData")))

