
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

name_sample <- c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3377T", "CHC3610T", "CHC3662T")


###########################
##### 0. DATA LOADING #####
###########################

### ATAC data 
project_all <- readRDS(file.path(resdir, "Save-ArchR-Project.rds"))
metadata <- getCellColData(project_all) # extract metadata

# Subset project to tumor cells of the sample of interest
cells_to_keep <- rownames(metadata)[which(!is.na(metadata$RNA_subtypes) & metadata$Sample%in%name_sample)]
project_s <- project_all[cells_to_keep,]


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
UMAP_df$RNA_subtypes <- factor(UMAP_df$RNA_subtypes, levels=c("H", "H+LP", "LP", "M")) # order subtypes for plot
UMAP_df <- UMAP_df[sample(rownames(UMAP_df)),] # shuffle randomly the cell position in the table to have random overlap of the points in the grahp

### Plot
ggplot(UMAP_df, aes(x=UMAP_1, y=UMAP_2, col=RNA_subtypes))+
  geom_point(size=0.4)+
  theme_classic()+
  scale_color_manual(values=c("#E8AAD6", "#BC53CD", "#510787", "#6F3322"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")
ggsave(file.path(resdir, "UMAP_tumor_cells_clean/UMAP_all_tumors_RNA_subtypes.png"), width=4, height=4)
ggplot(UMAP_df, aes(x=UMAP_1, y=UMAP_2, col=RNA_subtypes))+
  geom_point(size=0.4)+
  theme_classic()+
  scale_color_manual(values=c("#E8AAD6", "#BC53CD", "#510787", "#6F3322"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
ggsave(file.path(resdir, "UMAP_tumor_cells_clean/UMAP_all_tumors_RNA_subtypes_with_legend.png"), width=4, height=4)

save(UMAP_df, file=file.path(resdir, "UMAP_tumor_cells_clean/UMAP_all_tumors_RNA_subtypes.RData"))

