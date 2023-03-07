
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ArchR)
library(pheatmap)

addArchRGenome("hg38") # load the genome for each R session
# addArchRThreads(threads = 16) # takes automatically half of the threads available


######################
##### PARAMETERS #####
######################

resdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged"

names_samples <- c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")


###########################
##### 0. DATA LOADING #####
###########################

project_all <- readRDS(file.path(resdir, "Save-ArchR-Project.rds"))


############################################
##### 1. DIMENSIONALITY REDUCTION: LSI #####
############################################

project_all <- addIterativeLSI(
  ArchRProj = project_all,
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

project_all <- addHarmony(
  ArchRProj = project_all,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force=T
)


#####################################
##### 3. CLUSTER IDENTIFICATION #####
#####################################

project_all <- addClusters(
  input = project_all,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force=T
)



#################################
##### 4. UMAP VISUALIZATION #####
#################################

project_all <- addUMAP(
  ArchRProj = project_all, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force=T
)

plotEmbedding(ArchRProj = project_all, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

# there can be 2 warning messages regarding invalid uid --> seemingly not problematic

UMAP_df <- getEmbedding(ArchRProj = project_all, embedding = "UMAP", returnDF = TRUE)
save(UMAP_df, file=file.path(resdir, "UMAP_coordinates.RData"))

pdf(file.path(resdir, "ArchR_results.pdf"))
plotEmbedding(ArchRProj = project_all, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
plotEmbedding(ArchRProj = project_all, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotEmbedding(ArchRProj = project_all, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAP")
dev.off()

###########################
##### 5. SAVE PROJECT #####
###########################

metadata <- getCellColData(project_all)
save(metadata, file=file.path(resdir, "metadata_all.RData"))

saveArchRProject(ArchRProj = project_all, outputDirectory = resdir, load = TRUE)



