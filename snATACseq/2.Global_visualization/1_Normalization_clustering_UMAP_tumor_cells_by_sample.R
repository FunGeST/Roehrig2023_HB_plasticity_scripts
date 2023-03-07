
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

name_sample <- c("CHC3133T", "CHC3377T", "CHC3662T")[3]


###########################
##### 0. DATA LOADING #####
###########################

### ATAC data 
project_all <- readRDS(file.path(resdir, "Save-ArchR-Project.rds"))
metadata <- getCellColData(project_all) # extract metadata

# Subset project to tumor cells of the sample of interest
cells_to_keep <- rownames(metadata)[which(!is.na(metadata$RNA_subtypes) & metadata$Sample==name_sample)]
project_s <- project_all[cells_to_keep,]

### CNV clusters)
CNV_clusters <- geco.load(file.path(cnv_dir, name_sample, "Final_clusters/final_CNV_clusters.RData"))
CNV_clusters$cells <- paste0(gsub("_", "#", CNV_clusters$cells), "-1") # take ArchR cell nomenclature
CNV_clusters <- CNV_clusters[match(rownames(project_s), CNV_clusters$cells),] # restrict to cells in ArchR

if(name_sample=="CHC3133T"){
  CNV_clusters[which(CNV_clusters$final_CNV_cluster == "+1q, 2, 4, 5, 6, 8, 12, 14, 16, 17, 19, 20"), "final_CNV_cluster"] <- "C1_3133T"
  CNV_clusters[which(CNV_clusters$final_CNV_cluster == "+1q, 2, 3, 4, 5, 6, 8, 12, 14, 16, 17, 19, 20"), "final_CNV_cluster"] <- "C2_3133T"
}else if(name_sample=="CHC3377T"){
   CNV_clusters[which(CNV_clusters$final_CNV_cluster == "+2q, 6, 12, 17, 20, X"), "final_CNV_cluster"] <- "C1_3377T"
   CNV_clusters[which(CNV_clusters$final_CNV_cluster == "+2q, 6, 12, 17, 20, X, 8"), "final_CNV_cluster"] <- "C2_3377T"
   CNV_clusters[which(CNV_clusters$final_CNV_cluster == "+2q, 6, 12, 17, 20, X, 1q"), "final_CNV_cluster"] <- "C3_3377T"
}else if(name_sample=="CHC3662T"){
   CNV_clusters[which(CNV_clusters$final_CNV_cluster == "+2q, 17, 20"), "final_CNV_cluster"] <- "C1_3662T"
   CNV_clusters[which(CNV_clusters$final_CNV_cluster == "+2q, 17, 20, 1q"), "final_CNV_cluster"] <- "C2_3662T"
   CNV_clusters[which(CNV_clusters$final_CNV_cluster == "2q, 17, 1q, 2p, 8, 12"), "final_CNV_cluster"] <- "C3_3662T"
}


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
if(identical(UMAP_df$cells, CNV_clusters$cells)){UMAP_df$CNV_clusters <- CNV_clusters$final_CNV_cluster} # Add CNV cluster info
if(name_sample=="CHC3133T"){UMAP_df$CNV_clusters <- factor(UMAP_df$CNV_clusters, levels=c("C1_3133T", "C2_3133T"))} # order CNV clusters for plot
# Select colors
if(name_sample=="CHC3133T"){col_sample = c("#00B52B", "#97DE08")
}else if(name_sample=="CHC3377T"){col_sample = c("#27DAED", "#119324", "#236CAF")
}else if(name_sample=="CHC3662T"){col_sample = c("#F559E1", "#BE57F7", "#5B53FB")}

### Plot
ggplot(UMAP_df, aes(x=UMAP_1, y=UMAP_2, col=CNV_clusters))+
  geom_point(size=1)+
  theme_classic()+
  scale_color_manual(values=col_sample)+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")
ggsave(file.path(resdir, paste0("UMAP_tumor_cells_clean/UMAP_", name_sample, "_CNV_clusters.png")), width=4, height=4)
ggplot(UMAP_df, aes(x=UMAP_1, y=UMAP_2, col=CNV_clusters))+
  geom_point(size=1)+
  theme_classic()+
  scale_color_manual(values=col_sample)+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
ggsave(file.path(resdir, paste0("UMAP_tumor_cells_clean/UMAP_", name_sample, "_CNV_clusters_with_legend.png")), width=4, height=4)

save(UMAP_df, file=file.path(resdir, paste0("UMAP_tumor_cells_clean/UMAP_", name_sample, "_CNV_clusters.RData")))

