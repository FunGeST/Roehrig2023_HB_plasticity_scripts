
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ArchR)
library(pheatmap)

addArchRGenome("hg38") # load the genome for each R session
# addArchRThreads(threads = 16) # takes automatically half of the threads available

geco.load <- function (filename) # File to be loaded
{
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error - function geco.load : file ", filename, 
            " not found"))
  NULL
}


######################
##### PARAMETERS #####
######################

archrdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/"
resdir = file.path("/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/ATAC_tracks_per_subtype")

markers <- c("SLCO1B1", "APOB", "FRAS1", "FREM1", "PROM1", "TWIST2", "COL5A1")

###########################
##### 0. DATA LOADING #####
###########################

project_all <- readRDS(file.path(archrdir, "Save-ArchR-Project.rds"))
metadata <- getCellColData(project_all) # extract metadata


#####################################################
##### 1. ATAC TRACKS OF KNOWN GENES PER SUBTYPE #####
#####################################################

p <- plotBrowserTrack(
  ArchRProj = project_all, 
  groupBy = "RNA_subtypes", 
  useGroups = c("H", "H+LP", "LP", "M"),
  geneSymbol = markers, 
  upstream = 25000,
  downstream = 25000,
  tileSize=250,
  facetbaseSize=5,
  normMethod = "ReadsInTSS",
  log2Norm=T, 
  baseSize=10,
  pal=c("#E8AAD6", "#BC53CD", "#510787", "#6F3322")) ### /!\ Check that the coordinates given by ArchR are consistent with the hg38 coordinates of the genes found online (genecards, ensembl...)

plotPDF(plotList = p, 
        name = "ATAC_track_marker_genes.pdf", 
        ArchRProj = project_all, 
        addDOC = FALSE, width = 5, height = 5) # plotPDF generates the plots in the output directory of project_sample 







