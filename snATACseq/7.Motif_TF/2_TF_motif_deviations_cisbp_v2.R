
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
resdir = file.path("/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Motif_TF")
cisbp_dir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/0.Public_annotation/Cisbp_database/" # folder containing cisbp v2.00 results from MEME software


###########################
##### 0. DATA LOADING #####
###########################

project_all <- readRDS(file.path(archrdir, "Save-ArchR-Project.rds"))

### Add imputation weights
project_all <- addImputeWeights(project_all)

# Load cisbp v2.00 motifs
PWM_cisbp_v2 <- geco.load(file.path(cisbp_dir, "PWM_list_cisbp_v2.00.RData")) 


######################################
##### 1. SELECT BACKGROUND PEAKS #####
######################################

### Check that motif are present in the project
project_all <- addMotifAnnotations(ArchRProj = project_all, motifPWMs=PWM_cisbp_v2, name = "Motif_cisbp_v2.00") # takes ~5 min

# the function samples peaks based on similarity in GC-content and number of fragments across all samples using the Mahalanobis distance.
project_all <- addBgdPeaks(project_all, method = "ArchR", force=T) # we need to use method = ArchR and not method = chromVAR (default) because we work with non-fixed width peaks


#######################################################################
##### 2. COMPUTE PER CELL DEVIATIONS ACROSS ALL MOTIF ANNOTATIONS #####
#######################################################################

# create for each arrow file a deviation matrix called MotifMatrix
# takes some time: ~50 min

project_all <- addDeviationsMatrix(
  ArchRProj = project_all, 
  peakAnnotation = "Motif_cisbp_v2.00",
  force = TRUE
)

MotifMatrix <- getMatrixFromProject(
  ArchRProj = project_all,
  useMatrix = "Motif_cisbp_v2.00Matrix")
save(MotifMatrix, file=file.path(resdir, "MotifMatrix_cisbp_v2.00.RData"))
#From the above snapshot of the DataFrame, you can see that the seqnames of the MotifMatrix are not chromosomes. Typically, in matrices like the TileMatrix, PeakMatrix, and GeneScoreMatrix, we store the chromosome information in seqnames. The MotifMatrix does not have any corresponding position information but, instead, stores both the "devations" and "z-scores" from chromVAR into the same matrix using two different seqnames - deviations and z. This becomes important later on if you try to use the MotifMatrix (which is of class Sparse.Assays.Matrix) in functions such as getMarkerFeatures(). In these types of operations, ArchR will expect you to subset MotifMatrix to one of the two seqnames (i.e. select either z-scores or deviations to perform calculations).




