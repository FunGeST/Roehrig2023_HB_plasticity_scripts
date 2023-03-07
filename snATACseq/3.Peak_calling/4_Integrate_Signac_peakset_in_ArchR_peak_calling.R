
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ArchR)
library(Seurat)
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

datadir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged"
peak_dir = file.path(datadir, "/Signac/2.Peak_calling")


###########################
##### 0. DATA LOADING #####
###########################

# Load ArchR project
project_all <- readRDS(file.path(datadir, "Save-ArchR-Project.rds"))
atac_metadata <- getCellColData(project_all) # extract ArchR cell metadata

# Load Signac peak set
peak_set_signac <- geco.load(file.path(peak_dir, "combined_peaks_all_samples.RData"))


###################################################
##### 1. ADD SIGNAC PEAK SET TO ARCHR PROJECT #####
###################################################

addPeakSet(
  ArchRProj = project_all,
  peakSet = peak_set_signac,
  force = TRUE
) # WARNING: might generate a NULL peakset --> if so just run the function line per line, it works then

# histogram of peak width

pdf(file.path(datadir, "Peak_size_histogram.pdf"))
hist(width(peak_set_signac), xlab="peak size", breaks=200)
abline(v=median(width(peak_set_signac)), col="black", lty=2)
text(x=median(width(peak_set_signac)), y=10000, paste0("median peak size = ", median(width(peak_set_signac)), "bp"), pos = 4)
abline(v=mean(width(peak_set_signac)), col="blue", lty=2)
text(x=mean(width(peak_set_signac)), y=8000, paste0("mean peak size = ", round(mean(width(peak_set_signac)), digits=1), "bp"), pos = 4, col="blue")
hist(width(peak_set_signac), xlab="peak size", breaks=200, main="zoom", xlim=c(0,1000))
dev.off()


######################################
##### 2. ADD PEAK CALLING MATRIX #####
######################################

project_all <- addPeakMatrix(project_all)
getAvailableMatrices(project_all) # display which matrices are available: now there must be a "PeakMatrix"


################
##### SAVE #####
################

saveArchRProject(ArchRProj = project_all, outputDirectory = datadir, load = TRUE)


