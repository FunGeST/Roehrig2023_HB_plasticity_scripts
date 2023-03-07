

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
resdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Differential_peaks/"
rna_dir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/"

# define thresholds for differential vs marker peaks study 
threshold_loose_logFC = log2(1.5)
threshold_loose_FDR = 0.01

threshold_medium_logFC = log2(2)
threshold_medium_FDR = 1e-3

threshold_strict_logFC = log2(3)
threshold_strict_FDR = 1e-12


###########################
##### 0. DATA LOADING #####
###########################

# Load ArchR project
project_all <- readRDS(file.path(datadir, "Save-ArchR-Project.rds"))
atac_metadata <- getCellColData(project_all) # extract ArchR cell metadata

# Load RNA subtype annotations
seurat_sample_all <- geco.load(file.path(rna_dir, "Seurat_object_analyses.RData"))
rna_metadata <- seurat_sample_all@meta.data
rownames(rna_metadata) <- paste0(gsub("_", "#", rownames(rna_metadata)), "-1") # replace by ArchR nomenclature

# Add RNA subtypes info to ATAC ArchR project
atac_metadata$RNA_subtypes <- NA
atac_metadata[intersect(rownames(atac_metadata), rownames(rna_metadata)),"RNA_subtypes"] <- as.character(rna_metadata[intersect(rownames(atac_metadata), rownames(rna_metadata)),"H_LP_M_clean_groups_1"]) # 10939 common cells between ArchR (all cells) and Seurat (tumor cells)
project_all$RNA_subtypes <- atac_metadata$RNA_subtypes # add to ArchR project

# Define hepatocytes from ATAC data
atac_metadata$ATAC_hepatocytes <- NA
atac_metadata[which(atac_metadata$Sample %in% c("CHC2959N", "CHC3377N") & atac_metadata$Clusters %in% c("C9", "C10")),"ATAC_hepatocytes"] <- "hepatocytes"
project_all$ATAC_hepatocytes <- atac_metadata$ATAC_hepatocytes

atac_metadata$tumor_normal <- NA
atac_metadata[which(atac_metadata$ATAC_hepatocytes == "hepatocytes"),"tumor_normal"] <- "normal"
atac_metadata[which(!is.na(atac_metadata$RNA_subtypes)),"tumor_normal"] <- "tumor"
project_all$tumor_normal <- atac_metadata$tumor_normal

# Create column that combines all cell info (in order to remove NAs for marker peaks analysis)
atac_metadata$RNA_subtypes_hepatocytes <- atac_metadata$RNA_subtypes
atac_metadata[which(atac_metadata$ATAC_hepatocytes == "hepatocytes"),"RNA_subtypes_hepatocytes"] <- "normal"
atac_metadata[which(is.na(atac_metadata$RNA_subtypes_hepatocytes)),"RNA_subtypes_hepatocytes"] <- "other"
project_all$RNA_subtypes_hepatocytes <- atac_metadata$RNA_subtypes_hepatocytes

# 1518 hepatocytes (ATAC), 10939 tumor cells (RNA X ATAC)
# 2526 H tumor cells, 4013 H+LP tumor cells, 2436 LP tumor cells, 1964 M tumor cells (RNA X ATAC)


############################################
##### 1. COMPUTE DIFFERENTIAL ANALYSIS #####
############################################

##### WARNING: HERE FOR DIFFERENTIAL PEAKS WE CHOOSE TO KEEP CELLS WITH INTERSECTION ATAC (ARCHR) X RNA (SIGNAC) INSTEAD OF TAKING ALL ATAC CELLS FOR CONSISTENCY
# EXCEPT FOR HEPATOCYTES WHICH ARE CLUSTERS C9 (mostly CHC3377N) AND C10 (mostly CHC2959N) FROM ATAC DATA TO HAVE MORE CELLS

# in differential analyses ArchR matches the cells to study with background cells by reducing bias (TSS enrichment, log10(frags): https://github.com/GreenleafLab/ArchR/issues/470)

##### 1. Hepatocytes vs tumor cells #####
#----------------------------------------

diffPeaks_NT_vs_T_SE <- getMarkerFeatures(
  ArchRProj = project_all, 
  useMatrix = "PeakMatrix", 
  groupBy = "tumor_normal",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "tumor",
  bgdGroups = "normal",
  maxCells=1000 # choose high number of max cells to have as many cells as possible considered in differential analysis
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

save(diffPeaks_NT_vs_T_SE, file=file.path(resdir, "1.Tables/diffPeaks_NT_vs_T_summarized_exp.RData")) 

# Convert to a dataframe and remove the NA values for FDR (there can be few % of NA FDR because mean=0 and meanBGD=0 too)

diffPeaks_NT_vs_T <- getMarkers(diffPeaks_NT_vs_T_SE, cutOff = "FDR < 10 ") # choose unreachable thresholds to keep all peaks
diffPeaks_NT_vs_T <- as.data.frame(diffPeaks_NT_vs_T[[1]]) # initially diffPeaks_NT_vs_T is a list
# is sorted by FDR
save(diffPeaks_NT_vs_T, file=file.path(resdir, "1.Tables/diffPeaks_NT_vs_T.RData")) 


##### 2. H tumor cells vs LP tumor cells #####
#---------------------------------------------

diffPeaks_H_vs_LP_SE <- getMarkerFeatures(
  ArchRProj = project_all, 
  useMatrix = "PeakMatrix", 
  groupBy = "RNA_subtypes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "LP",
  bgdGroups = "H",
  maxCells=1000 # choose high number of max cells to have as many cells as possible considered in differential analysis
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

save(diffPeaks_H_vs_LP_SE, file=file.path(resdir, "1.Tables/diffPeaks_H_vs_LP_summarized_exp.RData")) 

# Convert to a dataframe and remove the NA values for FDR (there can be few % of NA FDR because mean=0 and meanBGD=0 too)

diffPeaks_H_vs_LP <- getMarkers(diffPeaks_H_vs_LP_SE, cutOff = "FDR < 10 ") # choose unreachable thresholds to keep all peaks
diffPeaks_H_vs_LP <- as.data.frame(diffPeaks_H_vs_LP[[1]]) # initially diffPeaks_NT_vs_T is a list
# is sorted by FDR
save(diffPeaks_H_vs_LP, file=file.path(resdir, "1.Tables/diffPeaks_H_vs_LP.RData")) 


##### 3. H/H+LP/LP tumor cells vs M tumor cells #####
#----------------------------------------------------

diffPeaks_H_HLP_LP_vs_M_SE <- getMarkerFeatures(
  ArchRProj = project_all, 
  useMatrix = "PeakMatrix", 
  groupBy = "RNA_subtypes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "M",
  bgdGroups = c("H", "H+LP", "LP"),
  maxCells=1000 # choose high number of max cells to have as many cells as possible considered in differential analysis
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

save(diffPeaks_H_HLP_LP_vs_M_SE, file=file.path(resdir, "1.Tables/diffPeaks_H_HLP_LP_vs_M_summarized_exp.RData")) 

# Convert to a dataframe and remove the NA values for FDR (there can be few % of NA FDR because mean=0 and meanBGD=0 too)

diffPeaks_H_HLP_LP_vs_M <- getMarkers(diffPeaks_H_HLP_LP_vs_M_SE, cutOff = "FDR < 10 ") # choose unreachable thresholds to keep all peaks
diffPeaks_H_HLP_LP_vs_M <- as.data.frame(diffPeaks_H_HLP_LP_vs_M[[1]]) # initially diffPeaks_NT_vs_T is a list
# is sorted by FDR
save(diffPeaks_H_HLP_LP_vs_M, file=file.path(resdir, "1.Tables/diffPeaks_H_HLP_LP_vs_M.RData")) 


##### 4. Hepatocytes vs H tumor cells #####
#------------------------------------------

diffPeaks_NT_vs_H_SE <- getMarkerFeatures(
  ArchRProj = project_all, 
  useMatrix = "PeakMatrix", 
  groupBy = "RNA_subtypes_hepatocytes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "H",
  bgdGroups = "normal",
  maxCells=1000 # choose high number of max cells to have as many cells as possible considered in differential analysis
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

save(diffPeaks_NT_vs_H_SE, file=file.path(resdir, "1.Tables/diffPeaks_NT_vs_H_summarized_exp.RData")) 

# Convert to a dataframe and remove the NA values for FDR (there can be few % of NA FDR because mean=0 and meanBGD=0 too)

diffPeaks_NT_vs_H <- getMarkers(diffPeaks_NT_vs_H_SE, cutOff = "FDR < 10 ") # choose unreachable thresholds to keep all peaks
diffPeaks_NT_vs_H <- as.data.frame(diffPeaks_NT_vs_H[[1]]) # initially diffPeaks_NT_vs_T is a list
# is sorted by FDR
save(diffPeaks_NT_vs_H, file=file.path(resdir, "1.Tables/diffPeaks_NT_vs_H.RData")) 


##### 5. Hepatocytes vs H+LP tumor cells #####
#---------------------------------------------

diffPeaks_NT_vs_H_LP_SE <- getMarkerFeatures(
  ArchRProj = project_all, 
  useMatrix = "PeakMatrix", 
  groupBy = "RNA_subtypes_hepatocytes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "H+LP",
  bgdGroups = "normal",
  maxCells=1000 # choose high number of max cells to have as many cells as possible considered in differential analysis
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

save(diffPeaks_NT_vs_H_LP_SE, file=file.path(resdir, "1.Tables/diffPeaks_NT_vs_H_LP_summarized_exp.RData")) 

# Convert to a dataframe and remove the NA values for FDR (there can be few % of NA FDR because mean=0 and meanBGD=0 too)

diffPeaks_NT_vs_H_LP <- getMarkers(diffPeaks_NT_vs_H_LP_SE, cutOff = "FDR < 10 ") # choose unreachable thresholds to keep all peaks
diffPeaks_NT_vs_H_LP <- as.data.frame(diffPeaks_NT_vs_H_LP[[1]]) # initially diffPeaks_NT_vs_T is a list
# is sorted by FDR
save(diffPeaks_NT_vs_H_LP, file=file.path(resdir, "1.Tables/diffPeaks_NT_vs_H_LP.RData")) 


##### 6. Hepatocytes vs LP tumor cells #####
#-------------------------------------------

diffPeaks_NT_vs_LP_SE <- getMarkerFeatures(
  ArchRProj = project_all, 
  useMatrix = "PeakMatrix", 
  groupBy = "RNA_subtypes_hepatocytes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "LP",
  bgdGroups = "normal",
  maxCells=1000 # choose high number of max cells to have as many cells as possible considered in differential analysis
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

save(diffPeaks_NT_vs_LP_SE, file=file.path(resdir, "1.Tables/diffPeaks_NT_vs_LP_summarized_exp.RData")) 

# Convert to a dataframe and remove the NA values for FDR (there can be few % of NA FDR because mean=0 and meanBGD=0 too)

diffPeaks_NT_vs_LP <- getMarkers(diffPeaks_NT_vs_LP_SE, cutOff = "FDR < 10 ") # choose unreachable thresholds to keep all peaks
diffPeaks_NT_vs_LP <- as.data.frame(diffPeaks_NT_vs_LP[[1]]) # initially diffPeaks_NT_vs_T is a list
# is sorted by FDR
save(diffPeaks_NT_vs_LP, file=file.path(resdir, "1.Tables/diffPeaks_NT_vs_LP.RData")) 


##### 7. Hepatocytes vs M tumor cells #####
#----------------------------------------

diffPeaks_NT_vs_M_SE <- getMarkerFeatures(
  ArchRProj = project_all, 
  useMatrix = "PeakMatrix", 
  groupBy = "RNA_subtypes_hepatocytes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "M",
  bgdGroups = "normal",
  maxCells=1000 # choose high number of max cells to have as many cells as possible considered in differential analysis
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

save(diffPeaks_NT_vs_M_SE, file=file.path(resdir, "1.Tables/diffPeaks_NT_vs_M_summarized_exp.RData")) 

# Convert to a dataframe and remove the NA values for FDR (there can be few % of NA FDR because mean=0 and meanBGD=0 too)

diffPeaks_NT_vs_M <- getMarkers(diffPeaks_NT_vs_M_SE, cutOff = "FDR < 10 ") # choose unreachable thresholds to keep all peaks
diffPeaks_NT_vs_M <- as.data.frame(diffPeaks_NT_vs_M[[1]]) # initially diffPeaks_NT_vs_T is a list
# is sorted by FDR
save(diffPeaks_NT_vs_M, file=file.path(resdir, "1.Tables/diffPeaks_NT_vs_M.RData")) 



##### OPTIONAL #####

##########################################
##### 2. IDENTIFY DIFFERENTIAL PEAKS #####
##########################################

# The idea is to identify several thresholds to see during analyses if there are striking differences

##### 1. Plot metrics to define threhsolds #####
#-----------------------------------------------

pdf(file.path(resdir, "2.Plots/", paste0("Volcano_plots_loose_thresholds_FDR_", threshold_loose_FDR, "_log2FC_", threshold_loose_logFC, ".pdf")))
### Mean/Log2fc plot
markerPlot(seMarker = diffPeaks_H_vs_LP_SE, name = "LP", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "MA")
### Volcano plot
markerPlot(seMarker = diffPeaks_H_vs_LP_SE, name = "LP", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "Volcano")
### Mean/Log2fc plot
markerPlot(seMarker = diffPeaks_H_HLP_LP_vs_M_SE, name = "M", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "MA")
### Volcano plot
markerPlot(seMarker = diffPeaks_H_HLP_LP_vs_M_SE, name = "M", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "Volcano")
### Mean/Log2fc plot
markerPlot(seMarker = diffPeaks_NT_vs_T_SE, name = "tumor", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "MA")
### Volcano plot
markerPlot(seMarker = diffPeaks_NT_vs_T_SE, name = "tumor", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "Volcano")
dev.off()

pdf(file.path(resdir, "2.Plots/", paste0("Volcano_plots_medium_thresholds_FDR_", threshold_medium_FDR, "_log2FC_", threshold_medium_logFC, ".pdf")))
### Mean/Log2fc plot
markerPlot(seMarker = diffPeaks_H_vs_LP_SE, name = "LP", cutOff = paste0("FDR <= ",threshold_medium_FDR, " & abs(Log2FC) >= ", threshold_medium_logFC), plotAs = "MA")
### Volcano plot
markerPlot(seMarker = diffPeaks_H_vs_LP_SE, name = "LP", cutOff = paste0("FDR <= ",threshold_medium_FDR, " & abs(Log2FC) >= ", threshold_medium_logFC), plotAs = "Volcano")
### Mean/Log2fc plot
markerPlot(seMarker = diffPeaks_H_HLP_LP_vs_M_SE, name = "M", cutOff = paste0("FDR <= ",threshold_medium_FDR, " & abs(Log2FC) >= ", threshold_medium_logFC), plotAs = "MA")
### Volcano plot
markerPlot(seMarker = diffPeaks_H_HLP_LP_vs_M_SE, name = "M", cutOff = paste0("FDR <= ",threshold_medium_FDR, " & abs(Log2FC) >= ", threshold_medium_logFC), plotAs = "Volcano")
### Mean/Log2fc plot
markerPlot(seMarker = diffPeaks_NT_vs_T_SE, name = "tumor", cutOff = paste0("FDR <= ",threshold_medium_FDR, " & abs(Log2FC) >= ", threshold_medium_logFC), plotAs = "MA")
### Volcano plot
markerPlot(seMarker = diffPeaks_NT_vs_T_SE, name = "tumor", cutOff = paste0("FDR <= ",threshold_medium_FDR, " & abs(Log2FC) >= ", threshold_medium_logFC), plotAs = "Volcano")
dev.off()

pdf(file.path(resdir, "2.Plots/", paste0("Volcano_plots_strict_thresholds_FDR_", threshold_strict_FDR, "_log2FC_", threshold_strict_logFC, ".pdf")))
### Mean/Log2fc plot
markerPlot(seMarker = diffPeaks_H_vs_LP_SE, name = "LP", cutOff = paste0("FDR <= ",threshold_strict_FDR, " & abs(Log2FC) >= ", threshold_strict_logFC), plotAs = "MA")
### Volcano plot
markerPlot(seMarker = diffPeaks_H_vs_LP_SE, name = "LP", cutOff = paste0("FDR <= ",threshold_strict_FDR, " & abs(Log2FC) >= ", threshold_strict_logFC), plotAs = "Volcano")
### Mean/Log2fc plot
markerPlot(seMarker = diffPeaks_H_HLP_LP_vs_M_SE, name = "M", cutOff = paste0("FDR <= ",threshold_strict_FDR, " & abs(Log2FC) >= ", threshold_strict_logFC), plotAs = "MA")
### Volcano plot
markerPlot(seMarker = diffPeaks_H_HLP_LP_vs_M_SE, name = "M", cutOff = paste0("FDR <= ",threshold_strict_FDR, " & abs(Log2FC) >= ", threshold_strict_logFC), plotAs = "Volcano")
### Mean/Log2fc plot
markerPlot(seMarker = diffPeaks_NT_vs_T_SE, name = "tumor", cutOff = paste0("FDR <= ",threshold_strict_FDR, " & abs(Log2FC) >= ", threshold_strict_logFC), plotAs = "MA")
### Volcano plot
markerPlot(seMarker = diffPeaks_NT_vs_T_SE, name = "tumor", cutOff = paste0("FDR <= ",threshold_strict_FDR, " & abs(Log2FC) >= ", threshold_strict_logFC), plotAs = "Volcano")
dev.off()


################
##### SAVE #####
################

#saveArchRProject(ArchRProj = project_all, outputDirectory = datadir, load = TRUE)

