
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

archrdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/"
datadir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Differential_peaks/1.Tables"
resdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Peaks_study"

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
project_all <- readRDS(file.path(archrdir, "Save-ArchR-Project.rds"))

# Load differential tables
diffPeaks_NT_vs_T <- geco.load(file.path(datadir, "diffPeaks_NT_vs_T.RData")) # T is >0, NT is <0
diffPeaks_NT_vs_H <- geco.load(file.path(datadir, "diffPeaks_NT_vs_H.RData")) # H is >0, NT is <0
diffPeaks_NT_vs_H_LP <- geco.load(file.path(datadir, "diffPeaks_NT_vs_H_LP.RData")) # H+LP is >0, NT is <0
diffPeaks_NT_vs_LP <- geco.load(file.path(datadir, "diffPeaks_NT_vs_LP.RData")) # LP is >0, NT is <0
diffPeaks_NT_vs_M <- geco.load(file.path(datadir, "diffPeaks_NT_vs_M.RData")) # M is >0, NT is <0
diffPeaks_H_vs_LP <- geco.load(file.path(datadir, "diffPeaks_H_vs_LP.RData")) # LP is >0, H is <0
diffPeaks_H_HLP_LP_vs_M <- geco.load(file.path(datadir, "diffPeaks_H_HLP_LP_vs_M.RData")) # M is >0, H/H+LP/LP is <0


###########################################
##### 1. NUMBER OF DIFFERENTIAL PEAKS #####
###########################################

colors_bars <- rep(c("#CD2626", "#1873CD"), times=7)

for(threshold in c("loose_thresholds", "medium_thresholds", "strict_thresholds")){
  if(threshold == "loose_thresholds"){FDR_threshold <- threshold_loose_FDR; Log2FC_threshold <- threshold_loose_logFC}
  if(threshold == "medium_thresholds"){FDR_threshold <- threshold_medium_FDR; Log2FC_threshold <- threshold_medium_logFC}
  if(threshold == "strict_thresholds"){FDR_threshold <- threshold_strict_FDR; Log2FC_threshold <- threshold_strict_logFC}
  
  
  ##### 1. Create table containing nb of differential peaks for each comparison #####
  #----------------------------------------------------------------------------------
  
  df_nb_diff_peaks <- as.data.frame(matrix(NA, nrow=2, ncol=7))
  rownames(df_nb_diff_peaks) <- c("pos_diff_peaks", "neg_diff_peaks")
  colnames(df_nb_diff_peaks) <- c("NT_vs_T", "NT_vs_H", "NT_vs_H_LP", "NT_vs_LP", "NT_vs_M", "H_vs_LP", "H_HLP_LP_vs_M")
  df_nb_diff_peaks$pos_or_neg <- rownames(df_nb_diff_peaks)
  
  df_nb_diff_peaks["pos_diff_peaks","NT_vs_T"] <- length(which(diffPeaks_NT_vs_T$FDR <= FDR_threshold & diffPeaks_NT_vs_T$Log2FC >= Log2FC_threshold))
  df_nb_diff_peaks["neg_diff_peaks","NT_vs_T"] <- -length(which(diffPeaks_NT_vs_T$FDR <= FDR_threshold & diffPeaks_NT_vs_T$Log2FC <= -Log2FC_threshold))
  df_nb_diff_peaks["pos_diff_peaks","NT_vs_H"] <- length(which(diffPeaks_NT_vs_H$FDR <= FDR_threshold & diffPeaks_NT_vs_H$Log2FC >= Log2FC_threshold))
  df_nb_diff_peaks["neg_diff_peaks","NT_vs_H"] <- -length(which(diffPeaks_NT_vs_H$FDR <= FDR_threshold & diffPeaks_NT_vs_H$Log2FC <= -Log2FC_threshold))
  df_nb_diff_peaks["pos_diff_peaks","NT_vs_H_LP"] <- length(which(diffPeaks_NT_vs_H_LP$FDR <= FDR_threshold & diffPeaks_NT_vs_H_LP$Log2FC >= Log2FC_threshold))
  df_nb_diff_peaks["neg_diff_peaks","NT_vs_H_LP"] <- -length(which(diffPeaks_NT_vs_H_LP$FDR <= FDR_threshold & diffPeaks_NT_vs_H_LP$Log2FC <= -Log2FC_threshold))
  df_nb_diff_peaks["pos_diff_peaks","NT_vs_LP"] <- length(which(diffPeaks_NT_vs_LP$FDR <= FDR_threshold & diffPeaks_NT_vs_LP$Log2FC >= Log2FC_threshold))
  df_nb_diff_peaks["neg_diff_peaks","NT_vs_LP"] <- -length(which(diffPeaks_NT_vs_LP$FDR <= FDR_threshold & diffPeaks_NT_vs_LP$Log2FC <= -Log2FC_threshold))
  df_nb_diff_peaks["pos_diff_peaks","NT_vs_M"] <- length(which(diffPeaks_NT_vs_M$FDR <= FDR_threshold & diffPeaks_NT_vs_M$Log2FC >= Log2FC_threshold))
  df_nb_diff_peaks["neg_diff_peaks","NT_vs_M"] <- -length(which(diffPeaks_NT_vs_M$FDR <= FDR_threshold & diffPeaks_NT_vs_M$Log2FC <= -Log2FC_threshold))
  df_nb_diff_peaks["pos_diff_peaks","H_vs_LP"] <- length(which(diffPeaks_H_vs_LP$FDR <= FDR_threshold & diffPeaks_H_vs_LP$Log2FC >= Log2FC_threshold))
  df_nb_diff_peaks["neg_diff_peaks","H_vs_LP"] <- -length(which(diffPeaks_H_vs_LP$FDR <= FDR_threshold & diffPeaks_H_vs_LP$Log2FC <= -Log2FC_threshold))
  df_nb_diff_peaks["pos_diff_peaks","H_HLP_LP_vs_M"] <- length(which(diffPeaks_H_HLP_LP_vs_M$FDR <= FDR_threshold & diffPeaks_H_HLP_LP_vs_M$Log2FC >= Log2FC_threshold))
  df_nb_diff_peaks["neg_diff_peaks","H_HLP_LP_vs_M"] <- -length(which(diffPeaks_H_HLP_LP_vs_M$FDR <= FDR_threshold & diffPeaks_H_HLP_LP_vs_M$Log2FC <= -Log2FC_threshold))
  # put minus before nb of negative peaks to center the barplot on 0 
  
  df_nb_diff_peaks2 <- reshape2::melt(df_nb_diff_peaks) # prepare table for barplot
  df_nb_diff_peaks2$variable_pos_or_neg <- factor(paste(df_nb_diff_peaks2$variable, df_nb_diff_peaks2$pos_or_neg, sep="_"), levels=c(paste(df_nb_diff_peaks2$variable, df_nb_diff_peaks2$pos_or_neg, sep="_")))
  
  ##### 2. Generate barplots #####
  #-------------------------------
  
  pdf(file.path(resdir, paste0("Nb_of_diff_peaks_", threshold, ".pdf")))
  print(ggplot(df_nb_diff_peaks2, aes(x=variable, y=value, fill=variable_pos_or_neg))+
    geom_bar(position="stack", stat="identity")+
    theme_classic()+
    scale_fill_manual(values=colors_bars)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  save(df_nb_diff_peaks, file=file.path(resdir, paste0("Nb_of_diff_peaks_", threshold, ".RData")))

}
