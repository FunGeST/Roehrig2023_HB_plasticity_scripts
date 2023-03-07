
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
grn_dir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/4.Integration/2.Gene_regulatory_network/Final_GRN/H_LP_M"
resdir = datadir

# define thresholds for differential vs marker peaks study 
threshold_medium_logFC = log2(2)
threshold_medium_FDR = 1e-3

# Thresholds for peak2genes
FDR_threshold_peaksgene = 0.01 # for correlation > 0.4 all FDR are far below 0.01 anyway
varCutOffATAC= varCutOffRNA = 0.25 # default ArchR tfhreshold
cor_threshold=0.4


###########################
##### 0. DATA LOADING #####
###########################

# Load differential tables
diffPeaks_H_vs_LP <- geco.load(file.path(datadir, "diffPeaks_H_vs_LP.RData")) # LP is >0, H is <0
diffPeaks_H_vs_LP$peak <- paste(diffPeaks_H_vs_LP$seqnames, diffPeaks_H_vs_LP$start, diffPeaks_H_vs_LP$end, sep="_"); rownames(diffPeaks_H_vs_LP) <- diffPeaks_H_vs_LP$peak # add peak name
diffPeaks_H_HLP_LP_vs_M <- geco.load(file.path(datadir, "diffPeaks_H_HLP_LP_vs_M.RData")) # M is >0, H/H+LP/LP is <0
diffPeaks_H_HLP_LP_vs_M$peak <- paste(diffPeaks_H_HLP_LP_vs_M$seqnames, diffPeaks_H_HLP_LP_vs_M$start, diffPeaks_H_HLP_LP_vs_M$end, sep="_"); rownames(diffPeaks_H_HLP_LP_vs_M) <- diffPeaks_H_HLP_LP_vs_M$peak # add peak name

# Nearest gene
nearest_gene <- geco.load(file.path(archrdir, "dist_nearest_gene_per_peak.Rdata"))
rownames(nearest_gene) <- nearest_gene$peak

# Peaks linked to genes
peaks2genes <- geco.load(file.path(archrdir, "Peaks_to_genes/peaks2genes_all_correlations_df.RData"))
peaks2genes_sig <- peaks2genes[which(peaks2genes$FDR <= FDR_threshold_peaksgene & peaks2genes$VarQATAC >= varCutOffATAC & peaks2genes$VarQRNA >= varCutOffRNA & peaks2genes$Correlation >= cor_threshold),]
peaks2genes_sig <- peaks2genes_sig %>% group_by(peak) %>% summarise(gene=paste0(gene, collapse=", ")) %>% as.data.frame() # collapse when several genes associated to peak
peaks2genes_sig$peak <- gsub(":", "_", gsub("-", "_", peaks2genes_sig$peak)) # change peak nomenclature so that everything matches
rownames(peaks2genes_sig) <- peaks2genes_sig$peak

# GRN TF motifs in peaks
motifs_in_peaks <- geco.load(file.path(grn_dir, "top_TF_motifs_in_peaks_linked_to_genes.RData"))
motifs_in_peaks$peaks <- gsub(":", "_", gsub("-", "_", motifs_in_peaks$peaks)) # change peak nomenclature so that everything matches
rownames(motifs_in_peaks) <- motifs_in_peaks$peaks


##########################################################
##### CREATE ANNOTATION TABLE FOR DIFFERENTIAL PEAKS #####
##########################################################

annot_diff_peaks <- data.frame(peak=nearest_gene$peak, LP_vs_H_logFC=NA, LP_vs_H_FDR=NA, M_vs_H_LP_logFC=NA, M_vs_H_LP_FDR=NA, nearest_gene=NA, dist_nearest_TSS=NA, genes_linked=NA, TF_GRN=NA, stringsAsFactors = F)
rownames(annot_diff_peaks) <- annot_diff_peaks$peak

### Add peak info
annot_diff_peaks <- annot_diff_peaks %>% separate(col=peak, into=c("chr", "start", "end"), sep="_", remove = FALSE) %>% as.data.frame()
annot_diff_peaks$start <- as.numeric(annot_diff_peaks$start); annot_diff_peaks$end <- as.numeric(annot_diff_peaks$end) # convert coordinates to numeric

### Add differential info
annot_diff_peaks[intersect(annot_diff_peaks$peak, diffPeaks_H_vs_LP$peak), "LP_vs_H_logFC"] <- diffPeaks_H_vs_LP[intersect(annot_diff_peaks$peak, diffPeaks_H_vs_LP$peak), "Log2FC"]
annot_diff_peaks[intersect(annot_diff_peaks$peak, diffPeaks_H_vs_LP$peak), "LP_vs_H_FDR"] <- diffPeaks_H_vs_LP[intersect(annot_diff_peaks$peak, diffPeaks_H_vs_LP$peak), "FDR"]
annot_diff_peaks[intersect(annot_diff_peaks$peak, diffPeaks_H_HLP_LP_vs_M$peak), "M_vs_H_LP_logFC"] <- diffPeaks_H_HLP_LP_vs_M[intersect(annot_diff_peaks$peak, diffPeaks_H_HLP_LP_vs_M$peak), "Log2FC"]
annot_diff_peaks[intersect(annot_diff_peaks$peak, diffPeaks_H_HLP_LP_vs_M$peak), "M_vs_H_LP_FDR"] <- diffPeaks_H_HLP_LP_vs_M[intersect(annot_diff_peaks$peak, diffPeaks_H_HLP_LP_vs_M$peak), "FDR"]

### Nearest gene
annot_diff_peaks[intersect(annot_diff_peaks$peak, nearest_gene$peak), "nearest_gene"] <- nearest_gene[intersect(annot_diff_peaks$peak, nearest_gene$peak), "gene"]
annot_diff_peaks[intersect(annot_diff_peaks$peak, nearest_gene$peak), "dist_nearest_TSS"] <- nearest_gene[intersect(annot_diff_peaks$peak, nearest_gene$peak), "dist"]

### gene linkes to peaks
annot_diff_peaks[intersect(annot_diff_peaks$peak, peaks2genes_sig$peak), "genes_linked"] <- peaks2genes_sig[intersect(annot_diff_peaks$peak, peaks2genes_sig$peak), "gene"]

### associated TF (from GRN)
annot_diff_peaks[intersect(annot_diff_peaks$peak, motifs_in_peaks$peaks), "TF_GRN"] <- motifs_in_peaks[intersect(annot_diff_peaks$peak, motifs_in_peaks$peaks), "Motifs_in_peaks_peak2genelink"]

save(annot_diff_peaks, file=file.path(resdir, "Differential_peaks_annotations.RData"))
openxlsx::write.xlsx(annot_diff_peaks, file.path(resdir, "Differential_peaks_annotations.xlsx"))



 