
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(Signac)
library(Seurat)
library(future)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
set.seed(1234)

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

name_samples <- c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3377T", "CHC3610T", "CHC3662T")

macs2_path <- "/home/amelie/miniconda2/envs/r_env_4.1/bin/macs2"
resdir <- "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Signac/2.Peak_calling"

common_peaks_definition = c("reduce", "disjoin")[1] # when we merge the peaks from the different samples, do we want to merge overlapping peaks (reduce) or separate them based on the overlapping peak ends (disjoin)?
# why choose "disjoin"? because with "reduce" we might over-merge groups of peaks where there would be for example 2 real peaks in some subtypes and 3 in others --> "disjoin" enables to keep these subtilities  
# why choose "reduce"? because with "disjoin" the number of very small peaks increases and sometimes they have no relevance (ex: large peak scattered in 4-5 small peaks)


###########################
##### 0. DATA LOADING #####
###########################

peaks_CHC2959T <- geco.load(file.path(resdir, "CHC2959T/peaks_granges_per_subtype.RData"))
peaks_CHC2960T <- geco.load(file.path(resdir, "CHC2960T/peaks_granges_per_subtype.RData"))
peaks_CHC3133T <- geco.load(file.path(resdir, "CHC3133T/peaks_granges_per_subtype.RData"))
peaks_CHC3377T <- geco.load(file.path(resdir, "CHC3377T/peaks_granges_per_subtype.RData"))
peaks_CHC3610T <- geco.load(file.path(resdir, "CHC3610T/peaks_granges_per_subtype.RData"))
peaks_CHC3662T <- geco.load(file.path(resdir, "CHC3662T/peaks_granges_per_subtype.RData"))

###  retrieve gene annotation

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # v86 = hg38 
# seqlevelsStyle(annotation) <- "UCSC" # this function does not work??
seqlevels(annotation) <- paste0("chr", seqlevels(annotation)) # to replace the line above


##########################################
##### 1. CREATION OF COMMON PEAK SET #####
##########################################

if(common_peaks_definition=="reduce"){
  combined.peaks <- reduce(x = c(peaks_CHC2959T, peaks_CHC2960T, peaks_CHC3133T, peaks_CHC3377T, peaks_CHC3610T, peaks_CHC3662T)) 
}else if(common_peaks_definition=="disjoin"){
  combined.peaks <- disjoin(x = c(peaks_CHC2959T, peaks_CHC2960T, peaks_CHC3133T, peaks_CHC3377T, peaks_CHC3610T, peaks_CHC3662T)) 
}

# Filter out bad peaks based on length

peakwidths <- width(combined.peaks) # determine widths of the peaks
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20] # keep peaks with > 20 bp and < 10000 bp but seemingly there are not such peaks in our dataset (min=200bp, max=3824 bp)
combined.peaks # 158707 peaks in total

save(combined.peaks, file=file.path(resdir, "combined_peaks_all_samples.RData"))

# convert granges to bed

combined.peaks_df <- data.frame(chr=seqnames(combined.peaks),
                 start=start(combined.peaks),
                 end=end(combined.peaks))
combined.peaks_df$name <- paste(combined.peaks_df$chr, combined.peaks_df$start, combined.peaks_df$end, sep="_")
combined.peaks_df$score=1

write.table(combined.peaks_df, file=file.path(resdir, "combined_peaks_all_samples.bed"), quote=F, sep="\t", row.names=F, col.names=F)



