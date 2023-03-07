
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

#####################
##### FUNCTIONS #####
#####################

enrich_chrom_states_table <- function(diffPeaks_pos_gr, diffPeaks_neg_gr, chrom_states_adult_liver, peak_set){
  # Look at the enrichment between all peaks and expected results for whole genome
  # We cannot use the PeakAnnoEnrichment function from AchR because it does enrichment on nb of peaks and assumes all peaks are the same size which is not our case here --> might bias results
  
  ##### 1. create recap table to contain results #####
  #---------------------------------------------------
  recap_cst18 <- data.frame(state= c("1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD", "5_Tx", "6_TxWk", "7_EnhG1", "8_EnhG2", "9_EnhA1", "10_EnhA2", "11_EnhWk", "12_ZNF/Rpts", "13_Het",
                                     "14_TssBiv", "15_EnhBiv", "16_ReprPC", "17_ReprPCWk", "18_Quies"), prop_diff_peaks_pos=NA, prop_diff_peaks_neg=NA, prop_all_peaks=NA, stringsAsFactors = F)
  rownames(recap_cst18) <- recap_cst18$state
  
  ##### 2. Proportion in differential peaks >0 #####
  #-------------------------------------------------
  # for each state compute nb of bp in peaks that intersect the state coordinates
  for(s in unique(chrom_states_adult_liver$state)){
    chrom_states_adult_liver_s <- chrom_states_adult_liver %>% filter(state==s)
    gr_state <- GenomicRanges::makeGRangesFromDataFrame(chrom_states_adult_liver_s, keep.extra.columns = F, ignore.strand = T, starts.in.df.are.0based	= T) # create Granges
    # get intersections in all peaks, pos diff peaks, neg diff peaks
    gr_intersect_all <- GenomicRanges::intersect(peak_set, gr_state); gr_intersect_pos <- GenomicRanges::intersect(diffPeaks_pos_gr[[1]], gr_state); gr_intersect_neg <- GenomicRanges::intersect(diffPeaks_neg_gr[[1]], gr_state) # get the intersected nucleotides between peaks and state coordinates 
    # if we do findOverlaps(peak_set, gr_state) we find that some peaks have several chromatin states linked
    recap_cst18[s,"prop_all_peaks"] <- sum(width(gr_intersect_all)) ; recap_cst18[s,"prop_diff_peaks_pos"] <- sum(width(gr_intersect_pos)) ; recap_cst18[s,"prop_diff_peaks_neg"] <- sum(width(gr_intersect_neg)) # how many bp in total in peaks intersect the state
  }
   #compute proportion
  recap_cst18$prop_all_peaks <- 100*recap_cst18$prop_all_peaks/sum(recap_cst18$prop_all_peaks) # divide by total nb of bp that are linked to one feature: since we consider bp and not peaks each bp is linked to 1 feature --> no overlap 
  recap_cst18$prop_diff_peaks_pos <- 100*recap_cst18$prop_diff_peaks_pos/sum(recap_cst18$prop_diff_peaks_pos) 
  recap_cst18$prop_diff_peaks_neg <- 100*recap_cst18$prop_diff_peaks_neg/sum(recap_cst18$prop_diff_peaks_neg) 
  
  ##### 3. Compute enrichment of diff peaks vs all peaks #####
  #-----------------------------------------------------------
  recap_cst18$enrichment_pos <- recap_cst18$prop_diff_peaks_pos / recap_cst18$prop_all_peaks; recap_cst18$enrichment_neg <- recap_cst18$prop_diff_peaks_neg / recap_cst18$prop_all_peaks
  return(recap_cst18) 
}

enrich_chrom_states_plots <- function(recap_cst18, resdir, comp, threshold){
  
  recap_cst18$state <- factor(recap_cst18$state, levels=unique(recap_cst18$state))
  
  pdf(file.path(resdir, comp, paste0("Diff_study_chromatin_states_", threshold, ".pdf"))) 
  # barplot
  print(ggplot(recap_cst18, aes(x=state, y=enrichment_pos, fill=state))+
          geom_bar(stat="identity")+
          theme_classic()+
          geom_hline(yintercept = 1, col="black", linetype="dashed")+
          ylab("Enrichment bp diff peaks >0 vs all peaks")+
          scale_fill_manual(values=c("darkred", "red"," red1", "red2", "forestgreen", "darkgreen", "darkolivegreen1", "darkolivegreen2", "orange", "orange1", "yellow2", "turquoise", "lightslateblue", "palevioletred4", "olivedrab3" ,"grey50", "grey80", "grey95")))
  print(ggplot(recap_cst18, aes(x=state, y=enrichment_neg, fill=state))+
          geom_bar(stat="identity")+
          theme_classic()+
          geom_hline(yintercept = 1, col="black", linetype="dashed")+
          ylab("Enrichment bp diff peaks <0 vs all peaks")+
          scale_fill_manual(values=c("darkred", "red"," red1", "red2", "forestgreen", "darkgreen", "darkolivegreen1", "darkolivegreen2", "orange", "orange1", "yellow2", "turquoise", "lightslateblue", "palevioletred4", "olivedrab3" ,"grey50", "grey80", "grey95")))
  # target plot (ggplot)
  print(ggplot(recap_cst18, aes(x=state, y=enrichment_pos, fill=state))+
          geom_bar(stat="identity")+
          theme_classic()+
          coord_polar()+
          geom_hline(yintercept = 1, col="black", linetype="dashed")+
          ylab("Enrichment bp diff peaks >0 vs all peaks")+
          theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
          scale_fill_manual(values=c("darkred", "red"," red1", "red2", "forestgreen", "darkgreen", "darkolivegreen1", "darkolivegreen2", "orange", "orange1", "yellow2", "turquoise", "lightslateblue", "palevioletred4", "olivedrab3" ,"grey50", "grey80", "grey95")))
  print(ggplot(recap_cst18, aes(x=state, y=enrichment_neg, fill=state))+
          geom_bar(stat="identity")+
          theme_classic()+
          coord_polar()+
          geom_hline(yintercept = 1, col="black", linetype="dashed")+
          ylab("Enrichment bp diff peaks <0 vs all peaks")+
          theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
          scale_fill_manual(values=c("darkred", "red"," red1", "red2", "forestgreen", "darkgreen", "darkolivegreen1", "darkolivegreen2", "orange", "orange1", "yellow2", "turquoise", "lightslateblue", "palevioletred4", "olivedrab3" ,"grey50", "grey80", "grey95")))
  dev.off()
  
}


######################
##### PARAMETERS #####
######################

archrdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/"
datadir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Differential_peaks/1.Tables"
resdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Peaks_study/Chromatin_states"

annot_dir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/0.Public_annotation"

# define thresholds for differential vs marker peaks study 
threshold_loose_logFC = log2(1.5)
threshold_loose_FDR = 0.01

threshold_medium_logFC = log2(2)
threshold_medium_FDR = 1e-3

threshold_strict_logFC = log2(3)
threshold_strict_FDR = 1e-12

# choose one comparison to study 
comp = c("NT_vs_T", "NT_vs_H", "NT_vs_H_LP", "NT_vs_LP", "NT_vs_M", "H_vs_LP", "H_HLP_LP_vs_M")[6] 


###########################
##### 0. DATA LOADING #####
###########################

# Load ArchR project
project_all <- readRDS(file.path(archrdir, "Save-ArchR-Project.rds"))
peak_set <- getPeakSet(project_all)

# Load differential tables
if(comp=="NT_vs_T"){diffPeaks_SE <- geco.load(file.path(datadir, "diffPeaks_NT_vs_T_summarized_exp.RData")) # T is >0, NT is <0
}else if(comp=="NT_vs_H"){diffPeaks_SE <- geco.load(file.path(datadir, "diffPeaks_NT_vs_H_summarized_exp.RData")) # LP is >0, H is <0
}else if(comp=="NT_vs_H_LP"){diffPeaks_SE <- geco.load(file.path(datadir, "diffPeaks_NT_vs_H_LP_summarized_exp.RData")) # LP is >0, H is <0
}else if(comp=="NT_vs_LP"){diffPeaks_SE <- geco.load(file.path(datadir, "diffPeaks_NT_vs_LP_summarized_exp.RData")) # LP is >0, H is <0
}else if(comp=="NT_vs_M"){diffPeaks_SE <- geco.load(file.path(datadir, "diffPeaks_NT_vs_M_summarized_exp.RData")) # LP is >0, H is <0
}else if(comp=="H_vs_LP"){diffPeaks_SE <- geco.load(file.path(datadir, "diffPeaks_H_vs_LP_summarized_exp.RData")) # LP is >0, H is <0
}else if(comp=="H_HLP_LP_vs_M"){diffPeaks_SE <- geco.load(file.path(datadir, "diffPeaks_H_HLP_LP_vs_M_summarized_exp.RData"))} # M is >0, H/H+LP/LP is <0

# Load annotation for peak enrichment
chrom_states_adult_liver <- read.table(file.path(annot_dir, "E066_18_core_K27ac_mnemonics.bed"), sep="\t", header=F); colnames(chrom_states_adult_liver) <- c("chr", "start", "end", "state") # E066 = normal adult liver 
chrom_states_adult_liver <- chrom_states_adult_liver %>% filter(chr %in% paste0("chr", 1:22)) # restrict to chromosomes 1 to 22 = like in ArchR
# chrom_states_HepG2 <- read.table(file.path(annot_dir, "E118_18_core_K27ac_mnemonics.bed"), sep="\t", header=F); colnames(chrom_states_HepG2) <- c("chr", "start", "end", "state")
# here we only use the adult liver data because cell line HepG2 brings too many biases 



###############################
##### 1. GENOMIC LOCATION #####
###############################

for(threshold in c("loose_thresholds", "medium_thresholds", "strict_thresholds")){
  if(threshold == "loose_thresholds"){FDR_threshold <- threshold_loose_FDR; Log2FC_threshold <- threshold_loose_logFC}
  if(threshold == "medium_thresholds"){FDR_threshold <- threshold_medium_FDR; Log2FC_threshold <- threshold_medium_logFC}
  if(threshold == "strict_thresholds"){FDR_threshold <- threshold_strict_FDR; Log2FC_threshold <- threshold_strict_logFC}
  
# retrieve differential peaks
  
  diffPeaks_pos_gr <- getMarkers(diffPeaks_SE, cutOff = paste0("FDR <= ", FDR_threshold, " & Log2FC >= ", Log2FC_threshold), returnGR = TRUE)
  diffPeaks_neg_gr <- getMarkers(diffPeaks_SE, cutOff = paste0("FDR <= ", FDR_threshold, " & Log2FC <= ", -Log2FC_threshold), returnGR = TRUE)

  # perform enrichment function
  
  recap_cst18 <- enrich_chrom_states_table(diffPeaks_pos_gr, diffPeaks_neg_gr, chrom_states_adult_liver, peak_set)
  save(recap_cst18, file=file.path(resdir, comp, paste0("Diff_study_chromatin_states_", threshold, ".RData")))
  enrich_chrom_states_plots(recap_cst18, resdir, comp, threshold)
  
}
  
  
 