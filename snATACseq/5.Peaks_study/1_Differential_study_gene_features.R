
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

enrich_gene_features_table <- function(diffPeaks_pos_gr, diffPeaks_neg_gr, promoter_body, peak_set){
  # Look at the enrichment between all peaks and expected results for whole genome
  # We cannot use the PeakAnnoEnrichment function from AchR because it does enrichment on nb of peaks and assumes all peaks are the same size which is not our case here --> might bias results
  
  ##### 1. create recap table to contain results #####
  #---------------------------------------------------
    recap_gene <- data.frame(feature= c("promoter", "body", "out"), prop_diff_peaks_pos=NA, prop_diff_peaks_neg=NA, prop_all_peaks=NA, stringsAsFactors = F); rownames(recap_gene) <- recap_gene$feature
  
  ##### 2. Proportion in differential peaks >0 #####
  #-------------------------------------------------
  # for each state compute nb of bp in peaks that intersect the state coordinates
  for(f in unique(promoter_body$gene_feature)){ # only for promoter and body
    promoter_body_f <- promoter_body %>% filter(gene_feature==f)
    gr_feature <- GenomicRanges::makeGRangesFromDataFrame(promoter_body_f, keep.extra.columns = F, ignore.strand = T, starts.in.df.are.0based	= T) # create Granges
    # get intersections in all peaks, pos diff peaks, neg diff peaks
    gr_intersect_all <- GenomicRanges::intersect(peak_set, gr_feature); gr_intersect_pos <- GenomicRanges::intersect(diffPeaks_pos_gr[[1]], gr_feature); gr_intersect_neg <- GenomicRanges::intersect(diffPeaks_neg_gr[[1]], gr_feature) # get the intersected nucleotides between peaks and state coordinates 
    # if we do findOverlaps(peak_set, gr_state) we find that some peaks have several chromatin states linked
    recap_gene[f,"prop_all_peaks"] <- sum(width(gr_intersect_all)) ; recap_gene[f,"prop_diff_peaks_pos"] <- sum(width(gr_intersect_pos)) ; recap_gene[f,"prop_diff_peaks_neg"] <- sum(width(gr_intersect_neg)) # how many bp in total in peaks intersect the state
  }
  # add out of gene
  recap_gene["out","prop_all_peaks"] <- sum(width(peak_set)) - recap_gene["promoter","prop_all_peaks"] - recap_gene["body","prop_all_peaks"] # the remaining nucleotides in peaks are out of genes
  recap_gene["out","prop_diff_peaks_pos"] <- sum(width(diffPeaks_pos_gr[[1]])) - recap_gene["promoter","prop_diff_peaks_pos"] - recap_gene["body","prop_diff_peaks_pos"] 
  recap_gene["out","prop_diff_peaks_neg"] <- sum(width(diffPeaks_neg_gr[[1]])) - recap_gene["promoter","prop_diff_peaks_neg"] - recap_gene["body","prop_diff_peaks_neg"] 
  #compute proportion
  recap_gene$prop_all_peaks <- 100*recap_gene$prop_all_peaks/sum(recap_gene$prop_all_peaks) # divide by total nb of bp that are linked to one feature: since we consider bp and not peaks each bp is linked to 1 feature --> no overlap 
  recap_gene$prop_diff_peaks_pos <- 100*recap_gene$prop_diff_peaks_pos/sum(recap_gene$prop_diff_peaks_pos) 
  recap_gene$prop_diff_peaks_neg <- 100*recap_gene$prop_diff_peaks_neg/sum(recap_gene$prop_diff_peaks_neg) 
  
  ##### 3. Compute enrichment of diff peaks vs all peaks #####
  #-----------------------------------------------------------
  recap_gene$enrichment_pos <- recap_gene$prop_diff_peaks_pos / recap_gene$prop_all_peaks; recap_gene$enrichment_neg <- recap_gene$prop_diff_peaks_neg / recap_gene$prop_all_peaks
  return(recap_gene) 
}

enrich_gene_features_plots <- function(recap_gene, resdir, comp, threshold){
  
  recap_gene$feature <- factor(recap_gene$feature, levels=unique(recap_gene$feature))
  
  pdf(file.path(resdir, comp, paste0("Diff_study_gene_features_", threshold, ".pdf"))) 
  # barplot
  print(ggplot(recap_gene, aes(x=feature, y=enrichment_pos, fill=feature))+
          geom_bar(stat="identity")+
          theme_classic()+
          geom_hline(yintercept = 1, col="black", linetype="dashed")+
          ylab("Enrichment bp diff peaks >0 vs all peaks")+
          scale_fill_manual(values=c("#0C255C", "#90D1F0", "#D981AC")))
  print(ggplot(recap_gene, aes(x=feature, y=enrichment_neg, fill=feature))+
          geom_bar(stat="identity")+
          theme_classic()+
          geom_hline(yintercept = 1, col="black", linetype="dashed")+
          ylab("Enrichment bp diff peaks <0 vs all peaks")+
          scale_fill_manual(values=c("#0C255C", "#90D1F0", "#D981AC")))
  # target plot (ggplot)
  print(ggplot(recap_gene, aes(x=feature, y=enrichment_pos, fill=feature))+
          geom_bar(stat="identity")+
          theme_classic()+
          coord_polar()+
          geom_hline(yintercept = 1, col="black", linetype="dashed")+
          ylab("Enrichment bp diff peaks >0 vs all peaks")+
          theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
          scale_fill_manual(values=c("#0C255C", "#90D1F0", "#D981AC")))
  print(ggplot(recap_gene, aes(x=feature, y=enrichment_neg, fill=feature))+
          geom_bar(stat="identity")+
          theme_classic()+
          coord_polar()+
          geom_hline(yintercept = 1, col="black", linetype="dashed")+
          ylab("Enrichment bp diff peaks <0 vs all peaks")+
          theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
          scale_fill_manual(values=c("#0C255C", "#90D1F0", "#D981AC")))
  dev.off()
  
}


######################
##### PARAMETERS #####
######################

archrdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/"
datadir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Differential_peaks/1.Tables"
resdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Peaks_study/Gene_features"

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
promoter_body <- read.table(file.path(annot_dir, "promoter_gene_body_ArchR.bed"), sep="\t", header=F); colnames(promoter_body) <- c("chr", "start", "end", "strand", "gene_feature", "gene")


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
  
  recap_gene <- enrich_gene_features_table(diffPeaks_pos_gr, diffPeaks_neg_gr, promoter_body, peak_set)
  save(recap_gene, file=file.path(resdir, comp, paste0("Diff_study_gene_features_", threshold, ".RData")))
  enrich_gene_features_plots(recap_gene, resdir, comp, threshold)
  
}
  
