
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


#####################
##### FUNCTIONS #####
#####################

computeEnrichment_custom <- function(positions = NULL, compare = NULL, background = NULL, r2=r2,r3=r3){
  # We cannot use the PeakAnnoEnrichment function from AchR because it does enrichment on nb of peaks and assumes all peaks are the same size which is not our case here --> might bias results
  
  ##### 1. create recap table to contain results #####
  #---------------------------------------------------
  recap_annot <- data.frame(annot = names(positions), prop_sel=NA, prop_all=NA, stringsAsFactors = F)
  rownames(recap_annot) <- recap_annot$annot
  
  ##### 2. Proportion in SE peaks #####
  #------------------------------------
  # for each annotation compute nb of bp in peaks that intersect the annotation coordinates
  #i=1
  for(a in recap_annot$annot){
    #print(i)
    positions_a <- positions[[a]]
    strand(positions_a) <- "*" # the strand does not matter
    # get intersections in all peaks, selected peaks
    gr_intersect_all <- GenomicRanges::intersect(r2, positions_a)
    gr_intersect_sel <- GenomicRanges::intersect(r3[compare], positions_a) # get the intersected nucleotides between peaks and state coordinates 
    recap_annot[a,"prop_all"] <- sum(width(gr_intersect_all))
    recap_annot[a,"prop_sel"] <- sum(width(gr_intersect_sel)) # how many bp in total in peaks intersect the annotation
    
    #i=i+1 # increment
  }
  #compute proportion
  recap_annot$prop_all <- 100*recap_annot$prop_all/sum(width(r2)) # divide by total nb of bp (here warning: annotations can overlap each other --> use the global sets)
  recap_annot$prop_sel <- 100*recap_annot$prop_sel/sum(width(r3[compare])) 
  
  ##### 3. Compute enrichment of diff peaks vs all peaks #####
  #-----------------------------------------------------------
  recap_annot$enrichment <- recap_annot$prop_sel / recap_annot$prop_all
  return(recap_annot) 
}

peakAnnoEnrichment_custom <- function (seMarker = NULL, ArchRProj = NULL, positions = NULL,
          cutOff = "FDR <= 0.1 & Log2FC >= 0.5", background = "all",
          logFile = createLogFile("peakAnnoEnrichment"))
{
  if (metadata(seMarker)$Params$useMatrix != "PeakMatrix") {
    stop("Only markers identified from PeakMatrix can be used!")
  }
  
  # object for whole peak set
  r2 <- getPeakSet(ArchRProj)
  pr2 <- paste(seqnames(r2), start(r2), end(r2), sep = "_")
  mcols(r2) <- NULL
  # object for seMarker = peaks to study
  r3 <- GRanges(rowData(seMarker)$seqnames, IRanges(rowData(seMarker)$start,
                                                    rowData(seMarker)$end))
  pr3 <- paste(seqnames(r3), start(r3), end(r3), sep = "_")
  mcols(r3) <- NULL
  
  if (length(which(pr2 %ni% pr3)) != 0) {
    stop("Peaks from seMarker do not match peakSet in ArchRProj!")
  }
  assayNames <- names(SummarizedExperiment::assays(seMarker)) # colnames of seMarker
  for (an in assayNames) {
    eval(parse(text = paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['",
                             an, "']]")))
  }
  passMat <- eval(parse(text = cutOff)) # restrict to peaks with the wanted thresholds
  passMat[is.na(passMat)] <- FALSE
  for (an in assayNames) {
    eval(parse(text = paste0("rm(", an, ")")))
  }
  background = seq_len(length(r2))
  enrichList <- lapply(seq_len(ncol(seMarker)), function(x) {
    idx <- which(passMat[, x])
    computeEnrichment_custom(positions = positions, compare = idx, background = background, r2=r2, r3=r3)
  }) %>% SimpleList # takes ~20 min
  print("test2")
  names(enrichList) <- colnames(seMarker)
  assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x) {
    d <- lapply(seq_along(enrichList), function(y) {
      enrichList[[y]][names(positions), x, drop = FALSE]
    }) %>% Reduce("cbind", .)
    colnames(d) <- names(enrichList)
    d
  }) %>% SimpleList
  names(assays) <- names(enrichList[[1]])
  assays <- rev(assays)
  return(assays)
}

plot_motif_enrich <- function(df, labels){ggplot(df, aes(rank, enrichment, color = enrichment)) + 
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      max.overlaps=Inf,
      data = df[match(labels, df$TF), ], aes(x = rank, y = enrichment, label = TF), 
      size = 6,
      nudge_x = 2,
      color = "black"
    ) + theme_ArchR() + 
    ylab("Enrichment") + 
    xlab("Rank Sorted features Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))+
    theme(axis.text=element_text(size=13), axis.title=element_text(size=15))}


######################
##### PARAMETERS #####
######################

archrdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/"
datadir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Differential_peaks/1.Tables"
resdir = file.path("/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Motif_TF")
cisbp_dir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/0.Public_annotation/Cisbp_database/" # folder containing cisbp v2.00 results from MEME software

# choose one comparison to study 
comp = c("NT_vs_T", "H_vs_LP", "H_HLP_LP_vs_M")[2] 

# define thresholds for differential peaks study 
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
metadata <- getCellColData(project_all) # get cell metadata

# Extract peak matrix
peak_matrix_se <- getMatrixFromProject(ArchRProj = project_all, useMatrix = "PeakMatrix", binarize=F) # is a summarized experiment. Raw insertion counts. There was a ceiling of max 4 counts per peak when creating the peak matrix (cf ArchR doc https://github.com/GreenleafLab/ArchR/discussions/943)
peak_matrix <- peak_matrix_se@assays@data[["PeakMatrix"]] # extract the matrix slot
r1 <- rowRanges(peak_matrix_se) # rowRanges will give the exact order in the matrix (do that for each extracted matrix from ArchR project because the peak order can change)
peak_names <- paste0(seqnames(r1), ":", start(r1), "-", end(r1))
rownames(peak_matrix) <- peak_names

# Load differential tables
if(comp=="NT_vs_T"){diffPeaks_SE <- geco.load(file.path(datadir, "diffPeaks_NT_vs_T_summarized_exp.RData")) # T is >0, NT is <0
}else if(comp=="H_vs_LP"){diffPeaks_SE <- geco.load(file.path(datadir, "diffPeaks_H_vs_LP_summarized_exp.RData")) # LP is >0, H is <0
}else if(comp=="H_HLP_LP_vs_M"){diffPeaks_SE <- geco.load(file.path(datadir, "diffPeaks_H_HLP_LP_vs_M_summarized_exp.RData"))} # M is >0, H/H+LP/LP is <0

# Load cisbp v2.00 motifs
PWM_cisbp_v2 <- geco.load(file.path(cisbp_dir, "PWM_list_cisbp_v2.00.RData")) 

# Load annotation files 
positions_v1 <- readRDS("/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Annotations/Motif_cisbp-Positions-In-Peaks.rds")
names(positions_v1) <- paste0(gsub("\\_.*", "", names(positions_v1)), "_v1") # change nomenclature of TF names: remove the number after the underscore, add the version of cisbp
positions_v2 <- readRDS("/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Annotations/Motif_cisbp_v2.00-Positions-In-Peaks.rds")
names(positions_v2) <- paste0(gsub("\\_.*", "", names(positions_v2)), "_v2") # change nomenclature of TF names: remove the number after the underscore, add the version of cisbp
positions <- c(positions_v1, positions_v2) # combine cisbp v1 and v2 since some TF have better results sometimes in v1 sometimes in v2
# comprises 1244 TF in total (691 common v1/v2, 179 exclusive v1, 374 exclusive v2)
  

#####################################################
##### 1. MOTIF ENRICHMENT IN DIFFERENTIAL PEAKS #####
#####################################################

##### 1. Annotate peaks #####
#----------------------------
# use the v2.00 version from CisBP (custom annotation)
#project_all <- addMotifAnnotations(ArchRProj = project_all, motifPWMs=PWM_cisbp_v2, name = "Motif_cisbp_v2.00") # takes ~2h

##### 2. Perform enrichment in differential peaks #####
#------------------------------------------------------
for(threshold in c("loose_thresholds", "medium_thresholds", "strict_thresholds")){
  if(threshold == "loose_thresholds"){FDR_threshold <- threshold_loose_FDR; Log2FC_threshold <- threshold_loose_logFC}
  if(threshold == "medium_thresholds"){FDR_threshold <- threshold_medium_FDR; Log2FC_threshold <- threshold_medium_logFC}
  if(threshold == "strict_thresholds"){FDR_threshold <- threshold_strict_FDR; Log2FC_threshold <- threshold_strict_logFC}
  
  # Get motifs in differential peaks >0 and <0 for the comparison, takes some time
  motifs_pos <- peakAnnoEnrichment_custom(
    seMarker = diffPeaks_SE,
    ArchRProj = project_all,
    positions = positions,
    cutOff = paste0("FDR <= ", FDR_threshold, " & Log2FC >= ", Log2FC_threshold)
  )
  
  motifs_neg <- peakAnnoEnrichment_custom(
    seMarker = diffPeaks_SE,
    ArchRProj = project_all,
    positions = positions,
    cutOff = paste0("FDR <= ", FDR_threshold, " & Log2FC <= ", -Log2FC_threshold)
  )
  
  # Rank TF motifs by enrichment
  df_pos <- data.frame(TF = rownames(motifs_pos[["enrichment"]]), enrichment = motifs_pos[["enrichment"]]); colnames(df_pos) <- c("TF", "enrichment")
  df_pos <- df_pos[order(df_pos$enrich, decreasing = TRUE),]
  df_pos$rank <- seq_len(nrow(df_pos))
  df_neg <- data.frame(TF = rownames(motifs_neg[["enrichment"]]), enrichment = motifs_neg[["enrichment"]]); colnames(df_neg) <- c("TF", "enrichment")
  df_neg <- df_neg[order(df_neg$enrich, decreasing = TRUE),]
  df_neg$rank <- seq_len(nrow(df_neg))
  save(df_pos, file=file.path(resdir, comp, paste0("ranked_TF_motif_enrichment_pos_cisbp_v1_v2_", threshold, ".RData")))
  save(df_neg, file=file.path(resdir, comp, paste0("ranked_TF_motif_enrichment_neg_cisbp_v1_v2_", threshold, ".RData")))
  
  # Remove info of v1/v2 for the plot and keep the best enrichment (sometimes it will be for v1, sometimes for v2)
  df_pos$TF <- gsub("\\_.*", "", df_pos$TF); df_pos <- df_pos[!duplicated(df_pos$TF), ];  df_pos$rank <- seq_len(nrow(df_pos))
  df_neg$TF <- gsub("\\_.*", "", df_neg$TF); df_neg <- df_neg[!duplicated(df_neg$TF), ];  df_neg$rank <- seq_len(nrow(df_neg))
  
  # Select labels to display on plot
  if(comp=="H_HLP_LP_vs_M"){
    labels_pos <- c("PURG", "NR5A1", "TWIST1", "ZFPM2", "LEF1", "RUNX2", "DZIP1", "TSHZ3", "DACH2", "DACH1", "CSRNP3", "ZMAT4", "AEBP1", "FOSL2", "RAG1", "ZIC2", "TRERF1",
                    "NFATC4", "PRDM2", "ZIC5", "NFE2L3", "AHRR", "ESR1", "RORA", "ZFP28", "PRRX2", "TRPS1", "ALX3", "ZIC1", "ZIC4")
    labels_pos <- intersect(labels_pos, df_pos[which(df_pos$rank <= 300),"TF"]) # keep TF in the top 300
    labels_neg <- c("L3MBTL4", "HNF4A", "SALL1", "HNF1A", "ONECUT2", "FOXA1", "NR1I2", "ONECUT1", "THRB", "PRDM16", "FOXA3", "HNF4G", "CREB3L3", "FOXP4", "FOXP2", "GATA4",
                    "PEG3", "NR1H4", "ATF5", "NR1I3", "CEBPA", "FOXP1", "USF3", "RXRA", "ZBED3", "CEBPG", "NR2F6", "ZXDC", "NCOA3", "ZKSCAN1", "PPARG", "PURA", "REPIN1",
                    "CENPX", "BCL6", "SLC2A4RG", "ZFPM1", "HHEX", "SALL4", "GLMP", "RREB1", "IKZF5", "KLF15", "TSHZ2","FOXQ1", "HNF1B")
    labels_neg <- intersect(labels_neg, df_neg[which(df_neg$rank <= 300),"TF"]) # keep TF in the top 300
  }else if(comp=="H_vs_LP"){
    labels_pos <- c("ST18", "CBX2", "SOX4", "PCGF2", "HIF3A", "AC092835.1", "MYT1", "EBF4", "ZFP37", "SALL4", "RFX3", "ZFHX2", "TEAD4", "PITX1", "TEAD2", "DACH1", "LEF1",
                    "HOXA3", "NFE2L3", "GTF2IRD1", "ETV4", "MAZ", "ETV5", "OTX1", "DACH2", "VAX2", "TIGD4", "TET3", "NFXL1", "MNX1", "SOX9", "HMGA1", "ZFP30", "DLX1", 
                    "FOXM1", "HOXB7")
    labels_pos <- intersect(labels_pos, df_pos[which(df_pos$rank <= 300),"TF"]) # keep TF in the top 300
    labels_neg <- c("AR", "ZMAT1", "CEBPD", "EPAS1", "JAZF1", "CEBPB", "HLF", "NPAS3", "IKZF1", "SP140", "NR1H4", "NFIL3", "BNC2", "SATB1", "SP100", "BHLHE40", "NR3C1",
                    "NR3C2", "RORA", "ARNTL", "KLF4", "ZBTB16", "STAT3", "AHR", "SP140L", "SOX7", "ATF5", "RORC", "HIVEP2", "MITF", "ATF3", "HIVEP1", "TFEC", "ZBTB38", 
                    "STAT4", "IKZF5", "BACH1", "NFIA", "IRF5", "HIVEP3", "PPARA", "CREB5", "MLXIPL", "ZFPM2", "PPARG", "KDM2A", "ZBTB21", "CEBPA")
    labels_neg <- intersect(labels_neg, df_neg[which(df_neg$rank <= 300),"TF"]) # keep TF in the top 300
  }

  # Show TF motifs ranked, color by enrichment
  pdf(file.path(resdir, comp, paste0("ranked_TF_motif_enrichment_cisbp_v1_v2_", threshold, ".pdf")))
  plot_motif_enrich(df_pos, labels_pos) + ggtitle("Differential peaks >0")
  plot_motif_enrich(df_neg, labels_neg) + ggtitle("Differential peaks <0")
  dev.off()
 

}

###################
##### . SAVE #####
###################

#saveArchRProject(ArchRProj = project_all)


