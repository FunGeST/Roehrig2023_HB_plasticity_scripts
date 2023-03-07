
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(Signac)
library(Seurat)
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

rna_dir <- "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict"
macs2_path <- "/home/amelie/miniconda2/envs/r_env_4.1/bin/macs2"
datadir_all <- "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Signac/1.Preprocess"
frag_dir_all <- "/mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/"
resdir_all <- "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Signac/2.Peak_calling"

name_samples <- c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3377T", "CHC3610T", "CHC3662T")


###########################
##### 0. DATA LOADING #####
###########################

sample_rna <- geco.load(file.path(rna_dir, "Seurat_object_analyses.RData"))

###  retrieve gene annotation for TSS enrichment

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # v86 = hg38 
# seqlevelsStyle(annotation) <- "UCSC" # this function does not work??
seqlevels(annotation) <- paste0("chr", seqlevels(annotation)) # to replace the line above


####################
##### FOR LOOP #####
####################

for(name_sample in name_samples){
  frag_dir <- file.path(frag_dir_all, name_sample, "/outs/")
  datadir <- file.path(datadir_all, name_sample)
  resdir <- file.path(resdir_all, name_sample); if(!dir.exists(resdir)){dir.create(resdir)}
  
  ###########################
  ##### 0. DATA LOADING #####
  ###########################
  
  sample_atac <- geco.load(file.path(datadir, "Signac_object_QC.RData"))
  
  # restrict snRNAseq metadata to tumor cells of the corresponding sample
  metadata_rna <- sample_rna@meta.data[which(sample_rna@meta.data$orig.ident==name_sample),]; rownames(metadata_rna) <- paste0(gsub(paste0(name_sample, "_"), "", rownames(metadata_rna)), "-1")
  metadata_rna$H_LP_M_clean_groups_1 <- as.character(metadata_rna$H_LP_M_clean_groups_1) # otherwise it is a factor
  metadata_rna[which(metadata_rna$SCT_snn_res.0.3==15),"H_LP_M_clean_groups_1"] <- "M_GLI3_cluster"
  
  ##########################################
  ##### 1. RESTRICT TO RNA TUMOR CELLS #####
  ##########################################
  
  sample_atac$subtype <- NA # create tumor subtype column for ATAC data
  sample_atac@meta.data[intersect(colnames(sample_atac), rownames(metadata_rna)),"subtype"] <- metadata_rna[intersect(colnames(sample_atac), rownames(metadata_rna)),"H_LP_M_clean_groups_1"] # keep the same tumor cells/same subtype as snRNAseq
  sample_atac <- sample_atac[,which(!is.na(sample_atac$subtype))] # restrict to these cells
  
  ################################################
  ##### 2. PEAK CALLING ON INDIVIDUAL SAMPLE #####
  ################################################
  
  ##### 1. Call Peaks using MACS2 #####
  #------------------------------------
  
  peaks <- CallPeaks(object = sample_atac, group.by = "subtype", macs2.path = macs2_path)
  # duree: ~20 min
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse") # like ArchR 
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  
  # remove sex chromosomes (X and Y, while ArchR keeps X) cf https://www.nature.com/articles/s41467-021-21583-9
  chroms_to_keep <- which(seqnames(peaks)%in%paste0("chr", 1:22))
  peaks <- peaks[chroms_to_keep,]
  
  
  ##### 2. Create peak/cell counts matrix #####
  #--------------------------------------------
  
  # quantify counts in each peak. Duree: ~15 min
  macs2_counts <- FeatureMatrix(fragments = Fragments(sample_atac), features = peaks, cells = colnames(sample_atac)) # matrix of peaks vs cells, I think it does the sum of insertions using fragment file

  ##### 3. Create Signac object #####
  #----------------------------------
  
  # Create chromatin assay to assign fragments 
  sample_atac2 <- CreateChromatinAssay(macs2_counts, fragments = file.path(frag_dir, "atac_fragments.tsv.gz"), annotation=annotation) # in annotation file we have much more ranges than the number of peaks (normal)
  
  # Create signac object
  sample_atac_macs2 <- CreateSeuratObject(sample_atac2, assay = "ATAC", meta.data=sample_atac@meta.data)
  Idents(sample_atac_macs2) <- sample_atac_macs2$subtype # change object identities to have 1 track per subtype
  
  
  ##### 4. Visualize coverage on example genes #####
  #-------------------------------------------------
  
  # plot Tn5 insertion events, the tracks are normalized using a per-group scaling factor = nb of cells in group x mean sequencing depth for that group of cells
  # normalization to account for differences in nb of cells and potential differences in sequencing depth between groups
  
  pdf(file.path(resdir, "coverage_per_subtype_examples.pdf"))
  CoveragePlot(object = sample_atac_macs2, region = "HNF4A", ranges = peaks, ranges.title = "MACS2", extend.upstream = 2000, extend.downstream = 2000)
  CoveragePlot(object = sample_atac_macs2, region = "LIN28B", ranges = peaks, ranges.title = "MACS2", extend.upstream = 2000, extend.downstream = 2000)
  CoveragePlot(object = sample_atac_macs2, region = "TBX5", ranges = peaks, ranges.title = "MACS2", extend.upstream = 2000, extend.downstream = 2000)
  dev.off()
  
  ################
  ##### SAVE #####
  ################
  
  save(peaks, file=file.path(resdir, "peaks_granges_per_subtype.RData"))
  save(macs2_counts, file=file.path(resdir, "peak_cell_counts_matrix.RData"))
  
  save(sample_atac_macs2, file=file.path(resdir, "Signac_object_MACs2_peaks.RData"))
}

