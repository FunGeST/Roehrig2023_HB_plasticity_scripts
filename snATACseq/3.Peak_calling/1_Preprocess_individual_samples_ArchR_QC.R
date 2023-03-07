
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
library(MASS)
library(viridis)
set.seed(1234)

geco.load <- function (filename) # File to be loaded
{
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error - function geco.load : file ", filename, 
            " not found"))
  NULL
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

######################
##### PARAMETERS #####
######################

name_sample = c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")[1]

datadir <- paste0("/mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/", name_sample, "/outs/")
resdir <- file.path("/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Signac/1.Preprocess", name_sample); if(!dir.exists(resdir)){dir.create(resdir)}

archr_dir <- paste0("/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/", name_sample)


###########################
##### 0. DATA LOADING #####
###########################

### Load unnormalized peak counts matrix: used for QC

counts <- Read10X(data.dir=file.path(datadir, "filtered_feature_bc_matrix")) # will load both RNA and ATAC data from cellranger ARC, returns a list of 2 objects
counts <- counts[[2]] # retrieve ATAC data 

### Load metadata 

metadata_peaks <- read.table(file=file.path(datadir, "atac_peak_annotation.tsv"), header=T, sep="\t")
metadata <- read.csv(file=file.path(datadir, "per_barcode_metrics.csv"), header=T, row.names=1)

### Retrieve gene annotation for TSS enrichment

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # we get a warning "The 2 combined objects have no sequence levels in common.", seems normal because occurs in do.call() in  GetGRangesFromEnsDb when the function transforms a list of regions for each chromosome into a single object --> there are no seq in common between chr1 and chr2 for example.... 
# seqlevelsStyle(annotation) <- "UCSC" # this function does not work??
seqlevels(annotation) <- paste0("chr", seqlevels(annotation)) # to replace the line above

### Load cells to keep after ArchR QC 

archr_cells <- geco.load(file.path(archr_dir, "metadata.RData"))
archr_cells <- gsub(paste0(name_sample, "#"), "", rownames(archr_cells)) # keep same nomenclature as cellranger


########################################
##### 1. CREATION OF SIGNAC OBJECT #####
########################################

### Create chromatin assay for Signac

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = NULL, # default is hg19 but according to developers this parameter is not compulsory for CreateChromatinAssay (https://github.com/timoast/signac/issues/232)
  fragments = file.path(datadir, "atac_fragments.tsv.gz"),
  min.cells = 0, # default is 10
  min.features = 0, # default is 200
  annotation=annotation
)

### Create Seurat object from chromatin assay

sample_atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = NULL
)



###################################
##### 2. QC METRICS BEFORE QC #####
###################################

##### 1. Compute metrics #####
#-----------------------------

# count total fragments per cell barcode

fragment_counts <- CountFragments(file.path(datadir, "atac_fragments.tsv.gz"))
rownames(fragment_counts) <- fragment_counts$CB
fragment_counts <- fragment_counts[colnames(sample_atac),] # restrict to the cells identified as non empty by cellranger
# fragment_counts containsn 5 columns: "CB" = cell barcode, "frequency_count" = total nb of fragments sequenced for the cell, "mononucleosome" and "nucleosome free"= total nb of fragments with 147bp<length<294bp or length<147bp,
# "reads_count"= total nb of reads sequenced for the cell
sample_atac$nb_fragments <- fragment_counts$frequency_count # add to signac object

# compute nucleosome signal score per cell
sample_atac <- NucleosomeSignal(object = sample_atac)

# compute TSS enrichment score per cell
sample_atac <- TSSEnrichment(object = sample_atac, fast = FALSE)



##### 2. Visualize metrics #####
#-------------------------------

pdf(file.path(resdir, "Metrics_before_QC.pdf"))

# TSS enrichment vs log10(nb fragments)
df_metrics <- sample_atac@meta.data[,c("nb_fragments", "TSS.enrichment")] # extract metrics of interest
df_metrics$log10_nb_fragments <- log10(df_metrics$nb_fragments)
df_metrics$density <- get_density(df_metrics$log10_nb_fragments, df_metrics$TSS.enrichment, n = 5000) # n: create a square n by n grid to compute density

ggplot(df_metrics) + geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + scale_color_viridis() 

# Nucleosome signal score
sample_atac$nucleosome_group <- ifelse(sample_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_atac, group.by = 'nucleosome_group')

# Violin plots 
VlnPlot(
  object = sample_atac,
  features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
  pt.size = 0.1,
  ncol = 5
)
dev.off()


##########################################################
##### 3. QC METRICS AFTER RESTRICTING TO ARCHR CELLS ##### 
##########################################################

##### 1. Restrict to ArchR cells #####
#-----------------------------------

# Intersection of cellranger cells in peak matrix and ArchR cells (QC has already been done in ArchR)

intersect_cellranger_ArchR <- intersect(archr_cells, colnames(counts)) # usually makes > 90% of ArchR cells, ~70% of cellranger cells
counts <- counts[,intersect_cellranger_ArchR]

sample_atac_qc <- sample_atac[,intersect_cellranger_ArchR]


##### 2. Visualize metrics #####
#-------------------------------

pdf(file.path(resdir, "Metrics_after_QC.pdf"))

# TSS enrichment vs log10(nb fragments)
df_metrics <- sample_atac_qc@meta.data[,c("nb_fragments", "TSS.enrichment")] # extract metrics of interest
df_metrics$log10_nb_fragments <- log10(df_metrics$nb_fragments)
df_metrics$density <- get_density(df_metrics$log10_nb_fragments, df_metrics$TSS.enrichment, n = 5000) # n: create a square n by n grid to compute density

ggplot(df_metrics) + geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + scale_color_viridis() 

# Violin plots 
VlnPlot(
  object = sample_atac_qc,
  features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
  pt.size = 0.1,
  ncol = 5
)
dev.off()



################
##### SAVE #####
################

save(sample_atac, file=file.path(resdir, "Signac_object_noQC.RData"))
save(sample_atac_qc, file=file.path(resdir, "Signac_object_QC.RData"))

