
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ArchR)
library(pheatmap)
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


### directory
resdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged"

### load ATAC data 
project_all <- readRDS(file.path(resdir, "Save-ArchR-Project.rds"))
peak_set <- getPeakSet(project_all) # extract peak set

### get gene annotation from ArchR
gene_annot <- getGenes(project_all)
# separate in + and - strand to define TSS
# TSS are the start coordinate for genes on positive strand and end coordinate for negative strand
gene_annot_TSS <- gene_annot
end(gene_annot_TSS[which(strand(gene_annot_TSS)=="+")]) <- start(gene_annot_TSS[which(strand(gene_annot_TSS)=="+")]) 
start(gene_annot_TSS[which(strand(gene_annot_TSS)=="-")]) <- end(gene_annot_TSS[which(strand(gene_annot_TSS)=="-")]) 

### Nearest gene for each peak
nearest_gene <- as.data.frame(GenomicRanges::nearest(peak_set, gene_annot, select="all", ignore.strand=FALSE))
colnames(nearest_gene) <- c("peak", "gene") # replace column names
nearest_gene$peak <- paste0(seqnames(peak_set), "_", start(peak_set), "_", end(peak_set))[nearest_gene$peak] # replace by peak name
nearest_gene$gene <- gene_annot$symbol[nearest_gene$gene] # replace by gene name
nearest_gene2 <- nearest_gene %>% group_by(peak) %>% summarise(gene=paste0(gene, collapse=", ")) %>% as.data.frame() # collapse when several genes associated to peak
save(nearest_gene2, file=file.path(resdir, "nearest_gene_per_peak.Rdata"))

### Distance to nearest TSS for each peak
dist_TSS <- as.data.frame(GenomicRanges::distanceToNearest(peak_set, gene_annot_TSS, ignore.strand=TRUE)) # ignore strand because we work with 1 bp for gene TSS
colnames(dist_TSS) <- c("peak", "gene", "dist") # replace column names
dist_TSS$peak <- paste0(seqnames(peak_set), "_", start(peak_set), "_", end(peak_set))[dist_TSS$peak] # replace by peak name
dist_TSS$gene <- gene_annot$symbol[dist_TSS$gene] # replace by gene name
save(dist_TSS, file=file.path(resdir, "dist_nearest_gene_per_peak.Rdata"))