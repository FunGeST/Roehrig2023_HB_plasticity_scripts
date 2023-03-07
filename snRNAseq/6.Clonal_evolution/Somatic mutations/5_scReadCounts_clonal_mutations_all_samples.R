
library(dplyr)
library(ggplot2)
library(infercnv)
library(phytools)
library(dendextend)
library(Seurat)


geco.load <- function (filename) # File to be loaded
{
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error - function geco.load : file ", filename, 
            " not found"))
  NULL
}

factoall <- function (d, ncmax = 10) 
{
  n <- ncol(d)
  for (i in 1:n) {
    if (is.factor(d[, i])) {
      d[, i] <- as.character(d[, i])
      na <- which(is.na(d[, i]))
      num <- suppressWarnings(as.numeric(d[, i]))
      nanum <- which(is.na(num))
      if (length(nanum) == length(na)) {
        #                int <- suppressWarnings(as.integer(d[, i]))
        #                naint <- which(is.na(int))
        #                nc <- nchar(num)
        #                if (length(naint) == length(nanum) & all(nc < ncmax)) {
        #                  d[, i] <- int
        #                }
        #                else {
        d[, i] <- num
        #                }
      }
    }
  }
  d
}

scale_fill_custom <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(c("red", "orange", "yellow", "grey20", "grey"), c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known")), 
    ...
  )
}


######################
##### PARAMETERS #####
######################

all_samples <- c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")
tumor_samples <- c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3377T", "CHC3610T", "CHC3662T")

datadir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/"
screadcounts_dir <- file.path("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/Somatic_mutations")
seurat_dir_tumor <- file.path("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict")
seurat_dir_all <- file.path("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/")
CNV_dir <- paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/", all_samples[1], "/Final_clusters")


###########################
##### 0. DATA LOADING #####
###########################

##### 1. WGS info #####
#----------------------

vcf_wgs <- geco.load(file.path(datadir, "Merged_pediatric_129s_89p_vcf_all_with_CCF_and_SBS_proba.RData")) #last version, contains all tumors but is in hg19
vcf_wgs <- vcf_wgs %>% filter(Sample %in% all_samples) %>% filter(Clonality == "clonal") %>% filter(FILTER == "PASS") # keep only clonal mutations for the samples of interest
vcf_wgs$mutation_hg38 <- paste(vcf_wgs$CHROM_hg38, vcf_wgs$POS_hg38, sep="_")
vcf_wgs <- vcf_wgs[which(vcf_wgs$CHROM == vcf_wgs$CHROM_hg38),] # remove potential liftover artifacts


##### 2. snMultiome info #####
#-----------------------------
### Load Seurat objects ###
sample_seurat_tumor <- geco.load(file.path(seurat_dir_tumor, "Seurat_object_analyses.RData")) # tumor cells
sample_seurat_all <- geco.load(file.path(seurat_dir_all, "Seurat_object_analyses_all_cells.RData")) # all cells

# add tumor info in sample_seurat_all
sample_seurat_all$tumor_cells <- "non tumor"
sample_seurat_all@meta.data[colnames(sample_seurat_tumor),"tumor_cells"] <- "tumor"

### Create merged table containing screadcount results for each sample ###
for(s1 in all_samples){
  if(s1=="CHC2959N"){ # 1st sample to study
    screadcounts_results <- read.table(file.path(screadcounts_dir, s1, "Screadcounts_output", paste0("scReadCounts_output_xf_25_", gsub("CHC", "", "CHC2959T"), "_mut.tsv")), sep="\t", header=T)
    screadcounts_results$ReadGroup <- gsub("-1", "", paste0(s1, "_", screadcounts_results$ReadGroup))
    screadcounts_results$sample <- rep(s1, times=nrow(screadcounts_results))
    screadcounts_results$CHROM <- paste0("chr", screadcounts_results$CHROM) # define mutation name
    screadcounts_results$mutation_hg38 <- paste0(screadcounts_results$CHROM, "_", screadcounts_results$POS) # define mutation name
    screadcounts_results <- screadcounts_results %>% filter(mutation_hg38 %in% vcf_wgs$mutation_hg38) # restrict to clonal mutations with filter=PASS
    screadcounts_results$sample_mutations <- rep("CHC2959T", times=nrow(screadcounts_results)) # annotate from which sample mutations come from
    for(s2 in tumor_samples[2:length(tumor_samples)]){
      screadcounts_results. <- read.table(file.path(screadcounts_dir, s1, "Screadcounts_output", paste0("scReadCounts_output_xf_25_", gsub("CHC", "", s2), "_mut.tsv")), sep="\t", header=T)
      screadcounts_results.$ReadGroup <- gsub("-1", "", paste0(s1, "_", screadcounts_results.$ReadGroup))
      screadcounts_results.$sample <- rep(s1, times=nrow(screadcounts_results.))
      screadcounts_results.$CHROM <- paste0("chr", screadcounts_results.$CHROM) # define mutation name
      screadcounts_results.$mutation_hg38 <- paste0(screadcounts_results.$CHROM, "_", screadcounts_results.$POS) # define mutation name
      screadcounts_results. <- screadcounts_results. %>% filter(mutation_hg38 %in% vcf_wgs$mutation_hg38) # restrict to clonal mutations with filter=PASS
      screadcounts_results.$sample_mutations <- rep(s2, times=nrow(screadcounts_results.)) # annotate from which sample mutations come from
      screadcounts_results <- rbind(screadcounts_results, screadcounts_results.)
    }
  }else{
    for(s2 in tumor_samples){
      screadcounts_results. <- read.table(file.path(screadcounts_dir, s1, "Screadcounts_output", paste0("scReadCounts_output_xf_25_", gsub("CHC", "", s2), "_mut.tsv")), sep="\t", header=T)
      screadcounts_results.$ReadGroup <- gsub("-1", "", paste0(s1, "_", screadcounts_results.$ReadGroup))
      screadcounts_results.$sample <- rep(s1, times=nrow(screadcounts_results.))
      screadcounts_results.$CHROM <- paste0("chr", screadcounts_results.$CHROM) # define mutation name
      screadcounts_results.$mutation_hg38 <- paste0(screadcounts_results.$CHROM, "_", screadcounts_results.$POS) # define mutation name
      screadcounts_results. <- screadcounts_results. %>% filter(mutation_hg38 %in% vcf_wgs$mutation_hg38) # restrict to clonal mutations with filter=PASS
      screadcounts_results.$sample_mutations <- rep(s2, times=nrow(screadcounts_results.)) # annotate from which sample mutations come from
      screadcounts_results <- rbind(screadcounts_results, screadcounts_results.)
    }
  }
}

screadcounts_results <- screadcounts_results[which(screadcounts_results$ReadGroup %in% colnames(sample_seurat_all)),] # screadcounts takes all cells from BAM --> we need to restrict to post-QC cells in Seurat object
save(screadcounts_results, file="C:/Users/amelie/Dropbox/Publication_multiome/Dossier_pour_Eric/Mutations/Results_screadcounts_table.RData")

# Annotate bulk mutations detected in snRNAseq
screadcounts_results_mut <- screadcounts_results[which(screadcounts_results$SNVCount>0),] # restrict to mutated positions
vcf_wgs$detected_in_snRNAseq <- ifelse(vcf_wgs$mutation_hg38 %in% screadcounts_results_mut$mutation_hg38, "yes", "no")

vcf_wgs. <- vcf_wgs[,c("Sequencing", "Sample", "Patient.identification", "mutation_hg38", "detected_in_snRNAseq", "CHROM", "POS", "CHROM_hg38", "POS_hg38", "REF", "ALT", "Hugo_Symbol", "Protein_Change", "FILTER")]
save(vcf_wgs., file="C:/Users/amelie/Dropbox/Publication_multiome/Dossier_pour_Eric/Mutations/Mutations_bulk_WGS.RData")



#######################################################
##### 1. TABLE OF SAMPLE MUTATIONS IN EACH SAMPLE #####
#######################################################

# Create table for all cells or tumor cells with UMAP coordinates
mutations_all_cells <- data.frame(cells = colnames(sample_seurat_all), UMAP_1 =sample_seurat_all@reductions$umap@cell.embeddings[colnames(sample_seurat_all), c("UMAP_1")], UMAP_2 =sample_seurat_all@reductions$umap@cell.embeddings[colnames(sample_seurat_all), c("UMAP_2")], stringsAsFactors = F)
mutations_tumor_cells <- data.frame(cells = colnames(sample_seurat_tumor), UMAP_1 =sample_seurat_tumor@reductions$umap@cell.embeddings[colnames(sample_seurat_tumor), c("UMAP_1")], UMAP_2 =sample_seurat_tumor@reductions$umap@cell.embeddings[colnames(sample_seurat_tumor), c("UMAP_2")], stringsAsFactors = F)
rownames(mutations_all_cells) <- mutations_all_cells$cells
rownames(mutations_tumor_cells) <- mutations_tumor_cells$cells

# Add mutation detection per cell
mutations_all_cells$mut_CHC2959T <- mutations_all_cells$mut_CHC2960T <- mutations_all_cells$mut_CHC3133T <- mutations_all_cells$mut_CHC3377T <- mutations_all_cells$mut_CHC3610T <- mutations_all_cells$mut_CHC3662T <- NA
mutations_tumor_cells$mut_CHC2959T <- mutations_tumor_cells$mut_CHC2960T <- mutations_tumor_cells$mut_CHC3133T <- mutations_tumor_cells$mut_CHC3377T <- mutations_tumor_cells$mut_CHC3610T <- mutations_tumor_cells$mut_CHC3662T <- NA

for(s in tumor_samples){
  # sum REF/ALT counts per cell for the clonal mutations of the sample under study
  recap_by_cell <- screadcounts_results %>% 
    filter(sample_mutations == s) %>%
    group_by(ReadGroup) %>% 
    summarise(SNV_total_counts = sum(SNVCount), Ref_total_counts = sum(RefCount), VAF = SNV_total_counts / (SNV_total_counts + Ref_total_counts)) %>%
    as.data.frame()
  rownames(recap_by_cell) <- recap_by_cell$ReadGroup
  mutations_all_cells[intersect(rownames(mutations_all_cells), recap_by_cell$ReadGroup), paste0("mut_", s)] <- recap_by_cell[intersect(rownames(mutations_all_cells), recap_by_cell$ReadGroup), "SNV_total_counts"]
  mutations_tumor_cells[intersect(rownames(mutations_tumor_cells), recap_by_cell$ReadGroup), paste0("mut_", s)] <- recap_by_cell[intersect(rownames(mutations_tumor_cells), recap_by_cell$ReadGroup), "SNV_total_counts"]
  
}
save(mutations_all_cells, file="C:/Users/amelie/Dropbox/Publication_multiome/Dossier_pour_Eric/Mutations/Mutations_screadcounts_all_cells.RData")
save(mutations_tumor_cells, file="C:/Users/amelie/Dropbox/Publication_multiome/Dossier_pour_Eric/Mutations/Mutations_screadcounts_tumor_cells.RData")

# Turn to binary tables (consider >=1 mutation = "1", else if no mutation: "0")
mutations_all_cells. <- mutations_all_cells
for(c in grep("mut_CHC", colnames(mutations_all_cells.), value=T)){
  mutations_all_cells.[,c] <- ifelse(mutations_all_cells[,c] == 0, 0, 1)
}
mutations_all_cells. <- mutations_all_cells.[sample(rownames(mutations_all_cells.)),] # randomize cells
mutations_tumor_cells. <- mutations_tumor_cells
for(c in grep("mut_CHC", colnames(mutations_tumor_cells.), value=T)){
  mutations_tumor_cells.[,c] <- ifelse(mutations_tumor_cells[,c] == 0, 0, 1)
}
mutations_tumor_cells. <- mutations_tumor_cells.[sample(rownames(mutations_tumor_cells.)),] # randomize cells

# Plot on UMAP
ggplot(mutations_all_cells., aes(x=UMAP_1, y=UMAP_2, col=mut_CHC2959T))+
  geom_point(size=1)+
  theme_classic()+
  scale_color_gradientn( colours = c("grey50", "red"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
ggplot(mutations_tumor_cells., aes(x=UMAP_1, y=UMAP_2, col=mut_CHC2959T))+
  geom_point(size=1)+
  theme_classic()+
  scale_color_gradientn( colours = c("grey50", "red"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))

