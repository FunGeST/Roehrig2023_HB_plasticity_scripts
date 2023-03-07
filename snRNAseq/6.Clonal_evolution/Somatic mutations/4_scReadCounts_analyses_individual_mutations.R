

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

name_sample <- c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")
mutation <- c("beta_cat", "APC", "ARID1A", "DDX3X")[1] # IRF2: no mutations found in bulk in those samples

datadir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/"
screadcounts_dir <- file.path("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/Somatic_mutations")
seurat_dir_tumor <- file.path("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict")
seurat_dir_all <- file.path("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/")
resdir <- file.path(screadcounts_dir, "Individual_mutations")


###########################
##### 0. DATA LOADING #####
###########################

### snMultiome info ###

### Create merged table containing screadcount results for the sample ###
# initialize
screadcounts_results <- read.table(file.path(screadcounts_dir, name_sample[1], "Screadcounts_output", paste0("scReadCounts_output_", mutation, "_mut.tsv")), sep="\t", header=T) # WARNING here we do not use the filtered BAM (xf=25) otherwise very few cells with mutation detection
screadcounts_results$ReadGroup <- gsub("-1", "", paste0(name_sample[1], "_", screadcounts_results$ReadGroup))
screadcounts_results$sample <- rep(name_sample[1], times=nrow(screadcounts_results))
# loop on the 2 NT samples
for(n in name_sample[2:length(name_sample)]){
  # Global output
  screadcounts_results. <- read.table(file.path(screadcounts_dir, n, "Screadcounts_output", paste0("scReadCounts_output_", mutation, "_mut.tsv")), sep="\t", header=T)
  if(nrow(screadcounts_results.)!=0){ # for some mutations there might not be any rows in the screadcount results
    screadcounts_results.$ReadGroup <- gsub("-1", "", paste0(n, "_", screadcounts_results.$ReadGroup))
    screadcounts_results.$sample <- rep(n, times=nrow(screadcounts_results.))
    screadcounts_results <- rbind(screadcounts_results, screadcounts_results.)  
  }
}

screadcounts_results$CHROM <- paste0("chr", screadcounts_results$CHROM)
screadcounts_results$mutation <- paste(screadcounts_results$CHROM, screadcounts_results$POS, sep="_") # create a column to easily retrieve the mutation position
length(unique(screadcounts_results$mutation)) 

save(screadcounts_results, file=file.path(resdir, paste0("screadcounts_results_", mutation, "_mutation_info.RData")))

### Load Seurat objects ###
sample_seurat_tumor <- geco.load(file.path(seurat_dir_tumor, "Seurat_object_analyses.RData")) # tumor cells
sample_seurat_all <- geco.load(file.path(seurat_dir_all, "Seurat_object_analyses_all_cells.RData")) # all cells

# add tumor info in sample_seurat_all
sample_seurat_all$tumor_cells <- "non tumor"
sample_seurat_all@meta.data[colnames(sample_seurat_tumor),"tumor_cells"] <- "tumor"


###################################################
##### 1. VISUALIZE MUTATIONS ON COMPLETE UMAP #####
###################################################

##### 1. Concatenate mutation info per cell #####
#------------------------------------------------

# Cells without mutation
recap_by_cell_no_mut <- screadcounts_results %>% 
  group_by(ReadGroup) %>% # we collapse all mutations because then cells with SNV = 0 will mean that the cell have not even one of the mutations (one cell can have 0 SNV counts for 1 mutation but 1 SNV count for another)
  summarise(SNV_total_counts = sum(SNVCount), Ref_total_counts = sum(RefCount), VAF = SNV_total_counts / (SNV_total_counts + Ref_total_counts))
recap_by_cell_no_mut <- as.data.frame(recap_by_cell_no_mut[which(recap_by_cell_no_mut$SNV_total_counts!=0 | recap_by_cell_no_mut$Ref_total_counts!=0 ),]) # remove rows that have 0 SNV and 0 Ref counts = not covered (~1% of the rows)
recap_by_cell_no_mut <- recap_by_cell_no_mut[which(recap_by_cell_no_mut$ReadGroup%in%colnames(sample_seurat_all)),] # screadcounts takes cells without QC -->  some of them are not in our postQC seurat object
recap_by_cell_no_mut <- recap_by_cell_no_mut %>% filter(SNV_total_counts == 0) # restrict to cells without mutation

# Cells with mutation
recap_by_cell_mut <- screadcounts_results %>% 
  group_by(ReadGroup, mutation) %>% # here we separate each mutation because we want to know which mutation is found in which cell 
  summarise(SNV_total_counts = sum(SNVCount), Ref_total_counts = sum(RefCount), VAF = SNV_total_counts / (SNV_total_counts + Ref_total_counts))
recap_by_cell_mut <- as.data.frame(recap_by_cell_mut[which(recap_by_cell_mut$SNV_total_counts!=0 | recap_by_cell_mut$Ref_total_counts!=0 ),]) # remove rows that have 0 SNV and 0 Ref counts = not covered (~1% of the rows)
recap_by_cell_mut <- recap_by_cell_mut[which(recap_by_cell_mut$ReadGroup%in%colnames(sample_seurat_all)),] # screadcounts takes cells without QC -->  some of them are not in our postQC seurat object
recap_by_cell_mut <- recap_by_cell_mut %>% filter(SNV_total_counts > 0) # restrict to cells with at least 1 mutation 

mut_in_cells <- reshape2::dcast(recap_by_cell_mut, ReadGroup ~ mutation, value.var = "SNV_total_counts") # recreate classical table from recap_by_cell_mut to see which mutations in which cells with the nb of SNV counts
for(m in colnames(mut_in_cells)[2:ncol(mut_in_cells)]){ # for each mutation column
  mut_in_cells[which(!is.na(mut_in_cells[,m])),m] <- m # replace each SNV count by the corresponding mutation
}
mut_in_cells$mutation_concat <- apply(mut_in_cells[2:ncol(mut_in_cells)], 1, function(x) paste(x[!is.na(x)], collapse=",")) # for each cell concatenate the mutation info
  

##### 2. Add mutation info to seurat object #####
#------------------------------------------------
sample_seurat_all$mut <- "not known"
sample_seurat_all@meta.data[recap_by_cell_no_mut$ReadGroup, "mut"] <- "no mutation" # assign cells covered in snRNA-seq but without mutation
sample_seurat_all@meta.data[mut_in_cells$ReadGroup, "mut"] <- mut_in_cells$mutation_concat # assign cells covered in snRNA-seq but without mutation

# Keep only mutations that are in >= 3 cells of the sample
samples_vs_mut <- table(sample_seurat_all$orig.ident, sample_seurat_all$mut) # compute how many cells in sample for each mutation
cells_to_set_0 <- c() # if there are <3 cells with the mutation, annotate the cell as "unknown"
for(s in rownames(samples_vs_mut)){
  for(m in colnames(samples_vs_mut)){
    # exclude cases where we expect 2 mutations in the same cell 
    if(s %in% c("CHC2959T", "CHC2960T") & m %in% c( "chr3_41224598", "chr3_41224598,chr3_41224606", "chr3_41224606")){ # leave as such
    }else if(s=="CHC3377T" & m %in% c("chr3_41224606", "chr3_41224606,chr3_41224607", "chr3_41224607")){ # leave as such
    }else{
      if(samples_vs_mut[s, m] < 3){
        mut_to_set_0 <- m
        sample_seurat_all@meta.data[which(sample_seurat_all@meta.data$orig.ident==s & sample_seurat_all@meta.data$mut==mut_to_set_0), "mut"] <- "not known"
      } 
    }
  }    
}


##### 3. Plot #####
#------------------
pdf(file.path(resdir, paste0("UMAP_all_cells_", mutation, "_mutations.pdf")))
print(DimPlot(subset(sample_seurat_all, subset=orig.ident%in%name_sample & selection_proposition_strict != "undefined"), group.by = "orig.ident"))
print(DimPlot(subset(sample_seurat_all, subset=orig.ident%in%name_sample & selection_proposition_strict != "undefined"), group.by = "mut", pt.size=1.5,
              order=c(unique(sample_seurat_all$mut)[which(!unique(sample_seurat_all$mut) %in% c("not known", "no mutation"))], "no mutation", "not known"),
              cols=c("grey95", "grey95", colorRampPalette(c("red", "yellow", "green", "blue"))(length(unique(sample_seurat_all$mut)[which(!unique(sample_seurat_all$mut) %in% c("not known", "no mutation"))]))))) + NoLegend()
print(DimPlot(subset(sample_seurat_all, subset=orig.ident%in%name_sample & selection_proposition_strict != "undefined"), group.by = "mut", pt.size=1.5,
              order=c(unique(sample_seurat_all$mut)[which(!unique(sample_seurat_all$mut) %in% c("not known", "no mutation"))], "no mutation", "not known"),
              cols=c("grey95", "grey95", colorRampPalette(c("red", "yellow", "green", "blue"))(length(unique(sample_seurat_all$mut)[which(!unique(sample_seurat_all$mut) %in% c("not known", "no mutation"))])))))
dev.off()


