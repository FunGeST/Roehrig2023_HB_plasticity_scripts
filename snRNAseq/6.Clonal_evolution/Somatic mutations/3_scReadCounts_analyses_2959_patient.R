

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



######################
##### PARAMETERS #####
######################

name_sample = c("CHC2959T", "CHC2960T")

datadir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/"
screadcounts_dir <- file.path("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/Somatic_mutations")
seurat_dir_tumor <- file.path("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict")
seurat_dir_tumor_cells_2959 <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/CHC2959T/"
seurat_dir_tumor_cells_2960 <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict/CHC2960T/"


###########################
##### 0. DATA LOADING #####
###########################

##### 1. WGS info #####
#----------------------

vcf_wgs <- geco.load(file.path(datadir, "Merged_pediatric_129s_89p_vcf_all_with_CCF_and_SBS_proba.RData")) #last version, contains all tumors but is in hg19
vcf_wgs <- vcf_wgs[which(vcf_wgs$CHROM == vcf_wgs$CHROM_hg38),] # remove potential liftover misconversions

vcf_wgs <- vcf_wgs %>% filter(Sample %in% c(name_sample)) 
vcf_wgs$mutation_hg19 <- paste(vcf_wgs$CHROM, vcf_wgs$POS, sep="_") # to compare hg38 positions to hg19. CHROM/POS are in hg19, CHROM_hg38/POS_hg38 are in hg38}
vcf_wgs$mutation_hg38 <- paste(vcf_wgs$CHROM_hg38, vcf_wgs$POS_hg38, sep="_") # to compare hg38 positions to hg19. CHROM/POS are in hg19, CHROM_hg38/POS_hg38 are in hg38}

vcf_clusters <- geco.load(file.path(datadir, "/vcf_with_clusters_all_HB_multi_WGS.RData")) # contains mutation clusters annotated by Théo, hg19
vcf_clusters <- vcf_clusters[grep("2959",vcf_clusters$Cluster),] # restrict to clusters of patient 2959 
vcf_clusters$mutation_hg19 <- paste(vcf_clusters$CHROM, vcf_clusters$POS, sep="_") # the table is in hg19

vcf_wgs <- vcf_wgs %>% filter(mutation_hg19 %in% vcf_clusters$mutation_hg19)  # restrict the vcf to mutations in Théo clusters
vcf_wgs <- merge(x=vcf_wgs, y=vcf_clusters[, c("mutation_hg19", "Cluster")], by.x="mutation_hg19", by.y="mutation_hg19", all.x=T) # add info of mutation clusters


##### 2. Screadcounts results #####
#----------------------------------

### Create merged table containing screadcount results for the sample ###
# initialize
screadcounts_results <- read.table(file.path(screadcounts_dir, name_sample[1], "Screadcounts_output", paste0("scReadCounts_output_xf_25_", gsub("CHC", "", name_sample[1]), "_mut.tsv")), sep="\t", header=T)
screadcounts_results$ReadGroup <- gsub("-1", "", paste0(name_sample[1], "_", screadcounts_results$ReadGroup))
screadcounts_results$sample <- rep(name_sample[1], times=nrow(screadcounts_results))
# loop 
for(n in name_sample){
  if(n==name_sample[1]){
    screadcounts_results. <- read.table(file.path(screadcounts_dir, n, "Screadcounts_output", paste0("scReadCounts_output_xf_25_", gsub("CHC", "", name_sample[2]), "_mut.tsv")), sep="\t", header=T)
    screadcounts_results.$ReadGroup <- gsub("-1", "", paste0(n, "_", screadcounts_results.$ReadGroup))
    screadcounts_results.$sample <- rep(n, times=nrow(screadcounts_results.))
    screadcounts_results <- rbind(screadcounts_results, screadcounts_results.)
  }else if(n==name_sample[2]){
    screadcounts_results. <- read.table(file.path(screadcounts_dir, n, "Screadcounts_output", paste0("scReadCounts_output_xf_25_", gsub("CHC", "", name_sample[1]), "_mut.tsv")), sep="\t", header=T)
    screadcounts_results.$ReadGroup <- gsub("-1", "", paste0(n, "_", screadcounts_results.$ReadGroup))
    screadcounts_results.$sample <- rep(n, times=nrow(screadcounts_results.))
    screadcounts_results <- rbind(screadcounts_results, screadcounts_results.)
    screadcounts_results. <- read.table(file.path(screadcounts_dir, n, "Screadcounts_output", paste0("scReadCounts_output_xf_25_", gsub("CHC", "", name_sample[2]), "_mut.tsv")), sep="\t", header=T)
    screadcounts_results.$ReadGroup <- gsub("-1", "", paste0(n, "_", screadcounts_results.$ReadGroup))
    screadcounts_results.$sample <- rep(n, times=nrow(screadcounts_results.))
    screadcounts_results <- rbind(screadcounts_results, screadcounts_results.)
  }
}

screadcounts_results$CHROM <- paste0("chr", screadcounts_results$CHROM)
screadcounts_results$mutation <- paste(screadcounts_results$CHROM, screadcounts_results$POS, sep="_") # create a column to easily retrieve the mutation position
length(unique(screadcounts_results$mutation)) 

### Annotate mutations clusters
screadcounts_results$mutation_cluster <- vcf_wgs[match(screadcounts_results$mutation, vcf_wgs$mutation_hg38),"Cluster"]

# Keep mutations with "PASS" filter
screadcounts_results$PASS <- ifelse(screadcounts_results$mutation %in% vcf_wgs[which(vcf_wgs$FILTER=="PASS"),"mutation_hg38"], "yes", "no")

save(screadcounts_results, file=file.path(screadcounts_dir,"Patient_2959/screadcounts_results_mutation_info.RData"))


##### 3. Load Seurat objects #####
#---------------------------------

### Load seurat object for the sample
sample_tumor_cells_2959 <- geco.load(paste0(seurat_dir_tumor_cells_2959, "/Seurat_object_analyses.RData"))
sample_tumor_cells_2960 <- geco.load(paste0(seurat_dir_tumor_cells_2960, "/Seurat_object_analyses.RData"))
# merge the samples
sample_tumor_cells <- merge(sample_tumor_cells_2959, sample_tumor_cells_2960) # no need to add cell ID, already present
# re normalize the merge + PCA + UMAP
sample_tumor_cells <- SCTransform(sample_tumor_cells, verbose = FALSE, variable.features.n=3000) 
sample_tumor_cells <- RunPCA(sample_tumor_cells, features = VariableFeatures(sample_tumor_cells), npcs=100)
sample_tumor_cells <- RunUMAP(sample_tumor_cells, dims=1:40, verbose=F, umap.method = "uwot") # use 40 dimensions (will work fine in any case)


#####################################################
##### 1. VISUALIZE MUTATIONS ON INDIVIDUAL UMAP #####
#####################################################
# Only the tumor cells of the samples of interest

# Choose which mutations to filter on
filter_by_cluster = c("2959_C1", "2959_C2", "2959_C3", "2959_C5", "2959_C8", "2959_C10")[2]
if(filter_by_cluster=="2959_C10"){col_cluster="black"}else if(filter_by_cluster=="2959_C2"){col_cluster="#F6493C"}else if(filter_by_cluster=="2959_C5"){col_cluster="#B79F00"
}else if(filter_by_cluster=="2959_C3"){col_cluster="#0070C0"}else if(filter_by_cluster=="2959_C8"){col_cluster="#7030A0"}else if(filter_by_cluster=="2959_C1"){col_cluster="#04BA38"}

##### 1. Concatenate mutation info per cell #####
#------------------------------------------------
recap_by_cell <- screadcounts_results %>% 
  filter(mutation_cluster %in% filter_by_cluster) %>%
  filter(PASS == "yes") %>% # keep only mutations that pass the criteria in WGS (to remove possible noise)
  group_by(ReadGroup) %>% 
  summarise(SNV_total_counts = sum(SNVCount), Ref_total_counts = sum(RefCount), VAF = SNV_total_counts / (SNV_total_counts + Ref_total_counts))

recap_by_cell <- recap_by_cell[which(recap_by_cell$SNV_total_counts!=0 | recap_by_cell$Ref_total_counts!=0 ),] # remove rows that have 0 SNV and 0 Ref counts = not covered (~1% of the rows)
recap_by_cell <- recap_by_cell[which(recap_by_cell$ReadGroup%in%colnames(sample_tumor_cells)),] # screadcounts takes cells without QC -->  some of them are not in our postQC seurat object

##### 2. Add SNVcount info per cell #####
#----------------------------------------
cells_mut_3_sup <- recap_by_cell[which(recap_by_cell$SNV_total_counts>=3),"ReadGroup"]$ReadGroup 
cells_mut_2 <- recap_by_cell[which(recap_by_cell$SNV_total_counts==2),"ReadGroup"]$ReadGroup
cells_mut_1 <- recap_by_cell[which(recap_by_cell$SNV_total_counts==1),"ReadGroup"]$ReadGroup 
cells_no_mut <- recap_by_cell[which(recap_by_cell$SNV_total_counts==0),"ReadGroup"]$ReadGroup 

sample_tumor_cells$mut <- "not known"
sample_tumor_cells@meta.data[cells_mut_1,"mut"] <- "mutation, 1 SNV read"
sample_tumor_cells@meta.data[cells_mut_2,"mut"] <- "mutation, 2 SNV reads"
sample_tumor_cells@meta.data[cells_mut_3_sup,"mut"] <- "mutation, >=3 SNV reads"
sample_tumor_cells@meta.data[cells_no_mut,"mut"] <- "no mutation"

##### 3. Add Refcount info per cell #####
#----------------------------------------
cells_ref_3_sup <- recap_by_cell[which(recap_by_cell$Ref_total_counts>=3),"ReadGroup"]$ReadGroup 
cells_ref_2 <- recap_by_cell[which(recap_by_cell$Ref_total_counts==2),"ReadGroup"]$ReadGroup
cells_ref_1 <- recap_by_cell[which(recap_by_cell$Ref_total_counts==1),"ReadGroup"]$ReadGroup 
cells_no_ref <- recap_by_cell[which(recap_by_cell$Ref_total_counts==0),"ReadGroup"]$ReadGroup 

sample_tumor_cells$ref <- "not known"
sample_tumor_cells@meta.data[cells_ref_1,"ref"] <- "ref nucleotide, 1 SNV read"
sample_tumor_cells@meta.data[cells_ref_2,"ref"] <- "ref nucleotide, 2 SNV reads"
sample_tumor_cells@meta.data[cells_ref_3_sup,"ref"] <- "ref nucleotide, >=3 SNV reads"
sample_tumor_cells@meta.data[cells_no_ref,"ref"] <- "no ref nucleotide"

##### 3. Plot #####
#------------------
pdf(file.path(screadcounts_dir, "Patient_2959", paste0("cluster_", filter_by_cluster, "_mutations_CHC2959T_CHC2960T.pdf")))
# remove doublets for better visualisation
print(DimPlot(sample_tumor_cells, group.by = "mut", order=c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known"),
              pt.size=1.5, cols=c("grey", "grey20", "yellow", "orange", "red")))
print(DimPlot(sample_tumor_cells, group.by = "mut", order=c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known"),
              pt.size=1.5, cols=c("grey", "grey20", "yellow", "orange", "red")) + NoLegend())
print(DimPlot(sample_tumor_cells, group.by = "mut", order=c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known"),
              pt.size=1.7, cols=c("grey", "grey", col_cluster, col_cluster, col_cluster)))
print(DimPlot(sample_tumor_cells, group.by = "mut", order=c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known"),
              pt.size=1.7, cols=c("grey", "grey", col_cluster,col_cluster, col_cluster)) + NoLegend())
print(DimPlot(sample_tumor_cells, group.by = "ref", order=c("ref nucleotide, >=3 SNV reads", "ref nucleotide, 2 SNV reads", "ref nucleotide, 1 SNV read", "no ref nucleotide", "not known"),
              pt.size=1.5, cols=c("grey", "grey20", "yellow", "orange", "red")))
print(DimPlot(sample_tumor_cells, group.by = "ref", order=c("ref nucleotide, >=3 SNV reads", "ref nucleotide, 2 SNV reads", "ref nucleotide, 1 SNV read", "no ref nucleotide", "not known"),
              pt.size=1.5, cols=c("grey", "grey20", "yellow", "orange", "red")) + NoLegend())

dev.off()

