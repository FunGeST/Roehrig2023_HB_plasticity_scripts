
library(dplyr)
library(Seurat)
library(ggplot2)
library(infercnv)


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

datadir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/4.Pseudo_CNV/0.Input_data"
resdir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/4.Pseudo_CNV/1.Results/All_samples_merged/"

########################
##### DATA LOADING #####
########################

### Import matrix of raw counts for all cells of all samples (cells are QC-filtered but genes are not filtered:35003 genes. The samples have been merged by Seurat but no other transformation)
counts_matrix <- geco.load(file.path(datadir, "counts_matrix_all_samples.RData"))

### Import annotation for all cells (indicates if the cells are normal = reference or other = to perform inferCNV on)
annotation <- geco.load(file.path(datadir, "cell_annotations.RData"))  # 21150 cells
rownames(annotation) <- annotation$cell
colnames(annotation) <- c("V1", "V2")

names <- rownames(annotation)
annotation <- annotation[nrow(annotation):1,]; rownames(annotation) <- names[length(names):1]
annotation$V1 <- NULL

### Import genomic position files to study the genome regions
genomic_position_file <- factoall(read.delim(file.path(datadir, "genomic_position_file_ensembl.txt"), header=F))  #60675 rows, contains ENSG ids
rownames(genomic_position_file) <- genomic_position_file$V1
genomic_position_file$V1 <- NULL
genomic_position_file[which(genomic_position_file$V2%in%c(1:22, "X", "Y", "MT")),"V2"] <- paste0("chr", genomic_position_file[which(genomic_position_file$V2%in%c(1:22, "X", "Y", "MT")),"V2"]) # add "chr" to the chromosome indexes
genomic_position_file <- genomic_position_file[which(genomic_position_file$V2%in%paste0("chr", c(1:22, "X", "Y"))),]  # keep only conventional chromosomes, 60579 rows

# we need to reorder the chromosomes because they are not perfectly ordered in the genomic_position_file:
# we do not just do arrange because it would remove the gene ordering performed previously 
genomic_position_file <- genomic_position_file[c(which(genomic_position_file$V2%in%paste0("chr", 1:7)), which(genomic_position_file$V2%in%paste0("chr", 8:9)), which(genomic_position_file$V2=="chr10"), which(genomic_position_file$V2=="chr11"),
                                                 which(genomic_position_file$V2%in%paste0("chr", 12:18)), which(genomic_position_file$V2=="chr19"), which(genomic_position_file$V2=="chr20"), which(genomic_position_file$V2=="chr21"),
                                                 which(genomic_position_file$V2=="chr22"), which(genomic_position_file$V2=="chrX"), which(genomic_position_file$V2=="chrY")),]

# load the cellranger file that gives the correspondance between gene names and ENSG
features_snRNAseq <- factoall(read.table(file.path(datadir, "features.tsv"), sep="\t", header=F))
features_snRNAseq <- features_snRNAseq[which(features_snRNAseq$V3=="Gene Expression"),] # features_snRNAseq contains both RNA and ATAC annotations, we keep only RNA
# warning: sometimes we have multiple ENSG id for the same gene name: 
ensg_to_filter <- c("ENSG00000285053", "ENSG00000261186", "ENSG00000234323", "ENSG00000261480", "ENSG00000286070", "ENSG00000269226", "ENSG00000271858") # keep the ENSG id that are given by HGNC and remove the other ENSG for the same gene 
features_snRNAseq <- features_snRNAseq[which(!features_snRNAseq$V1%in%ensg_to_filter),]
rownames(features_snRNAseq) <- features_snRNAseq$V2 # replace gene ENSG id by gene names in rownames
 

########################################
##### CNV ESTIMATION WITH INFERCNV #####
########################################

##### 1. Prepare the matrices based on genes #####
#-------------------------------------------------

# Conversion gene symbol --> gene id pour avoir concordance entre les fichiers

if(identical(rownames(counts_matrix), features_snRNAseq$V2)){rownames(counts_matrix) <- features_snRNAseq$V1} # conversion gene name --> gene id to have the same IDs for both files
genes_to_keep <- intersect(rownames(counts_matrix), rownames(genomic_position_file)) # restrict to the genes from count matrix that are also found in the ENSEMBL FILE (version 102 --> not the latest --> overlap is not 100%)

counts_matrix <- counts_matrix[genes_to_keep,] # = 34921 genes out of 35003
genomic_position_file <- genomic_position_file[genes_to_keep,]# it is not a problem if the genes are not ordered, this is performed in the first step of inferCNV (CreateInfercnvObject 
# but even if the positions in each chromosome are ordered in CreatedInfercnvObject, the chromosome order remains the same as the one in genomic_position_file --> this is why we needed to change it to have chroms in the proper order

##### 2. Create inferCNV object #####
#------------------------------------

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annotation,
                                    delim="\t",
                                    gene_order_file=genomic_position_file,
                                    ref_group_names=c("normal"),
                                    chr_exclude = c("chrMT", "chrM"))

##### 3. Run inferCNV #####
#--------------------------

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=resdir,  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster_by_group can be T because here all cells belong to the "other" annotation
                             denoise=T,
                             HMM=T
)
# inferCNV performs smoothing on windows of 101 genes

