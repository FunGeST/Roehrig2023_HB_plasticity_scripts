
library(dplyr)
library(Seurat)
library(ggplot2)
library(infercnv)
library(magrittr)

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


resdir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/4.Pseudo_CNV/0.Input_data/"

##############################
##### TABLE PREPARATION ######
##############################

##### Determination of reference cells for CNV = normal hepatocyte cells in CHC3981N #####
#-----------------------------------------------------------------------------------------

CHC2959N <- geco.load("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/2.Basic_analyses/SCTransform/No_regression/CHC2959N/With_MT_genes/Seurat_object_analyses.RData")
normal_cells_for_CNV <- rownames(CHC2959N@meta.data)[which(CHC2959N@meta.data$seurat_clusters%in%c(0,2))] # cluster dans CHC2959N qui d'après nous correspond aux hépatocytes = cellules normales. J'ai enlevé le cluster 7 car il y a un fort %MT et ce sont peut être des doublets 
normal_cells_for_CNV <- paste0("CHC2959N_", normal_cells_for_CNV)
save(normal_cells_for_CNV, file=file.path(resdir, "normal_cells_for_CNV_CHC2959N_clusters_0_2.RData"))


##### Count matrix tables #####
#------------------------------

# InferCNV demande que les cellules soient déjà filtrées mais on garde tous les gènes

list_samples <- list()
i=1

for(s in c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")){
  sample_noQC <- geco.load(paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/1.Preprocess/Filtered_genes_gtf_files_5/", s, "/Seurat_object_noQC.RData"))
  sample_QC <- geco.load(paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/1.Preprocess/Filtered_genes_gtf_files_5/", s, "/Seurat_object_QC_with_MT_genes.RData"))
  sample <- sample_noQC[,intersect(colnames(sample_noQC), colnames(sample_QC))]
  list_samples[[i]] <- sample
  i <- i+1
}

CHC2959N <- list_samples[[1]]
CHC2959T <- list_samples[[2]]
CHC2960T <- list_samples[[3]]
CHC3133T <- list_samples[[4]]
CHC3377N <- list_samples[[5]]
CHC3377T <- list_samples[[6]]
CHC3610T <- list_samples[[7]]
CHC3662T <- list_samples[[8]]

all_samples <- merge(CHC2959N, y = c(CHC2959T, CHC2960T, CHC3133T, CHC3377N, CHC3377T, CHC3610T, CHC3662T), add.cell.ids = c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T"), project = "HB_multiome")
counts_matrix <-  GetAssayData(all_samples, slot="counts")

save(counts_matrix, file=file.path(resdir, "counts_matrix_all_samples.RData"))


##### Annotation table #####
#---------------------------

annotation <- factoall(data.frame(cell=colnames(counts_matrix), type=NA))

annotation <- annotation %>% mutate(type=case_when(cell %in% normal_cells_for_CNV ~ "normal",
                                                   !cell %in% normal_cells_for_CNV~ "other"))

save(annotation, file=file.path(resdir, "cell_annotations.RData"))


##### Gene ordering table #####
#------------------------------

### Conversion GTF to txt with script from inferCNV (reference genome = ENSEMBL)
python /mnt/MDWorks6/Workdir_amelie/scRNAseq/HB_tumors/Analyses/2.Seurat_analyses/SCTransform/Without_regression/pseudo_CNV/gtf_to_position_file.py Homo_sapiens.GRCh38.102.gtf genomic_position_file_ensembl.txt
# Homo_sapiens.GRCh38.102.gtf is found on the ENSEMBL website

# if I want to redo with cellranger genes.gtf:
# python /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/4.Pseudo_CNV/0.Input_data/gtf_to_position_file.py /mnt/MD1200_5/Workdir_amelie/Cellranger_arc_new_reference_package/refdata_cellranger_arc_new_genes_gtf_file_5/genes.gtf genomic_position_file_ensembl.txt






