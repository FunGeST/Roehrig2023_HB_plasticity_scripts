
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ggplot2)
#library(geco.utils)
library(Matrix)
library(Seurat)
library(stringr)
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


######################
##### PARAMETERS #####
######################

sample_name="CHC2959N"

input.directory <- paste0("/mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/", sample_name, "/outs/filtered_feature_bc_matrix")
output.directory <- paste0("/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/1.Preprocess/Filtered_genes_gtf_files_5/", sample_name)
if(!dir.exists(output.directory)){dir.create(output.directory)}



########################
##### DATA LOADING #####
########################

#if(packageVersion("Seurat") == "3.0.2"){sample.data <- Read10X(data.dir=input.directory, gene.column=1) 
#}else if(packageVersion("Seurat") == "3.2.3"){sample.data <- Read10X(data.dir=input.directory, gene.column=1, strip.suffix=T)}
sample.data <- Read10X(data.dir=input.directory, gene.column=1, strip.suffix=T) # list of 2 matrices. We have to load the matrix with ENSG and not gene names because there are 5 gene names duplicates. strip.suffix gets rid of the "-1" suffix in cellranger barcodes
sample.data.RNA <- sample.data[[1]] # 35010 genes

cellranger_features <- factoall(read.table(file.path(input.directory, "features.tsv.gz"), sep="\t", header=F))
ensg_to_filter <- c("ENSG00000285053", "ENSG00000261186", "ENSG00000234323", "ENSG00000261480", "ENSG00000286070", "ENSG00000269226", "ENSG00000271858") # keep the ENSG id that are given by HGNC and remove the other ENSG for the same gene 

sample.data.RNA <- sample.data.RNA[which(!rownames(sample.data.RNA)%in%ensg_to_filter),] # remove the genes in expression dataset
rownames(cellranger_features) <- cellranger_features$V1 
rownames(sample.data.RNA) <- cellranger_features[rownames(sample.data.RNA),"V2"] # change gene ID to gene name

sample <- CreateSeuratObject(counts = sample.data.RNA, project=sample_name, min.cells=0, min.features=0) # Create Seurat object from expression. 35003 genes. do not set filter on cells or genes, this will be done after


#######################
##### QC ANALYSIS #####
#######################

##### Observation of main features, plots before QC #####
#--------------------------------------------------------

sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern="^MT-") # add mitochondrial percentage 

sample$ratio_UMI_genes <- sample$nCount_RNA/sample$nFeature_RNA
sample$log_ratio_genes_UMI <- log10(sample$nFeature_RNA)/log10(sample$nCount_RNA)

sample.qc <- FetchData(sample, vars=c("nFeature_RNA","nCount_RNA","percent.mt", "ratio_UMI_genes", "log_ratio_genes_UMI")) # retrieve data as data frame for histograms

pdf(file.path(output.directory, "Before_QC.pdf"))

ggplot(data=sample.qc)+
  geom_histogram(aes(x=ratio_UMI_genes),bins=100)+
  xlab("UMI counts/nb of genes")
ggplot(data=sample.qc)+
  geom_histogram(aes(x=log_ratio_genes_UMI),bins=100)+
  xlab("log10(nb of genes)/log10(UMI counts)")+
  geom_vline(xintercept = 0.8)

FeatureScatter(sample, feature1="nCount_RNA", feature2="percent.mt")+
  xlab("UMI counts")+
  ylab("Mitochondrial %")+
  geom_hline(yintercept = 5)+
  geom_vline(xintercept = 1000)

FeatureScatter(sample, feature1="nFeature_RNA", feature2="percent.mt")+
  xlab("Nb of genes")+
  ylab("Mitochondrial %")+
  geom_hline(yintercept = 5)+
  geom_vline(xintercept = 500)

FeatureScatter(sample, feature1="nCount_RNA", feature2="nFeature_RNA")+
  xlab("UMI counts")+
  ylab("Nb of genes")+
  geom_hline(yintercept = 500)+
  geom_vline(xintercept = 1000)
ggplot(data=sample.qc, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "orange") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_hline(yintercept = 500)+
  geom_vline(xintercept = 1000)

VlnPlot(sample, features = c("nFeature_RNA"), ncol=3)+ 
  ylab("Nb of genes")+
  geom_hline(yintercept = 500)
VlnPlot(sample, features = c("nCount_RNA"), ncol=3)+
  ylab("UMI counts")+
  geom_hline(yintercept = 1000)
VlnPlot(sample, features = c("percent.mt"), ncol=3)+
  ylab("Mitochondrial %")+
  geom_hline(yintercept = 5)

ggplot(data=sample.qc)+
  geom_histogram(aes(x=nCount_RNA),bins=100)+
  xlab("UMI counts")+
  geom_vline(xintercept = 1000)
ggplot(data=sample.qc)+
  geom_histogram(aes(x=log10(nCount_RNA)),bins=150)+
  xlab("log10 UMI counts")
ggplot(data=sample.qc)+
  geom_histogram(aes(x=nFeature_RNA),bins=100)+
  xlab("Nb of genes")+
  geom_vline(xintercept = 500)
ggplot(data=sample.qc)+
  geom_histogram(aes(x=log10(nFeature_RNA)),bins=150)+
  xlab("log10 Nb of genes")
sample.qc %>% dplyr::mutate(rank_nCount_RNA = dense_rank(nCount_RNA)) %>%
  ggplot() + geom_point(aes(x=rank_nCount_RNA, y=nCount_RNA))+xlab("UMI counts rank")+ylab("UMI counts")+scale_y_log10()

dev.off()


##### Filtering #####
#--------------------

# Apply QC thresholds
sample_after_qc <- subset(sample, subset = percent.mt < 5 & nCount_RNA > 1000 & nFeature_RNA > 500)

# Keep genes with >=1 transcript in >= 3 non-outlier cells
sample_after_qc_mat <- as.matrix(GetAssayData(object = sample_after_qc, slot = "counts"))
selection <- apply(sample_after_qc_mat, 1, function(x) length(x[x>=1])>=3)
sample_after_qc <- sample_after_qc[selection,] # restrict to genes in the Seurat object


##### Plots after QC ##### 
#-------------------------

sample_after_qc[["percent.mt"]] <- PercentageFeatureSet(sample_after_qc, pattern="^MT-") # add mitochondrial percentage (not provided by CellRanger)
sample.qc <- FetchData(sample_after_qc, vars=c("nFeature_RNA","nCount_RNA","percent.mt", "ratio_UMI_genes", "log_ratio_genes_UMI")) # retrieve data as data frame for histograms

pdf(file.path(output.directory, "After_QC.pdf"))
ggplot(data=sample.qc)+
  geom_histogram(aes(x=ratio_UMI_genes),bins=100)+
  xlab("UMI counts/nb of genes")
ggplot(data=sample.qc)+
  geom_histogram(aes(x=log_ratio_genes_UMI),bins=100)+
  xlab("log10(nb of genes)/log10(UMI counts)")+
  geom_vline(xintercept = 0.8)

FeatureScatter(sample_after_qc, feature1="nCount_RNA", feature2="percent.mt")+
  xlab("UMI counts")+
  ylab("Mitochondrial %")+
  geom_hline(yintercept = 5)+
  geom_vline(xintercept = 1000)

FeatureScatter(sample_after_qc, feature1="nFeature_RNA", feature2="percent.mt")+
  xlab("Nb of genes")+
  ylab("Mitochondrial %")+
  geom_hline(yintercept = 5)+
  geom_vline(xintercept = 1700)

FeatureScatter(sample_after_qc, feature1="nCount_RNA", feature2="nFeature_RNA")+
  xlab("UMI counts")+
  ylab("Nb of genes")+
  geom_hline(yintercept = 500)+
  geom_vline(xintercept = 1000)
ggplot(data=sample.qc, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "orange") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_hline(yintercept = 500)+
  geom_vline(xintercept = 1000)

VlnPlot(sample_after_qc, features = c("nFeature_RNA"), ncol=3)+ # do.return transforms into ggplot2 graph
  ylab("Nb of genes")+
  geom_hline(yintercept = 500)
VlnPlot(sample_after_qc, features = c("nCount_RNA"), ncol=3)+
  ylab("UMI counts")+
  geom_hline(yintercept = 1000)
VlnPlot(sample_after_qc, features = c("percent.mt"))+
  ylab("Mitochondrial %")+
  geom_hline(yintercept = 5)

ggplot(data=sample.qc)+
  geom_histogram(aes(x=nCount_RNA),bins=100)+
  xlab("UMI counts")+
  geom_vline(xintercept = 1000)
ggplot(data=sample.qc)+
  geom_histogram(aes(x=log10(nCount_RNA)),bins=150)+
  xlab("log10 UMI counts")
ggplot(data=sample.qc)+
  geom_histogram(aes(x=nFeature_RNA),bins=100)+
  xlab("Nb of genes")+
  geom_vline(xintercept = 500)
ggplot(data=sample.qc)+
  geom_histogram(aes(x=log10(nFeature_RNA)),bins=150)+
  xlab("log10 Nb of genes")
sample.qc %>% dplyr::mutate(rank_nCount_RNA = dense_rank(nCount_RNA)) %>%
  ggplot() + geom_point(aes(x=rank_nCount_RNA, y=nCount_RNA))+xlab("UMI counts rank")+ylab("UMI counts")+scale_y_log10()

dev.off()



##### Final steps and save #####
#-------------------------------

#Remove MT genes

MT_genes <- grep("MT-", rownames(sample_after_qc), value=T)
sample_after_qc_no_MT <- sample_after_qc[which(!rownames(sample_after_qc)%in%MT_genes),]


save(sample_after_qc, file=file.path(output.directory, "Seurat_object_QC_with_MT_genes.RData"))
save(sample_after_qc_no_MT, file=file.path(output.directory, "Seurat_object_QC_without_MT_genes.RData"))
save(sample, file=file.path(output.directory, "Seurat_object_noQC.RData")) 


