

### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)


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

resdir <- paste0("/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/1.Preprocess/Filtered_genes_gtf_files_5/")
name_samples <- c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")

# Define thresholds for QC
minGenes = 500
minUMI = 1000
maxMT = 5

###########################
##### 0. DATA LOADING #####
###########################

### Before QC
allSamps_noQC <- c()
for (sample_name in name_samples) { 
  
  load(file.path(resdir, paste0(sample_name, "/Seurat_object_noQC.RData"))) 
  sample.qc <- FetchData(sample, vars=c("nFeature_RNA","nCount_RNA","percent.mt", "ratio_UMI_genes", "log_ratio_genes_UMI")) %>% mutate(Sample_Name = sample_name)
  
  allSamps_noQC <- allSamps_noQC %>% bind_rows(sample.qc )
  print(dim(allSamps_noQC))
}

### After QC
allSamps_QC <- c()
for (sample_name in name_samples) { 
  
  load(file.path(resdir, paste0(sample_name, "/Seurat_object_QC_with_MT_genes.RData"))) 
  sample_after_qc.qc <- FetchData(sample_after_qc, vars=c("nFeature_RNA","nCount_RNA","percent.mt", "ratio_UMI_genes", "log_ratio_genes_UMI")) %>% mutate(Sample_Name = sample_name)
  
  allSamps_QC <- allSamps_QC %>% bind_rows(sample_after_qc.qc )
  print(dim(allSamps_QC))
}


####################
##### 2. PLOTS #####
####################

##### 1. Nb of genes #####
#-------------------------

ggplot(allSamps_noQC) +
  aes(x = Sample_Name, y = nFeature_RNA) +
  geom_violin(fill = "#beeafa") +
  geom_jitter(width = 0.5, size = 0.1, alpha = 0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  facet_grid(. ~ Sample_Name, scales = "free") + theme(strip.text.x = element_blank()) +
  geom_hline(aes(yintercept = minGenes), color = "red") +
  ylab("Nb of genes") +
  ggtitle("Number of genes before QC")
ggsave(file.path(resdir, "Genes_number_before_QC.pdf"), width = 10, height = 4)

ggplot(allSamps_QC) +
  aes(x = Sample_Name, y = nFeature_RNA) +
  geom_violin(fill = "#beeafa") +
  geom_jitter(width = 0.5, size = 0.1, alpha = 0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  facet_grid(. ~ Sample_Name, scales = "free") + theme(strip.text.x = element_blank()) +
  geom_hline(aes(yintercept = minGenes), color = "red") +
  ylab("Nb of genes") +
  ggtitle("Number of genes after QC")
ggsave(file.path(resdir, "Genes_number_after_QC.pdf"),  width = 10, height = 4)


##### 2. Nb of UMI #####
#-----------------------

ggplot(allSamps_noQC) +
  aes(x = Sample_Name, y = nCount_RNA) +
  geom_violin(fill = "#ccfabe") +
  geom_jitter(width = 0.5, size = 0.1, alpha = 0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  facet_grid(. ~ Sample_Name, scales = "free") + theme(strip.text.x = element_blank()) +
  geom_hline(aes(yintercept = minUMI), color = "red") +
  ylab("Nb of UMI") +
  ggtitle("Number of UMI before QC", subtitle = "y axis cut at 50000") +
  ylim(0, 50000) +
  theme( plot.title = element_text(size = 12),
         plot.subtitle = element_text(size = 8))
ggsave(file.path(resdir, "UMI_number_before_QC.pdf"), width = 10, height = 4)

ggplot(allSamps_QC) +
  aes(x = Sample_Name, y = nCount_RNA) +
  geom_violin(fill = "#ccfabe") +
  geom_jitter(width = 0.5, size = 0.1, alpha = 0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  facet_grid(. ~ Sample_Name, scales = "free") + theme(strip.text.x = element_blank()) +
  geom_hline(aes(yintercept = minUMI), color = "red") +
  ylab("Nb of UMI") +
  ggtitle("Number of UMI after QC", subtitle = "y axis cut at 50000") +
  ylim(0, 50000) +
  theme( plot.title = element_text(size = 12),
         plot.subtitle = element_text(size = 8))
ggsave(file.path(resdir, "UMI_number_after_QC.pdf"), width = 10, height = 4)


##### 3. Mitochondrial % #####
#-----------------------------

ggplot(allSamps_noQC) +
  aes(x = Sample_Name, y = percent.mt) +
  geom_violin(fill = "#fabebe") +
  geom_jitter(width = 0.5, size = 0.1, alpha = 0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  facet_grid(. ~ Sample_Name, scales = "free") + theme(strip.text.x = element_blank()) +
  geom_hline(aes(yintercept = maxMT), color = "red") +
  ylab("Mitochondrial %") +
  ggtitle("Mitochondrial % before QC")
ggsave(file.path(resdir, "Mitochondrial_percentage_before_QC.pdf"), width = 10, height = 4)

ggplot(allSamps_QC) +
  aes(x = Sample_Name, y = percent.mt) +
  geom_violin(fill = "#fabebe") +
  geom_jitter(width = 0.5, size = 0.1, alpha = 0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  facet_grid(. ~ Sample_Name, scales = "free") + theme(strip.text.x = element_blank()) +
  geom_hline(aes(yintercept = maxMT), color = "red") +
  ylab("Mitochondrial %") +
  ggtitle("Mitochondrial % after QC")
ggsave(file.path(resdir, "Mitochondrial_percentage_after_QC.pdf"), width = 10, height = 4)


