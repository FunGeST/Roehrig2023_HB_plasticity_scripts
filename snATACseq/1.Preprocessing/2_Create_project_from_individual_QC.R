
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ArchR)

addArchRGenome("hg38") # load the genome for each R session
# addArchRThreads(threads = 16) # takes automatically half of the threads available


######################
##### PARAMETERS #####
######################

resdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged"

names_samples <- c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")

list_arrows <- paste0("/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/", names_samples, "/ArrowFiles/", names_samples, ".arrow")


#############################
##### 1. CREATE PROJECT #####
#############################

project_all <- ArchRProject(
  ArrowFiles = list_arrows, 
  outputDirectory = resdir,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# save project

saveArchRProject(ArchRProj = project_all, outputDirectory = resdir, load = TRUE)

# here we do not filter doublets (filterDoublets in ArchR) because we want to compare first with snRNAseq --> keep the info but we will restrict the clusters later


##### Check the different pieces of information contained in the arrow project #####
#-----------------------------------------------------------------------------------

metadata_atac <- getCellColData(project_all)
save(metadata_atac, file=file.path(resdir, "metadata_all_samples.RData"))

### nb of cells per sample

nb_cells_atac <- as.data.frame(table(project_all$Sample)) 
colnames(nb_cells_atac) <- c("Sample", "Nb_cells")
write.table(nb_cells_atac, file.path(resdir, "nb_cells_per_sample.txt"), sep="\t", col.names=T, row.names=F)


### quantile of TSS enrichment on all samples

quantile_atac_TSS_enrich <- as.data.frame(quantile(project_all$TSSEnrichment)) 
colnames(quantile_atac_TSS_enrich) <- "value"
write.table(quantile_atac_TSS_enrich, file.path(resdir, "quantile_TSS_enrichment_all_sample.txt"), sep="\t", col.names=T, row.names=T)


#####################################
##### 4. DISPLAY THE QC METRICS #####
#####################################

##### All samples #####
#----------------------

### nb of fragments vs TSS enrichment

df <- getCellColData(project_all, select = c("log10(nFrags)", "TSSEnrichment"))
pdf(file.path(resdir, "TSS_by_Unique_Frags.pdf"))
ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
dev.off()


### fragment size distribution DOES NOT WORK WHEN ALL SAMPLES?

pdf(file.path(resdir, "Fragment_size_distribution.pdf"))
plotFragmentSizes(ArchRProj = project_all)
dev.off()


### other QC metrics

for(QC_name in c("TSSEnrichment", "log10(nFrags)", "ReadsInTSS", "ReadsInPromoter", "ReadsInBlacklist", "PromoterRatio", "NucleosomeRatio", "log10(nMultiFrags)", "log10(nMonoFrags)", "log10(nFrags)", "log10(nDiFrags)", "DoubletScore", "DoubletEnrichment", "BlacklistRatio")){
  pdf(paste0(resdir, "/", QC_name, ".pdf"))
  print(plotGroups(
    ArchRProj = project_all, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = QC_name,
    plotAs = "ridges"
  ))
  print(plotGroups(
    ArchRProj = project_all, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = QC_name,
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
  ))
  dev.off()
}




