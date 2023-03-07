
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ArchR)

addArchRGenome("hg38") # load the genome for each R session
# addArchRThreads(threads = 16) # takes automatically half of the threads available


######################
##### PARAMETERS #####
######################

resdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses"

name_sample <- c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")[1]

if(!dir.exists(resdir)){dir.create(file.path(resdir, name_sample))}
input_bam <- paste0("/mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/", name_sample, "/outs/atac_fragments.tsv.gz")

if(name_sample=="CHC2959N"){
  TSS_threshold=7;nFrags_threshold=10^3.5
}else if(name_sample=="CHC2959T"){
  TSS_threshold=6;nFrags_threshold=10^3.7
}else if(name_sample=="CHC2960T"){
  TSS_threshold=8;nFrags_threshold=10^3.5
}else if(name_sample=="CHC3133T"){
  TSS_threshold=4;nFrags_threshold=10^3
}else if(name_sample=="CHC3377N"){
  TSS_threshold=4;nFrags_threshold=10^3
}else if(name_sample=="CHC3377T"){
  TSS_threshold=6.5;nFrags_threshold=10^3.5
}else if(name_sample=="CHC3610T"){
  TSS_threshold=4;nFrags_threshold=10^3.7
}else if(name_sample=="CHC3662T"){
  TSS_threshold=5;nFrags_threshold=10^3.5
}


#################################
##### 1. CREATE ARROW FILES #####
#################################

# WARNING: before check in home directory that there are no arrow files already present

#createarrowfiles takes ~45 min 
ArrowFiles <- createArrowFiles(
  inputFiles = input_bam,
  sampleNames = name_sample,
  QCDir=file.path(resdir, name_sample),
  minTSS = TSS_threshold, 
  minFrags = nFrags_threshold, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# ArrowFiles is a character vector of arrow file path .arrow, the files are automatically saved in my home
# it seems that it is not possible to run several createArrowFiles in different R sessions at the same time --> always resulted in hd5 problem


#############################
##### 2. INFER DOUBLETS #####
#############################

# takes ~5 min per sample
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)


###################################
##### 3. CREATE ARROW PROJECT #####
###################################

project_sample1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = file.path(resdir, name_sample),
  copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# save project

saveArchRProject(ArchRProj = project_sample1, outputDirectory = file.path(resdir, name_sample), load = TRUE)


##### Check the different pieces of information contained in the arrow project #####
#-----------------------------------------------------------------------------------

metadata_atac <- getCellColData(project_sample1) # extract metadata
save(metadata_atac, file=file.path(resdir, name_sample, "metadata.RData"))


#####################################
##### 4. DISPLAY THE QC METRICS #####
#####################################

### nb of fragments vs TSS enrichment

df <- getCellColData(project_sample1, select = c("log10(nFrags)", "TSSEnrichment"))
pdf(file.path(resdir, name_sample, "TSS_by_Unique_Frags.pdf"))
ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
)
dev.off()


### fragment size distribution

pdf(file.path(resdir, name_sample, "Fragment_size_distribution.pdf"))
plotFragmentSizes(ArchRProj = project_sample1)
dev.off()


### other QC metrics

for(QC_name in c("TSSEnrichment", "log10(nFrags)", "ReadsInTSS", "ReadsInPromoter", "ReadsInBlacklist", "PromoterRatio", "NucleosomeRatio", "log10(nMultiFrags)", "log10(nMonoFrags)", "log10(nFrags)", "log10(nDiFrags)", "DoubletScore", "DoubletEnrichment", "BlacklistRatio")){
  pdf(paste0(resdir, "/", name_sample, "/", QC_name, ".pdf"))
  print(plotGroups(
    ArchRProj = project_sample1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = QC_name,
    plotAs = "ridges"
  ))
  print(plotGroups(
    ArchRProj = project_sample1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = QC_name,
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
  ))
  dev.off()
}







