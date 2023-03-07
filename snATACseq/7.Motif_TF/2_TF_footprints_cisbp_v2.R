
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ArchR)
library(pheatmap)

addArchRGenome("hg38") # load the genome for each R session
# addArchRThreads(threads = 16) # takes automatically half of the threads available

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

archrdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/"
resdir = file.path("/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Motif_TF")

# Select motifs to look at
motifs <- c("HNF1A", "NR5A1", "CEBPB", "MAZ")


###########################
##### 0. DATA LOADING #####
###########################

project_all <- readRDS(file.path(archrdir, "Save-ArchR-Project.rds"))

### Add imputation weights
project_all <- addImputeWeights(project_all)

### Load object containing TF motif positions 
motifPositions <- readRDS(file.path(archrdir, "Annotations/Motif_cisbp_v2.00-Positions-In-Peaks.rds")) # GrangesList object
 
### Restrict to TF of interest
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE))); markerMotifs


######################################
##### 1. ADD PSEUDOBULK COVERAGE #####
######################################

# To accurately profile TF footprints, a large number of reads are required. Therefore, cells are grouped to create pseudo-bulk ATAC-seq profiles that can be then used for TF footprinting. 
project_all <- addGroupCoverages(ArchRProj = project_all, groupBy = "RNA_subtypes_hepatocytes") # quite fast
# use RNA_subtypes_hepatocytes and not RNA_subtypes because RNA_subtypes contains NAs while RNA_subtypes_hepatocytes will enable to create group for hepatocytes and other cells (that we can remove later)


####################################
##### 2. COMPUTE TF FOOTPRINTS #####
####################################

seFoot <- getFootprints(
  ArchRProj = project_all, 
  positions = motifPositions[markerMotifs], 
  groupBy = "RNA_subtypes_hepatocytes",
  useGroups = c("H", "LP", "M")
) # takes ~15 min


#################################
##### 3. PLOT TF FOOTPRINTS #####
#################################
# I had to deconstruct the function to make it work

##### 1. Arguments #####
#-----------------------
# PARAMATERS TO CHANGE
names="MAZ" # select TF name to study (must be one by one)
normMethod = c("Subtract", "Divide")[1] # method to choose for normalization
smoothWindow = 5 # smoothing window in bp to have cleaner results

seFoot = seFoot
pal = c("#E8AAD6", "#510787", "#6F3322")
flank = 250
flankNorm = 50
baseSize = 10
ArchRProj = project_all
plotName = paste0("Plot-Footprints-",normMethod)
height = 8
width = 6
addDOC = FALSE
force = FALSE
logFile = createLogFile("plotFootprints")

##### 2. Footprint function #####
#--------------------------------
name <- gsub("\\.pdf", "", plotName)
ArchRProj <- ArchR:::.validArchRProject(ArchRProj)
outDir <- file.path(getOutputDirectory(ArchRProj), "Plots")

dir.create(outDir, showWarnings = FALSE)
filename <- file.path(outDir, paste0(name, "_", names, ".pdf"))

name <- names
errorList <- list()
rowDF <- SummarizedExperiment::rowData(seFoot)
footMat <- ArchR:::.getAssay(seFoot[BiocGenerics::which(rowDF[, 2] == "footprint"), ], name)
biasMat <- ArchR:::.getAssay(seFoot[BiocGenerics::which(rowDF[, 2] == "bias"), ], name)
footDF <- rowDF[BiocGenerics::which(rowDF[, 2] == "footprint"),]
biasDF <- rowDF[BiocGenerics::which(rowDF[, 2] == "bias"),]
errorList$footMat <- footMat
errorList$biasMat <- biasMat
errorList$footDF <- footDF
errorList$biasDF <- biasDF
if (!is.null(smoothWindow)) {
  ArchR:::.logMessage("Applying smoothing window to footprint", logFile = logFile)
  footMat <- apply(footMat, 2, function(x) ArchR:::.centerRollMean(x, smoothWindow))
  biasMat <- apply(biasMat, 2, function(x) ArchR:::.centerRollMean(x, smoothWindow))
}
ArchR:::.logMessage("Normalizing by flanking regions", logFile = logFile)
idx <- which(abs(footDF$x) >= flank - flankNorm)
footMat <- t(t(footMat)/colMeans(footMat[idx, , drop = FALSE]))
biasMat <- t(t(biasMat)/colMeans(biasMat[idx, , drop = FALSE]))
errorList$footMatNorm <- footMat
errorList$biasMatNorm <- footMat
if (tolower(normMethod) == "none") {
  title <- ""
}else if (tolower(normMethod) == "subtract") {
  title <- "Tn5 Bias Subtracted\n"
  footMat <- footMat - biasMat
}else if (tolower(normMethod) == "divide") {
  title <- "Tn5 Bias Divided\n"
  footMat <- footMat/biasMat
}else {
  stop("normMethod not recognized!")
}
ArchR:::.logMessage(paste0("NormMethod = ", normMethod), logFile = logFile)
footMatMean <- ArchR:::.groupMeans(footMat, SummarizedExperiment::colData(seFoot)$Group)
footMatSd <- ArchR:::.groupSds(footMat, SummarizedExperiment::colData(seFoot)$Group)
biasMatMean <- ArchR:::.groupMeans(biasMat, SummarizedExperiment::colData(seFoot)$Group)
biasMatSd <- ArchR:::.groupSds(biasMat, SummarizedExperiment::colData(seFoot)$Group)
smoothFoot <- rowMaxs(apply(footMatMean, 2, function(x) ArchR:::.centerRollMean(x, 11)))
errorList$footMatMean <- footMatMean
errorList$footMatSd <- footMatSd
errorList$biasMatMean <- biasMatMean
errorList$biasMatSd <- biasMatSd
errorList$smoothFoot <- smoothFoot
plotIdx <- seq_len(nrow(footMatMean))
plotFootDF <- lapply(seq_len(ncol(footMatMean)), function(x) {
  data.frame(x = footDF$x, mean = footMatMean[, x], sd = footMatSd[,
                                                                   x], group = colnames(footMatMean)[x])[plotIdx, ,
                                                                                                         drop = FALSE]
}) %>% Reduce("rbind", .)
plotFootDF$group <- factor(paste0(plotFootDF$group), levels = unique(gtools::mixedsort(paste0(plotFootDF$group))))
plotBiasDF <- lapply(seq_len(ncol(biasMatMean)), function(x) {
  data.frame(x = biasDF$x, mean = biasMatMean[, x], sd = biasMatSd[,
                                                                   x], group = colnames(biasMatMean)[x])[plotIdx, ,
                                                                                                         drop = FALSE]
}) %>% Reduce("rbind", .)
plotBiasDF$group <- factor(paste0(plotBiasDF$group), levels = unique(gtools::mixedsort(paste0(plotBiasDF$group))))
errorList$plotFootDF <- plotFootDF
errorList$plotBiasDF <- plotBiasDF

if (is.null(pal)) {
  pal <- paletteDiscrete(values = gtools::mixedsort(SummarizedExperiment::colData(seFoot)$Group))
}
plotMax <- plotFootDF[order(plotFootDF$mean, decreasing = TRUE), ]
plotMax <- plotMax[abs(plotMax$x) > 20 & abs(plotMax$x) < 50, ]
plotMax <- plotMax[!duplicated(plotMax$group), ]
plotMax <- plotMax[seq_len(ceiling(nrow(plotMax)/4)), ]
plotMax$x <- 25
ggFoot <- ggplot(plotFootDF, aes(x = x, y = mean, color = group)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
  geom_line() + scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) + xlab("Distance to motif center (bp)") +
  coord_cartesian(expand = FALSE, ylim = c(quantile(plotFootDF$mean, 1e-04), 1.15 * quantile(smoothFoot, 0.999)), xlim = c(min(plotFootDF$x), max(plotFootDF$x))) +
  theme_classic()  + ggtitle(name) +
  guides(fill = FALSE) + guides(color = FALSE) + ylab(paste0(title, "Normalized Insertions")) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position = "none")
ggBias <- ggplot(plotBiasDF, aes(x = x, y = mean, color = group)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
  geom_line() + scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) + xlab("Distance to motif center (bp)") +
  coord_cartesian(expand = FALSE, ylim = c(quantile(plotBiasDF$mean, 1e-04), 1.05 * quantile(plotBiasDF$mean, 0.999)), xlim = c(min(plotBiasDF$x), max(plotBiasDF$x))) +
  theme_ArchR(baseSize = baseSize) + ylab("Tn5-Bias Normalized Insertions") +
  theme(legend.position = "bottom", legend.box.background = element_rect(color = NA))
gg <- ggAlignPlots(ggFoot, ArchR:::.ggSmallLegend(ggBias), sizes = c(2, 1), draw = FALSE)

pdf(filename, width = width, height = height, useDingbats = FALSE)
grid::grid.newpage()
grid::grid.draw(gg)
dev.off()
dev.off()


