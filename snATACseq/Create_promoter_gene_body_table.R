
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ArchR)
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


addArchRGenome("hg38") # load the genome for each R session
# addArchRThreads(threads = 16) # takes automatically half of the threads available

######################
##### PARAMETERS #####
######################

datadir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/"


###########################
##### 0. DATA LOADING #####
###########################

# Load ArchR project
project_all <- readRDS(file.path(datadir, "Save-ArchR-Project.rds"))
peak_set <- getPeakSet(project_all) # set of peaks in the ArchR project


##################################################
##### 1. CREATE PROMOTER/GENE BODY/OUT TABLE #####
##################################################

##### 1. Extract the genes from ArchR project #####
#--------------------------------------------------

gene_gr <- getGenes(project_all) # 25017 genes, granges object
# getGenes is the function that is used by ArchR in addGeneScoreMatrix to identify gene TSS


##### 2. Create table from gene granges object #####
#---------------------------------------------------

gene_annot_df <- as.data.frame(matrix(NA, ncol=5, nrow=length(gene_gr)))
colnames(gene_annot_df) <- c("gene", "chr", "start", "end", "strand")

gene_annot_df$gene <- gene_gr$symbol
gene_annot_df$chr <- as.character(seqnames(gene_gr))
gene_annot_df$start <- ranges(gene_gr)@start
gene_annot_df$end <- ranges(gene_gr)@start + ranges(gene_gr)@width -1
gene_annot_df$strand <- as.character(strand(gene_gr))
gene_annot_df <- gene_annot_df[which(!is.na(gene_annot_df$gene)),] # remove the NA genes = ~0.2%


##### 3. Create promoter/gene body/intergenic table #####
#--------------------------------------------------------

# define promoter as being TSS +/- 0/2000bp but in ArchR seems they take 5000bp ? (cf extendUpstream parameter in addGeneScoreMatrix)
for(i in 1:nrow(gene_annot_df))
{
  #print(i)
  tmp <- rbind(gene_annot_df[i,],gene_annot_df[i,])
  tmp$gene_feature <- c("promoter","body")
  if(gene_annot_df[i,"strand"]=="+"){
    tmp[1,c("start","end")] <- tmp[1,"start"]+c(-2000,0)
    tmp[2,c("start","end")] <- c(gene_annot_df[i,"start"],gene_annot_df[i,"end"])
  }else if(gene_annot_df[i,"strand"]=="-"){
    tmp[1,c("start","end")] <- tmp[1,"end"]+c(0,2000)
    tmp[2,c("start","end")] <- c(gene_annot_df[i,"start"],gene_annot_df[i,"end"])
  }
  
  if(i==1){new <- tmp}else{new <- rbind(new,tmp)}
}
prom_body_annot_df <- new

# convert to bed file
prom_body_annot_df <- prom_body_annot_df[,c("chr", "start", "end", "strand", "gene_feature", "gene")]


##### 4. Add the "out of gene" info ONLY ON THE EXPECTED PEAKS #####
#-------------------------------------------------------------------

# since we want the "out of gene" info for peaks in additioni to promoter and gene body, but cannot create an "out of gene" bed file, we look at which peaks have no intersection with promoter and gene body and 
# assign those peaks to "out of gene"

promoter_body <- prom_body_annot_df

# Create GRangesList to provide for peak enrichment: we have to separate each gene feature in its own granges object
promoter_body_gr <- list()
for(f in unique(promoter_body$gene_feature)){
  promoter_body. <- promoter_body %>% filter(gene_feature == f)
  gr_gene_feature <- GenomicRanges::makeGRangesFromDataFrame(promoter_body., keep.extra.columns = F, ignore.strand = F, starts.in.df.are.0based	= T)
  promoter_body_gr[f] <- gr_gene_feature
}

promoter_body_gr_list <- c(promoter=promoter_body_gr[["promoter"]], body=promoter_body_gr[["body"]])

# Add annotation to project
project_all <- addPeakAnnotations(ArchRProj = project_all, regions = promoter_body_gr_list, name = "gene_features", force=T)

# Add the column "out of gene in the RDS generated files before enrichment

get_matches <- getMatches(project_all, name="gene_features")
get_matches_mat <- get_matches@assays@data$matches
idx_out <- which(get_matches_mat[,"promoter"]==0 & get_matches_mat[,"body"]==0) # recover indexes of the peaks that do not match neither promoter nor gene body = out of gene
peaks_out_gr <- peak_set[idx_out] # restrict peak set to peaks out of genes

gene_annot_df_out <- as.data.frame(matrix(NA, ncol=6, nrow=length(peaks_out_gr)))
colnames(gene_annot_df_out) <- colnames(prom_body_annot_df)

gene_annot_df_out$gene <- "out_of_gene"
gene_annot_df_out$chr <- as.character(seqnames(peaks_out_gr))
gene_annot_df_out$start <- ranges(peaks_out_gr)@start
gene_annot_df_out$end <- ranges(peaks_out_gr)@start + ranges(peaks_out_gr)@width -1
gene_annot_df_out$strand <- as.character(strand(peaks_out_gr))
gene_annot_df_out$gene_feature <- "out"

# add gene_annot_df_out to promoter_body table

promoter_body <- rbind(promoter_body, gene_annot_df_out)

##### SAVE #####

write.table(prom_body_annot_df, file="/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/0.Public_annotation/promoter_gene_body_ArchR.bed", sep="\t", col.names=F, row.names=F, quote=F)



