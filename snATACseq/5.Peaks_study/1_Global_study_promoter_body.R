
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
resdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Peaks_study/Gene_features"
annot_dir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/0.Public_annotation"


###########################
##### 0. DATA LOADING #####
###########################

# Load ArchR project
project_all <- readRDS(file.path(datadir, "Save-ArchR-Project.rds"))
peak_set <- getPeakSet(project_all)

# Load annotation for peak enrichment
promoter_body <- read.table(file.path(annot_dir, "promoter_gene_body_ArchR.bed"), sep="\t", header=F); colnames(promoter_body) <- c("chr", "start", "end", "strand", "gene_feature", "gene")


############################
##### 1. GENE FEATURES #####
############################

# Look at the enrichment between all peaks and expected results for whole genome
# We cannot use the PeakAnnoEnrichment function from AchR because it does enrichment on nb of peaks and assumes all peaks are the same size which is not our case here --> might bias results

##### 1. create recap table to contain results #####
#---------------------------------------------------

recap_gene <- data.frame(feature= c("promoter", "body", "out"), prop_peaks=NA, prop_genome=NA, stringsAsFactors = F)
rownames(recap_gene) <- recap_gene$feature


##### 2. Proportion in peaks #####
#---------------------------------
# for each state compute nb of bp in peaks that intersect the state coordinates

for(f in unique(promoter_body$gene_feature)){ # only for promoter and body
  promoter_body_f <- promoter_body %>% filter(gene_feature==f)
  gr_feature <- GenomicRanges::makeGRangesFromDataFrame(promoter_body_f, keep.extra.columns = F, ignore.strand = T, starts.in.df.are.0based	= T) # create Granges
  gr_intersect <- GenomicRanges::intersect(peak_set, gr_feature) # get the intersected regions between peaks and state coordinates = the nucleotides present in the state, not the whole peaks (because some peaks are linked to > 1 chromatin state)
  # if we do findOverlaps(peak_set, gr_state) we find that some peaks have several chromatin states linked
  recap_gene[f,"prop_peaks"] <- sum(width(gr_intersect)) # how many bp in total in peaks intersect the state
}
recap_gene["out","prop_peaks"] <- sum(width(peak_set)) - recap_gene["promoter","prop_peaks"] - recap_gene["body","prop_peaks"] # the remaining nucleotides in peaks are out of genes
recap_gene$prop_peaks <- 100*recap_gene$prop_peaks/sum(recap_gene$prop_peaks) # divide by total nb of bp that are linked to one feature: since we consider bp and not peaks each bp is linked to 1 feature --> no overlap 



##### 3. Proportion in genome #####
#----------------------------------

for(f in recap_gene$feature[1:2]){ # not consider out of gene for the moment because is not included in the promoter body annotation table
  bp_genome <- sum(promoter_body[which(promoter_body$gene_feature==f),"end"] - promoter_body[which(promoter_body$gene_feature==f),"start"]) # compute total nb of bp of bins related to the gene feature
  # no "+1" when computing the bin size because we have bins like "6000-6500, 6500-6700, 6700-6900" etc
  recap_gene[f,"prop_genome"] <- 100*bp_genome/3200000000 # compute proportion by dividing by the total nb nucleotides in genome (since "out" = what is not included in the promoter body annot table)
}
recap_gene[which(recap_gene$feature=="out"),"prop_genome"] <- 100 - recap_gene[which(recap_gene$feature=="promoter"),"prop_genome"] - recap_gene[which(recap_gene$feature=="body"),"prop_genome"] 

# total sum of proportions must be =1


##### 4. Compute enrichment of peaks vs genome #####
#---------------------------------------------------

recap_gene$enrichment <- recap_gene$prop_peaks / recap_gene$prop_genome


##### 5. Generate plots #####
#----------------------------

recap_gene$feature <- factor(recap_gene$feature, levels=unique(recap_gene$feature))

pdf(file.path(resdir, "Global_study_gene_features.pdf")) 
# barplot
ggplot(recap_gene, aes(x=feature, y=prop_peaks, fill=feature))+
  geom_bar(stat="identity")+
  theme_classic()+
  ylab("% of bases in all peaks")+
  scale_fill_manual(values=c("#0C255C", "#90D1F0", "#D981AC"))
ggplot(recap_gene, aes(x=feature, y=prop_genome, fill=feature))+
  geom_bar(stat="identity")+
  theme_classic()+
  ylab("% of bases in genome")+
  scale_fill_manual(values=c("#0C255C", "#90D1F0", "#D981AC"))
ggplot(recap_gene, aes(x=feature, y=enrichment, fill=feature))+
  geom_bar(stat="identity")+
  theme_classic()+
  geom_hline(yintercept = 1, col="black", linetype="dashed")+
  ylab("Enrichment bp all peaks vs genome")+
  scale_fill_manual(values=c("#0C255C", "#90D1F0", "#D981AC"))
# target plot (ggplot)
ggplot(recap_gene, aes(x=feature, y=enrichment, fill=feature))+
  geom_bar(stat="identity")+
  theme_classic()+
  coord_polar()+
  geom_hline(yintercept = 1, col="black", linetype="dashed")+
  ylab("Enrichment bp all peaks vs genome")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_fill_manual(values=c("#0C255C", "#90D1F0", "#D981AC"))
dev.off()





