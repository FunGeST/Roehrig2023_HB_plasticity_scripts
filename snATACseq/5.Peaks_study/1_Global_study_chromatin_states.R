
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
resdir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/Final_analyses/All_samples_merged/Peaks_study/Chromatin_states"
annot_dir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/0.Public_annotation"


###########################
##### 0. DATA LOADING #####
###########################

# Load ArchR project
project_all <- readRDS(file.path(datadir, "Save-ArchR-Project.rds"))
peak_set <- getPeakSet(project_all) # retrieve peakset granges

# Load annotation for peak enrichment
chrom_states_adult_liver <- read.table(file.path(annot_dir, "E066_18_core_K27ac_mnemonics.bed"), sep="\t", header=F); colnames(chrom_states_adult_liver) <- c("chr", "start", "end", "state") # E066 = normal adult liver 
# chrom_states_HepG2 <- read.table(file.path(annot_dir, "E118_18_core_K27ac_mnemonics.bed"), sep="\t", header=F); colnames(chrom_states_HepG2) <- c("chr", "start", "end", "state")
# here we only use the adult liver data because cell line HepG2 brings too many biases 


###############################
##### 1. CHROMATIN STATES #####
###############################

# Look at the enrichment between all peaks and expected results for whole genome
# We cannot use the PeakAnnoEnrichment function from AchR because it does enrichment on nb of peaks and assumes all peaks are the same size which is not our case here --> might bias results

##### 1. create recap table to contain results #####
#---------------------------------------------------

recap_cst18 <- data.frame(state= c("1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD", "5_Tx", "6_TxWk", "7_EnhG1", "8_EnhG2", "9_EnhA1", "10_EnhA2", "11_EnhWk", "12_ZNF/Rpts", "13_Het",
                                   "14_TssBiv", "15_EnhBiv", "16_ReprPC", "17_ReprPCWk", "18_Quies"), prop_peaks=NA, prop_genome=NA, stringsAsFactors = F)
rownames(recap_cst18) <- recap_cst18$state

##### 2. Proportion in peaks #####
#---------------------------------
# for each state compute nb of bp in peaks that intersect the state coordinates

for(s in unique(chrom_states_adult_liver$state)){
  chrom_states_adult_liver_s <- chrom_states_adult_liver %>% filter(state==s)
  gr_state <- GenomicRanges::makeGRangesFromDataFrame(chrom_states_adult_liver_s, keep.extra.columns = F, ignore.strand = T, starts.in.df.are.0based	= T) # create Granges
  gr_intersect <- GenomicRanges::intersect(peak_set, gr_state) # get the intersected regions between peaks and state coordinates = the nucleotides present in the state, not the whole peaks (because some peaks are linked to > 1 chromatin state)
  # if we do findOverlaps(peak_set, gr_state) we find that some peaks have several chromatin states linked
  recap_cst18[s,"prop_peaks"] <- sum(width(gr_intersect)) # how many bp in total in peaks intersect the state
}
recap_cst18$prop_peaks <- 100*recap_cst18$prop_peaks/sum(recap_cst18$prop_peaks) # divide by total nb of bp that are linked to one state: since we consider bp and not peaks each bp is linked to 1 chromatin state --> no overlap 


##### 3. Proportion in genome #####
#----------------------------------

for(s in recap_cst18$state){
  bp_genome <- sum(chrom_states_adult_liver[which(chrom_states_adult_liver$state==s),"end"] - chrom_states_adult_liver[which(chrom_states_adult_liver$state==s),"start"]) # compute total nb of bp of bins related to the state
  # no "+1" when computing the bin size because we have bins like "6000-6500, 6500-6700, 6700-6900" etc
  recap_cst18[s,"prop_genome"] <- 100*bp_genome/(sum(chrom_states_adult_liver[,"end"] - chrom_states_adult_liver[,"start"])) # compute proportion by dividing by the total nb of bp in all bins
}
# total sum of proportions must be =1
# total nb of nuvleotides covered with cst18 annot file: 2865575419 = 90% of total nb of nucleotides in genome

##### 4. Compute enrichment of peaks vs genome #####
#---------------------------------------------------

recap_cst18$enrichment <- recap_cst18$prop_peaks / recap_cst18$prop_genome


##### 5. Generate plots #####
#----------------------------

recap_cst18$state <- factor(recap_cst18$state, levels=unique(recap_cst18$state))

pdf(file.path(resdir, "Global_study_adult_liver_cst18.pdf")) 
# barplot
ggplot(recap_cst18, aes(x=state, y=prop_peaks, fill=state))+
  geom_bar(stat="identity")+
  theme_classic()+
  ylab("% of bases in all peaks")+
  scale_fill_manual(values=c("darkred", "red"," red1", "red2", "forestgreen", "darkgreen", "darkolivegreen1", "darkolivegreen2", "orange", "orange1", "yellow2", "turquoise", "lightslateblue", "palevioletred4", "olivedrab3" ,"grey50", "grey80", "grey95"))
ggplot(recap_cst18, aes(x=state, y=prop_genome, fill=state))+
  geom_bar(stat="identity")+
  theme_classic()+
  ylab("% of bases in genome")+
  scale_fill_manual(values=c("darkred", "red"," red1", "red2", "forestgreen", "darkgreen", "darkolivegreen1", "darkolivegreen2", "orange", "orange1", "yellow2", "turquoise", "lightslateblue", "palevioletred4", "olivedrab3" ,"grey50", "grey80", "grey95"))
ggplot(recap_cst18, aes(x=state, y=enrichment, fill=state))+
  geom_bar(stat="identity")+
  theme_classic()+
  geom_hline(yintercept = 1, col="black", linetype="dashed")+
  ylab("Enrichment bp all peaks vs genome")+
  scale_fill_manual(values=c("darkred", "red"," red1", "red2", "forestgreen", "darkgreen", "darkolivegreen1", "darkolivegreen2", "orange", "orange1", "yellow2", "turquoise", "lightslateblue", "palevioletred4", "olivedrab3" ,"grey50", "grey80", "grey95"))
# target plot (ggplot)
ggplot(recap_cst18, aes(x=state, y=enrichment, fill=state))+
  geom_bar(stat="identity")+
  theme_classic()+
  coord_polar()+
  ylab("Enrichment bp all peaks vs genome")+
  geom_hline(yintercept = 1, col="black", linetype="dashed")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_fill_manual(values=c("darkred", "red"," red1", "red2", "forestgreen", "darkgreen", "darkolivegreen1", "darkolivegreen2", "orange", "orange1", "yellow2", "turquoise", "lightslateblue", "palevioletred4", "olivedrab3" ,"grey50", "grey80", "grey95"))
dev.off()



