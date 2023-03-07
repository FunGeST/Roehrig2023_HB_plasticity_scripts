
### A FAIRE SUR SERVEURS LINUX avec source /home/amelie/miniconda2/etc/profile.d/conda.sh et conda activate r_env_4.1

library(dplyr)
library(tidyr)
library(ArchR)
library(pheatmap)
library(chromVAR)
library(chromVARmotifs)
library(universalmotif)
library(TFBSTools)
library(seqLogo)

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

cisbp_dir = "/mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/2.snATACseq/0.Public_annotation/Cisbp_database/" # folder containing cisbp v2.00 results from MEME software


###########################
##### 0. DATA LOADING #####
###########################

### Load cisbp motifs v2.00 (from MEME software)
ppm <- universalmotif::read_meme(file.path(cisbp_dir, "Homo_sapiens.meme"), skip = 0, readsites = FALSE, readsites.meta = FALSE) # the file is a list containing probabilitiy matrices that will require being converted to PWM


############################################
##### 1. CREATE CUSTOM MOTIF REFERENCE #####
############################################
# The idea is that ArchR takes the older version of cisbp motifs but since then many TF have been added --> we want to recover those TF


### Convert PPMatrix (probabilities) ato PWMatrix (weights, what is used by ArchR)
pwm.list <- TFBSTools::PWMatrixList() # create empty list
for(i in 1:length(ppm)){
  pwm.matrix <- universalmotif::convert_motifs(ppm[[i]], "TFBSTools-PWMatrix") # convert_motifs transforms the PPM matrix to PWM with log2(PPM* sum(TFBSTools::bg(PPM)) / TFBSTools::bg(PPM))
  pwm.list[[i]] <- pwm.matrix
}

### Add names to each slot in the list
names_list <- c()
for(i in 1:length(pwm.list)){
  pwm_i <- pwm.list[[i]]
  names_list <- c(names_list, pwm_i@ID)
}
names(pwm.list) <- names_list

### For the 83/1065 TF that have names with "()" and "_" take same nomenclature as others
names(pwm.list) <- gsub("\\_.*","",names(pwm.list))
names(pwm.list) <- gsub("[\\(]", "", gsub("[\\)]", "", names(pwm.list)))

### Save
save(pwm.list, file=file.path(cisbp_dir, "PWM_list_cisbp_v2.00.RData"))

### Print Logo  
x=pwm.list[["HNF4A"]]
x <- (2^Matrix(x))*TFBSTools::bg(x)/sum(TFBSTools::bg(x))
seqLogo::seqLogo(x)


