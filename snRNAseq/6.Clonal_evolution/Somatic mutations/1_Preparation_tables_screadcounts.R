

library(dplyr)
library(ggplot2)
library(Seurat)


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

name_samples = c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3377T", "CHC3610T", "CHC3662T") # we will investigate tumor samples
name_samples1= c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3610T") # those samples were sequenced in WGS in a first project (2020)
name_samples2 = c("CHC3377T", "CHC3662T") # those samples were sequenced in WGS in a second project (2022)

resdir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data"


###########################
##### 0. DATA LOADING #####
###########################

### Load mutation table

# WGS table for CHC2959T, CHC2960T, CHC3133T, CHC3610T
vcf_all1 <- geco.load(file.path(resdir, "Tumors_203s_151p_vcf_filtered.RData")) # the file originates from hepato partage\GENOMIC_DATA\WXS\WGS\RData\Somatic_mutations\Tumors
# WGS table for CHC3377T, CHC3662T
vcf_all2 <- geco.load(file.path(resdir, "vcf_3662T_3377T_vcf_all_with_CCF.RData")) # the file originates from hepato partage\GEPELIN\Mutation_tables\New_samples
 

###############################################
##### 1. CREATE MUTATION TABLE PER SAMPLE #####
###############################################

for(name_sample in name_samples){
  if(name_sample%in%name_samples1){ 
    vcf_all. <- vcf_all1 %>% filter(Sample==name_sample) %>% select(CHROM, POS_hg38, REF, ALT) %>% arrange(CHROM, POS_hg38) # keep only mutations for the sample, keep positions, sort by position
    colnames(vcf_all.) <- c("CHROM", "POS", "REF", "ALT")
    vcf_all. <- vcf_all.[which(!is.na(vcf_all.$POS)),] # remove the possible few NAs that could have been introduced between hg19 and hg38
    write.table(vcf_all., paste0(resdir, "/mutation_table_", name_sample, ".tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  }else if(name_sample%in%name_samples2){ 
    vcf_all. <- vcf_all2 %>% filter(Sample==name_sample) %>% select(CHROM, POS_hg38, REF, ALT) %>% arrange(CHROM, POS_hg38) # keep only mutations for the sample, keep positions, sort by position
    colnames(vcf_all.) <- c("CHROM", "POS", "REF", "ALT")
    vcf_all. <- vcf_all.[which(!is.na(vcf_all.$POS)),] # remove the possible few NAs that could have been introduced between hg19 and hg38
    write.table(vcf_all., paste0(resdir, "/mutation_table_", name_sample, ".tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  }
}


###########################################################
##### 2. CREATE MUTATION TABLE FOR 2959T/60T COMBINED #####
###########################################################

vcf_59_60 <- vcf_all1 %>% filter(Sample%in%c("CHC2959T", "CHC2960T")) %>% select(CHROM, POS_hg38, REF, ALT) %>% arrange(CHROM, POS_hg38) # keep only mutations for the sample, keep positions, sort by position
colnames(vcf_59_60) <- c("CHROM", "POS", "REF", "ALT")
vcf_59_60 <- vcf_59_60[which(!is.na(vcf_59_60$POS)),] # remove the possible few NAs that could have been introduced between hg19 and hg38
vcf_59_60 <- unique(vcf_59_60) # remove duplicate rows (due to common mutations between CHC2959T and CHC2960T)
write.table(vcf_59_60, paste0(resdir, "/mutation_table_CHC2959T_CHC2960T_union.tsv"), sep="\t", col.names=T, row.names=F, quote=F)


