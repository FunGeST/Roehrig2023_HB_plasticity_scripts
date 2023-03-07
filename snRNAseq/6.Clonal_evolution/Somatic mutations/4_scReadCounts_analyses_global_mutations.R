

library(dplyr)
library(ggplot2)
library(infercnv)
library(phytools)
library(dendextend)
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

scale_fill_custom <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(c("red", "orange", "yellow", "grey20", "grey"), c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known")), 
    ...
  )
}


######################
##### PARAMETERS #####
######################

sample_to_study = "CHC3662T"
all_samples <- c("CHC2959N", "CHC2959T", "CHC2960T", "CHC3133T", "CHC3377N", "CHC3377T", "CHC3610T", "CHC3662T")
name_sample = c(sample_to_study, all_samples[which(all_samples!=sample_to_study)]) # the sammple under study must be first

datadir <- "F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/"
screadcounts_dir <- file.path("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/Somatic_mutations")
seurat_dir_tumor <- file.path("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/Tumor/2.Tumor_selection_strict")
seurat_dir_all <- file.path("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/6.Clean_cluster_definition_tumor_immune/")
CNV_dir <- paste0("F:/Amelie-Datas/Pediatric tumors/3.Multiome ATAC RNA/2.snRNAseq/9.Clonal_evolution/", name_sample[1], "/Final_clusters")


###########################
##### 0. DATA LOADING #####
###########################

##### 1. WGS info #####
#----------------------

vcf_wgs_hg19 <- geco.load(file.path(datadir, "Merged_pediatric_129s_89p_vcf_all_with_CCF_and_SBS_proba.RData")) #last version, contains all tumors but is in hg19
if(sample_to_study %in% c("CHC3662T", "CHC3377T")){vcf_wgs_hg38 <- geco.load(file.path(datadir, "/vcf_3662T_3377T_vcf_all_with_CCF.RData")) # this file contains the hg38 info + clonality etc
}else if(sample_to_study %in% c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3610T")){vcf_wgs_hg38 <- geco.load(file.path(datadir, "/Tumors_203s_151p_vcf_filtered.RData"))} # this file contains the hg38 info + clonality etc
head(colnames(vcf_wgs_hg38))
if(sample_to_study %in% c("CHC3662T", "CHC3377T")){vcf_wgs_hg38 <- arrange(vcf_wgs_hg38, CHROM, POS) # CHROM and POS correspond to hg19
}else if(sample_to_study %in% c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3610T")){vcf_wgs_hg38 <- arrange(vcf_wgs_hg38, CHROM, POS_hg19)}

vcf_wgs_hg19 <- vcf_wgs_hg19 %>% filter(Sample %in% c(name_sample[1])) 
vcf_wgs_hg19$mutation_hg19 <- paste(vcf_wgs_hg19$CHROM, vcf_wgs_hg19$POS, sep="_")
rownames(vcf_wgs_hg19) <- vcf_wgs_hg19$mutation_hg19

vcf_wgs_hg38 <- vcf_wgs_hg38 %>% filter(Sample%in% c(name_sample[1])) 
if(sample_to_study %in% c("CHC3662T", "CHC3377T")){vcf_wgs_hg38$mutation_hg19 <- paste(vcf_wgs_hg38$CHROM, vcf_wgs_hg38$POS, sep="_") # to compare hg38 positions to hg19. CHROM/POS are in hg19, CHROM_hg38/POS_hg38 are in hg38
}else if(sample_to_study %in% c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3610T")){vcf_wgs_hg38$mutation_hg19 <- paste(vcf_wgs_hg38$CHROM, vcf_wgs_hg38$POS_hg19, sep="_")} # to compare hg38 positions to hg19. CHROM/POS are in hg19, CHROM_hg38/POS_hg38 are in hg38}
rownames(vcf_wgs_hg38) <- vcf_wgs_hg38$mutation_hg19

if(identical(vcf_wgs_hg19$mutation_hg19, vcf_wgs_hg38$mutation_hg19)){# must be TRUE, the positions must correspond between hg19 and hg38
  if(sample_to_study %in% c("CHC3662T", "CHC3377T")){vcf_wgs_hg19$mutation_hg38 <- paste(vcf_wgs_hg38$CHROM_hg38, vcf_wgs_hg38$POS_hg38, sep="_")
  }else if(sample_to_study %in% c("CHC2959T", "CHC2960T", "CHC3133T", "CHC3610T")){vcf_wgs_hg19$mutation_hg38 <- paste(vcf_wgs_hg38$CHROM, vcf_wgs_hg38$POS_hg38, sep="_")}
} 

if(sample_to_study %in% c("CHC3662T", "CHC3377T")){vcf_wgs_hg38 <- vcf_wgs_hg38[which(vcf_wgs_hg38$CHROM == vcf_wgs_hg38$CHROM_hg38),]} # remove liftover misconversions = when chromosome hg19 is different from hg38 = ~1%
vcf_wgs_hg19 <- vcf_wgs_hg19[rownames(vcf_wgs_hg38),]


##### 2. snMultiome info #####
#-----------------------------

### Create merged table containing screadcount results for the sample ###
# initialize
screadcounts_results <- read.table(file.path(screadcounts_dir, name_sample[1], "Screadcounts_output", paste0("scReadCounts_output_xf_25_", gsub("CHC", "", name_sample[1]), "_mut.tsv")), sep="\t", header=T)
screadcounts_results$ReadGroup <- gsub("-1", "", paste0(name_sample[1], "_", screadcounts_results$ReadGroup))
screadcounts_results$sample <- rep(name_sample[1], times=nrow(screadcounts_results))
# loop on the 2 NT samples
for(n in name_sample[2:length(name_sample)]){
  # Global output
  screadcounts_results. <- read.table(file.path(screadcounts_dir, n, "Screadcounts_output", paste0("scReadCounts_output_xf_25_", gsub("CHC", "", name_sample[1]), "_mut.tsv")), sep="\t", header=T)
  screadcounts_results.$ReadGroup <- gsub("-1", "", paste0(n, "_", screadcounts_results.$ReadGroup))
  screadcounts_results.$sample <- rep(n, times=nrow(screadcounts_results.))
  screadcounts_results <- rbind(screadcounts_results, screadcounts_results.)
}

# Add ATAC screadcounts results for more reliability
#for(n in name_sample){
#  # Global output
#  screadcounts_results. <- read.table(file.path(screadcounts_dir, n, "Screadcounts_output", paste0("scReadCounts_output_", gsub("CHC", "", name_sample[1]), "_mut_ATAC.tsv")), sep="\t", header=T)
#  screadcounts_results.$ReadGroup <- gsub("-1", "", paste0(n, "_", screadcounts_results.$ReadGroup))
#  screadcounts_results.$sample <- rep(n, times=nrow(screadcounts_results.))
#  screadcounts_results <- rbind(screadcounts_results, screadcounts_results.)
#}
screadcounts_results$CHROM <- paste0("chr", screadcounts_results$CHROM)
screadcounts_results$mutation <- paste(screadcounts_results$CHROM, screadcounts_results$POS, sep="_") # create a column to easily retrieve the mutation position
length(unique(screadcounts_results$mutation)) 

### Annotate clonal/subclonal mutations
clonal_mut <- unique(vcf_wgs_hg19[which(vcf_wgs_hg19$Clonality=="clonal"),"mutation_hg38"]) 
subclonal_mut <- unique(vcf_wgs_hg19[which(vcf_wgs_hg19$Clonality=="subclonal"),"mutation_hg38"])

### Annotate cisplatin mutations
cisplat_mut_0.9 <- unique(vcf_wgs_hg19[which(vcf_wgs_hg19$SBS35.prob_by_clonality>=0.9),"mutation_hg38"]) 
cisplat_mut_0.75 <- unique(vcf_wgs_hg19[which(vcf_wgs_hg19$SBS35.prob_by_clonality>=0.75),"mutation_hg38"]) 
#cisplat_mut <- unique(vcf_wgs_hg19[which(vcf_wgs_hg19$SBS_cat3%in%c("CT_CCC", "CT_CCT")),"mutation_hg38"]) #  mutations probably related to SBS35

screadcounts_results$Clonality <- ifelse(screadcounts_results$mutation %in% clonal_mut, "clonal",
                                         ifelse(screadcounts_results$mutation %in% subclonal_mut, "subclonal", NA)) 
screadcounts_results$SBS35 <- ifelse(screadcounts_results$mutation %in% cisplat_mut_0.9, 0.9,
                                         ifelse(screadcounts_results$mutation %in% cisplat_mut_0.75, 0.75, 0))
screadcounts_results$PASS <- ifelse(screadcounts_results$mutation %in% vcf_wgs_hg19[which(vcf_wgs_hg19$FILTER=="PASS"),"mutation_hg38"], "yes", "no")
table(screadcounts_results$Clonality, screadcounts_results$SBS35) 

save(screadcounts_results, file=file.path(screadcounts_dir, name_sample[1], "Analyses/screadcounts_results_mutation_info.RData"))

### Load Seurat objects ###
sample_seurat_indiv <- geco.load(file.path(seurat_dir_tumor, name_sample[1], "Seurat_object_analyses.RData")) # tumor cells of the sample under study
sample_seurat_tumor <- geco.load(file.path(seurat_dir_tumor, "Seurat_object_analyses.RData")) # tumor cells
sample_seurat_all <- geco.load(file.path(seurat_dir_all, "Seurat_object_analyses_all_cells.RData")) # all cells

# add tumor info in sample_seurat_all
sample_seurat_all$tumor_cells <- "non tumor"
sample_seurat_all@meta.data[colnames(sample_seurat_tumor),"tumor_cells"] <- "tumor"

### Load CNV clusters ###
CNV_clusters <- geco.load(file.path(CNV_dir, "final_CNV_clusters.RData"))
rownames(CNV_clusters) <- CNV_clusters$cells
sample_seurat_indiv$CNV_clusters <- NA; sample_seurat_indiv@meta.data[CNV_clusters$cells,"CNV_clusters"] <- CNV_clusters$final_CNV_cluster


###################################################
##### 1. VISUALIZE MUTATIONS ON COMPLETE UMAP #####
###################################################

# Choose which mutations to filter on
filter_by_clonality = c("clonal", "subclonal")[c(1)]
filter_by_SBS35 = c(0, 0.75, 0.9)[c(1)] # 0 = not SBS35 related, 0.75 = SBS35 related with proba of 0.75, 0.9 = SBS35 related with proba of 0.9

##### 1. Concatenate mutation info per cell #####
#------------------------------------------------
recap_by_cell <- screadcounts_results %>% 
  filter(Clonality %in% filter_by_clonality) %>%
  filter(SBS35 >= filter_by_SBS35) %>%
  filter(PASS == "yes") %>% # keep only mutations that pass the criteria in WGS (to remove possible noise)
  group_by(ReadGroup) %>% 
  summarise(SNV_total_counts = sum(SNVCount), Ref_total_counts = sum(RefCount), VAF = SNV_total_counts / (SNV_total_counts + Ref_total_counts))

recap_by_cell <- recap_by_cell[which(recap_by_cell$SNV_total_counts!=0 | recap_by_cell$Ref_total_counts!=0 ),] # remove rows that have 0 SNV and 0 Ref counts = not covered (~1% of the rows)
recap_by_cell <- recap_by_cell[which(recap_by_cell$ReadGroup%in%colnames(sample_seurat_all)),] # screadcounts takes cells without QC -->  some of them are not in our postQC seurat object

##### 2. Concatenate info per mutation #####
#-------------------------------------------

#recap_by_mut <- data.frame(mutation=unique(screadcounts_results$mutation), CHC2959N=NA, CHC2959T=NA, CHC2960T=NA, CHC3133T=NA, CHC3377N=NA, CHC3377T=NA, CHC3610T=NA, CHC3662T=NA, stringsAsFactors = F)
#for(mut in recap_by_mut$mutation){
#  for(s in colnames(recap_by_mut)[2:ncol(recap_by_mut)]){
#    recap_by_mut[which(recap_by_mut$mutation==mut), s] <- sum(screadcounts_results[which(screadcounts_results$mutation==mut & screadcounts_results$sample==s),"SNVCount"])
#  }
#}
#recap_by_mut <- recap_by_mut[which(recap_by_mut$mutation %in% screadcounts_results[which(screadcounts_results$Clonality=="clonal" & screadcounts_results$PASS=="yes"),"mutation"]),]
#sum_recap_by_mut <- apply(recap_by_mut[,2:ncol(recap_by_mut)],1, sum); recap_by_mut <- recap_by_mut[which(sum_recap_by_mut!=0),]

# Mutations detected only in CHC3662T with sufficient reads
#mut_good <- c("chr5_149505378", "chr13_41045182", "chr18_7594287", "chr2_62009577", "chr1_43555655", "chr18_62768071", "chr10_122387057", "chr20_37218294", "chr7_7226298", "chr9_128521392", "chr20_3632893", "chr3_179624125", "chr5_108962951", "chr2_31638407", "chr2_205321965")
#mut_good <- recap_by_mut[which(recap_by_mut$CHC3662T > 0 & recap_by_mut$CHC2959N == 0 & recap_by_mut$CHC2959T==0 & recap_by_mut$CHC2960T==0 &
 #                                recap_by_mut$CHC3133T==0 & recap_by_mut$CHC3377N==0 & recap_by_mut$CHC3377T==0 & recap_by_mut$CHC3610T==0),"mutation"]
#recap_mut_good <- recap_by_mut[which(recap_by_mut$mutation %in% mut_good),]
#vcf_mut_good <- vcf_wgs_hg19[which(vcf_wgs_hg19$mutation_hg38 %in% mut_good),]
# Mutations detected more in other cells than CHC3662T
#mut_bad <- c("chr1_227289880", "chr19_9043971", "chr15_28761141", "chr8_22374906", "chr1_61082452", "chr19_7267235", "chr2_233697515", "chr1_213230724", "chr19_23804487", "chr1_2158111", "chr3_44510582", "chr7_27769709", "chr1_2158111", "chr12_58880157", "chr12_81183233", "chr6_138660192", "chr10_32378520")
#mut_bad <- recap_by_mut[which(recap_by_mut$CHC3662T == 0 & (recap_by_mut$CHC2959N > 0 | recap_by_mut$CHC2959T>0 | recap_by_mut$CHC2960T>0 |
#                                recap_by_mut$CHC3133T==0 | recap_by_mut$CHC3377N>0 | recap_by_mut$CHC3377T>0 | recap_by_mut$CHC3610T>0)),"mutation"]
#recap_mut_bad <- recap_by_mut[which(recap_by_mut$mutation %in% mut_bad),]
#vcf_mut_bad <- vcf_wgs_hg19[which(vcf_wgs_hg19$mutation_hg38 %in% mut_bad),]


##### 2. Add SNVcount info per cell #####
#----------------------------------------
cells_mut_3_sup <- recap_by_cell[which(recap_by_cell$SNV_total_counts>=3),"ReadGroup"]$ReadGroup 
cells_mut_2 <- recap_by_cell[which(recap_by_cell$SNV_total_counts==2),"ReadGroup"]$ReadGroup
cells_mut_1 <- recap_by_cell[which(recap_by_cell$SNV_total_counts==1),"ReadGroup"]$ReadGroup 
cells_no_mut <- recap_by_cell[which(recap_by_cell$SNV_total_counts==0),"ReadGroup"]$ReadGroup 

sample_seurat_all$mut <- "not known"
sample_seurat_all@meta.data[cells_mut_1,"mut"] <- "mutation, 1 SNV read"
sample_seurat_all@meta.data[cells_mut_2,"mut"] <- "mutation, 2 SNV reads"
sample_seurat_all@meta.data[cells_mut_3_sup,"mut"] <- "mutation, >=3 SNV reads"
sample_seurat_all@meta.data[cells_no_mut,"mut"] <- "no mutation"


##### 3. Add Refcount info per cell #####
#----------------------------------------
cells_ref_3_sup <- recap_by_cell[which(recap_by_cell$Ref_total_counts>=3),"ReadGroup"]$ReadGroup 
cells_ref_2 <- recap_by_cell[which(recap_by_cell$Ref_total_counts==2),"ReadGroup"]$ReadGroup
cells_ref_1 <- recap_by_cell[which(recap_by_cell$Ref_total_counts==1),"ReadGroup"]$ReadGroup 
cells_no_ref <- recap_by_cell[which(recap_by_cell$Ref_total_counts==0),"ReadGroup"]$ReadGroup 

sample_seurat_all$ref <- "not known"
sample_seurat_all@meta.data[cells_ref_1,"ref"] <- "ref nucleotide, 1 SNV read"
sample_seurat_all@meta.data[cells_ref_2,"ref"] <- "ref nucleotide, 2 SNV reads"
sample_seurat_all@meta.data[cells_ref_3_sup,"ref"] <- "ref nucleotide, >=3 SNV reads"
sample_seurat_all@meta.data[cells_no_ref,"ref"] <- "no ref nucleotide"


##### 4. Plot #####
#------------------
### Plot enrichment in sample tumor cells vs other cells
sample_seurat_all$expected_cells <- "unexpected cells"
sample_seurat_all@meta.data[rownames(sample_seurat_indiv@meta.data),"expected_cells"] <- "expected cells" # annotate cells in which we expect to see the mutation = tumor cells from the sample under study
df <- reshape2::melt(table(sample_seurat_all$mut, sample_seurat_all$expected_cells)); df$Var1 <- factor(df$Var1, levels=c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known"))
# Compute significance with chisquare test
sample_seurat_all$mut_red <-  sample_seurat_all$mut
sample_seurat_all@meta.data[grep("mutation, ", sample_seurat_all$mut_red),"mut_red"] <- "mutation >= 1 SNV read"  # reduce the info to "mutated" or "not mutated"
sample_seurat_indiv$CNV_clusters <- factor(sample_seurat_indiv$CNV_clusters, levels=unique(sample_seurat_indiv$CNV_clusters))
tt <- table(sample_seurat_all$mut_red, sample_seurat_all$expected_cells)
pval <- signif(chisq.test(tt[which(rownames(tt)!="not known"),])$p.value, digits=3)

# Prepare file names
if(length(filter_by_clonality) > 1){filter_by_clonality <- paste(filter_by_clonality, collapse="_")}else{filter_by_clonality <- paste0(filter_by_clonality, "_")}
if(filter_by_SBS35==0){filter_by_SBS35 <- ""}else{filter_by_SBS35 <- paste0("_SBS35_", filter_by_SBS35, "_")}
pdf(file.path(screadcounts_dir, name_sample[1], "Analyses", paste0(filter_by_clonality, filter_by_SBS35, "mutations_", paste(name_sample, collapse="_"), ".pdf")))
# remove doublets for better visualisation
print(DimPlot(subset(sample_seurat_all, subset=orig.ident%in%name_sample & selection_proposition_strict != "undefined"), group.by = "orig.ident"))
print(DimPlot(subset(sample_seurat_all, subset=orig.ident%in%name_sample & selection_proposition_strict != "undefined"),
              group.by = "mut", order=c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known"),
              pt.size=0.5, cols=c("grey", "grey20", "yellow", "orange", "red")))
print(DimPlot(subset(sample_seurat_all, subset=orig.ident%in%name_sample & selection_proposition_strict != "undefined"),
              group.by = "mut", order=c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known"),
              pt.size=0.5, cols=c("grey", "grey20", "yellow", "orange", "red")) + NoLegend())
print(DimPlot(subset(sample_seurat_all, subset=orig.ident%in%name_sample & selection_proposition_strict != "undefined"),
              group.by = "ref", order=c("ref nucleotide, >=3 SNV reads", "ref nucleotide, 2 SNV reads", "ref nucleotide, 1 SNV read", "no ref nucleotide", "not known"),
              pt.size=0.5, cols=c("grey", "grey20", "yellow", "orange", "red")))
print(DimPlot(subset(sample_seurat_all, subset=orig.ident%in%name_sample & selection_proposition_strict != "undefined"),
              group.by = "ref", order=c("ref nucleotide, >=3 SNV reads", "ref nucleotide, 2 SNV reads", "ref nucleotide, 1 SNV read", "no ref nucleotide", "not known"),
              pt.size=0.5, cols=c("grey", "grey20", "yellow", "orange", "red")) + NoLegend())
ggplot(df[which(df$Var1!="not known"),], aes(x=Var2, y=value, fill=Var1))+
  geom_bar(stat="identity", position="fill")+
  theme_classic()+
  scale_fill_custom()+
  xlab("") + ylab("% of cells") + ggtitle(paste0("mutation detection in covered cells, pval=", pval))
dev.off()


#####################################################
##### 2. VISUALIZE MUTATIONS ON INDIVIDUAL UMAP #####
#####################################################
# Only the tumor cells of the sample of interest

# Choose which mutations to filter on
filter_by_clonality = c("clonal", "subclonal")[c(1)]
filter_by_SBS35 = c(0, 0.75, 0.9)[c(2)] # 0 = not SBS35 related, 0.75 = SBS35 related with proba of 0.75, 0.9 = SBS35 related with proba of 0.9

##### 1. Concatenate mutation info per cell #####
#------------------------------------------------
recap_by_cell <- screadcounts_results %>% 
  filter(Clonality %in% filter_by_clonality) %>%
  filter(SBS35 >= filter_by_SBS35) %>%
  filter(PASS == "yes") %>% # keep only mutations that pass the criteria in WGS (to remove possible noise)
  group_by(ReadGroup) %>% 
  summarise(SNV_total_counts = sum(SNVCount), Ref_total_counts = sum(RefCount), VAF = SNV_total_counts / (SNV_total_counts + Ref_total_counts))

recap_by_cell <- recap_by_cell[which(recap_by_cell$SNV_total_counts!=0 | recap_by_cell$Ref_total_counts!=0 ),] # remove rows that have 0 SNV and 0 Ref counts = not covered (~1% of the rows)
recap_by_cell <- recap_by_cell[which(recap_by_cell$ReadGroup%in%colnames(sample_seurat_indiv)),] # screadcounts takes cells without QC -->  some of them are not in our postQC seurat object


##### 2. Add SNVcount info per cell #####
#----------------------------------------
cells_mut_3_sup <- recap_by_cell[which(recap_by_cell$SNV_total_counts>=3),"ReadGroup"]$ReadGroup 
cells_mut_2 <- recap_by_cell[which(recap_by_cell$SNV_total_counts==2),"ReadGroup"]$ReadGroup
cells_mut_1 <- recap_by_cell[which(recap_by_cell$SNV_total_counts==1),"ReadGroup"]$ReadGroup 
cells_no_mut <- recap_by_cell[which(recap_by_cell$SNV_total_counts==0),"ReadGroup"]$ReadGroup 

sample_seurat_indiv$mut <- "not known"
sample_seurat_indiv@meta.data[cells_mut_1,"mut"] <- "mutation, 1 SNV read"
sample_seurat_indiv@meta.data[cells_mut_2,"mut"] <- "mutation, 2 SNV reads"
sample_seurat_indiv@meta.data[cells_mut_3_sup,"mut"] <- "mutation, >=3 SNV reads"
sample_seurat_indiv@meta.data[cells_no_mut,"mut"] <- "no mutation"


##### 3. Add Refcount info per cell #####
#----------------------------------------
cells_ref_3_sup <- recap_by_cell[which(recap_by_cell$Ref_total_counts>=3),"ReadGroup"]$ReadGroup 
cells_ref_2 <- recap_by_cell[which(recap_by_cell$Ref_total_counts==2),"ReadGroup"]$ReadGroup
cells_ref_1 <- recap_by_cell[which(recap_by_cell$Ref_total_counts==1),"ReadGroup"]$ReadGroup 
cells_no_ref <- recap_by_cell[which(recap_by_cell$Ref_total_counts==0),"ReadGroup"]$ReadGroup 

sample_seurat_indiv$ref <- "not known"
sample_seurat_indiv@meta.data[cells_ref_1,"ref"] <- "ref nucleotide, 1 SNV read"
sample_seurat_indiv@meta.data[cells_ref_2,"ref"] <- "ref nucleotide, 2 SNV reads"
sample_seurat_indiv@meta.data[cells_ref_3_sup,"ref"] <- "ref nucleotide, >=3 SNV reads"
sample_seurat_indiv@meta.data[cells_no_ref,"ref"] <- "no ref nucleotide"

##### 3. Plot #####
#------------------
### Plot enrichment in CNV clusters
df <- reshape2::melt(table(sample_seurat_indiv$mut, sample_seurat_indiv$CNV_clusters)); df$Var1 <- factor(df$Var1, levels=c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known"))
# Compute significance with chisquare test
sample_seurat_indiv$mut_red <-  sample_seurat_indiv$mut
sample_seurat_indiv@meta.data[grep("mutation, ", sample_seurat_indiv$mut_red),"mut_red"] <- "mutation >= 1 SNV read"  # reduce the info to "mutated" or "not mutated"
sample_seurat_indiv$CNV_clusters <- factor(sample_seurat_indiv$CNV_clusters, levels=unique(sample_seurat_indiv$CNV_clusters))
tt <- table(sample_seurat_indiv$mut_red, sample_seurat_indiv$CNV_clusters)
# this time compute cluster vs cluster chisquare
pval <- c(); i=1; names_pval <- c()
for(c1 in levels(sample_seurat_indiv$CNV_clusters)[-length(levels(sample_seurat_indiv$CNV_clusters))]){
  i=i+1
  for(c2 in levels(sample_seurat_indiv$CNV_clusters)[i:length(levels(sample_seurat_indiv$CNV_clusters))]){
    pval <- c(pval, signif(chisq.test(tt[which(rownames(tt)!="not known"),c(c1,c2)])$p.value, digits=3))
    names_pval <- c(names_pval, paste0(c1, " VS ", c2))
  }
}
names(pval) <- names_pval

# Prepare file names
if(length(filter_by_clonality) > 1){filter_by_clonality <- paste(filter_by_clonality, collapse="_")}else{filter_by_clonality <- paste0(filter_by_clonality, "_")}
if(filter_by_SBS35==0){filter_by_SBS35 <- ""}else{filter_by_SBS35 <- paste0("_SBS35_", filter_by_SBS35, "_")}
pdf(file.path(screadcounts_dir, name_sample[1], "Analyses", paste0(filter_by_clonality, filter_by_SBS35, "mutations_", name_sample[1], ".pdf")))
# remove doublets for better visualisation
print(DimPlot(sample_seurat_indiv, group.by = "CNV_clusters", pt.size=1.2))
print(DimPlot(sample_seurat_indiv, group.by = "mut", order=c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known"),
              pt.size=1.5, cols=c("grey", "grey20", "yellow", "orange", "red")))
print(DimPlot(sample_seurat_indiv, group.by = "mut", order=c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known"),
              pt.size=1.5, cols=c("grey", "grey20", "yellow", "orange", "red")) + NoLegend())
print(DimPlot(sample_seurat_indiv, group.by = "mut", order=c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known"),
              pt.size=1.7, cols=c("grey", "grey", "red", "red", "red")))
print(DimPlot(sample_seurat_indiv, group.by = "mut", order=c("mutation, >=3 SNV reads", "mutation, 2 SNV reads", "mutation, 1 SNV read", "no mutation", "not known"),
              pt.size=1.7, cols=c("grey", "grey", "red", "red", "red")) + NoLegend())
print(DimPlot(sample_seurat_indiv, group.by = "ref", order=c("ref nucleotide, >=3 SNV reads", "ref nucleotide, 2 SNV reads", "ref nucleotide, 1 SNV read", "no ref nucleotide", "not known"),
              pt.size=1.5, cols=c("grey", "grey20", "yellow", "orange", "red")))
print(DimPlot(sample_seurat_indiv, group.by = "ref", order=c("ref nucleotide, >=3 SNV reads", "ref nucleotide, 2 SNV reads", "ref nucleotide, 1 SNV read", "no ref nucleotide", "not known"),
              pt.size=1.5, cols=c("grey", "grey20", "yellow", "orange", "red")) + NoLegend())
ggplot(df[which(df$Var1!="not known"),], aes(x=Var2, y=value, fill=Var1))+
  geom_bar(stat="identity", position="fill")+
  theme_classic()+
  scale_fill_custom()+
  xlab("") + ylab("% of cells") + ggtitle("mutation detection in covered cells")
dev.off()
save(pval, file=file.path(screadcounts_dir, name_sample[1], "Analyses", paste0(filter_by_clonality, filter_by_SBS35, "mutations_", name_sample[1], "_pval_mut_vs_nomut.RData")))


#################################################################
##### 3. DIFFERENTIAL EXPRESSION ON CISPLATIN MUTATED CELLS #####
#################################################################
# Try to identify specific genes of cells with cisplatin mutation in the sample vs tumor cells that do not have the cisplatin mutation

# Choose which mutations to filter on
filter_by_clonality = c("clonal")
filter_by_SBS35 = c(0.75, 0.9)[c(1)] # 0 = not SBS35 related, 0.75 = SBS35 related with proba of 0.75, 0.9 = SBS35 related with proba of 0.9

##### 1. Concatenate mutation info per cell #####
#------------------------------------------------
recap_by_cell <- screadcounts_results %>% 
  filter(Clonality %in% filter_by_clonality) %>%
  filter(SBS35 >= filter_by_SBS35) %>%
  filter(PASS == "yes") %>% # keep only mutations that pass the criteria in WGS (to remove possible noise)
  group_by(ReadGroup) %>% 
  summarise(SNV_total_counts = sum(SNVCount), Ref_total_counts = sum(RefCount), VAF = SNV_total_counts / (SNV_total_counts + Ref_total_counts))

recap_by_cell <- recap_by_cell[which(recap_by_cell$SNV_total_counts!=0 | recap_by_cell$Ref_total_counts!=0 ),] # remove rows that have 0 SNV and 0 Ref counts = not covered (~1% of the rows)
recap_by_cell <- recap_by_cell[which(recap_by_cell$ReadGroup%in%colnames(sample_seurat_indiv)),] # screadcounts takes cells without QC -->  some of them are not in our postQC seurat object

# identify cells with/without mutation
cells_with_mut <- recap_by_cell %>% filter(SNV_total_counts > 0) %>% pull(ReadGroup)
cells_without_mut <- recap_by_cell %>% filter(SNV_total_counts == 0) %>% pull(ReadGroup)
sample_seurat_indiv$cells_cisplat_mut <- ifelse(rownames(sample_seurat_indiv@meta.data) %in% cells_with_mut, "cisplatin mutation", 
                                                ifelse(rownames(sample_seurat_indiv@meta.data) %in% cells_without_mut, "no cisplatin mutation", NA))


##### 2. Perform differential expression between cells with/without cisplatin mutation #####
#-------------------------------------------------------------------------------------------
DefaultAssay(sample_seurat_indiv) <- "RNA" # the assay to use for differential expression is RNA and not SCT (recommended by satijalab)
Idents(sample_seurat_indiv) <- "cells_cisplat_mut"

diff_exp_cisplat_vs_no_cisplat <- FindMarkers(sample_seurat_indiv, ident.1 = "cisplatin mutation", ident.2 = "no cisplatin mutation", logfc.threshold=0, min.pct = 0) # findmarkers can compute logFC even when gene not detected in a population because avg logFC = differences of the log(mean(expm1(x)+1))
diff_exp_cisplat_vs_no_cisplat_sig <- diff_exp_cisplat_vs_no_cisplat %>% filter(p_val_adj <= 0.1)





