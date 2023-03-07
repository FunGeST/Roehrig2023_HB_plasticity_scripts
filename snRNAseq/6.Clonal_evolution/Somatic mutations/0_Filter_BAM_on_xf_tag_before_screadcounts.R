


### EXPLANATIONS ###
#-------------------

# https://www.10xgenomics.com/resources/analysis-guides/tutorial-navigating-10x-barcoded-bam-files

# Answer of Karen Lai from 10X Genomics for xf=0:
#It means that the read is not confidently mapped to the transcriptome, does not map to exactly one feature, and is not treated as a UMI count. 
#For example, the read could be mapped uniquely to the genome (with MAPQ 255), but mapped to multiple or no annotations in the transcriptome. 
#The 17 means it confidently mapped to the genome and a feature (i.e gene), and 25 is tagging the read that was used for UMI counting. Some of the reads with xf 17 can be duplicates. In short, the collapsing of the UMIs is done using the xf tag. 



### RNA BAM 
#----------

# filter on xf=25 to keep reliable reads
samtools view -bS /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3610T/outs/gex_possorted_bam.bam -d xf:25 > /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3610T/outs/gex_possorted_bam_xf_25.bam

# create index
samtools index /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3610T/outs/gex_possorted_bam_xf_25.bam  /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3610T/outs/gex_possorted_bam_xf_25.bam.bai



