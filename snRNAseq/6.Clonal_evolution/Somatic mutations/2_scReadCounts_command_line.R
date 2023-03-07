

### HERE WE TAKE AS INPUT THE WHOLE BAM FILE = ALL CELLS, NOT ONLY TUMOR CELLS (WE WILL RESTRICT LATER IN ANALYSES)

# screadcounts usually runs from 1min to 15 min

##### 1. Look for tumor mutations in the corresponding tumor #####
#-----------------------------------------------------------------

### CHC2959T
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC2959T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2959T/scReadCounts_output.tsv

### CHC2960T
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2960T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC2960T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2960T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2960T/scReadCounts_output.tsv

### CHC3133T
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3133T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC3133T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3133T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3133T/scReadCounts_output.tsv

### CHC3377T
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC3377T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3377T/scReadCounts_output.tsv

### CHC3610T
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3610T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC3610T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3610T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3610T/scReadCounts_output.tsv

### CHC3662T
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3662T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC3662T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3662T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3662T/scReadCounts_output.tsv


##### 2. Look for tumor mutations in the corresponding non-tumor samples #####
#-----------------------------------------------------------------------------
# = use non-tumor samples as negative controls

### CHC2959N
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC2959T_CHC2960T_union.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2959N/scReadCounts_output_2959T_2960T.tsv

### CHC3377N
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377N/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC3377T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377N/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3377N/scReadCounts_output.tsv


##### 3. CHC2959T/60T: look for trunk and specific mutations in each sample #####
#--------------------------------------------------------------------------------

### CHC2959T
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC2959T_CHC2960T_union.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2959T/scReadCounts_output_2959T_2960T.tsv

### CHC2960T
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2960T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC2959T_CHC2960T_union.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2960T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2960T/scReadCounts_output_2959T_2960T.tsv


##### 4. Look for beta catenin mutations #####
#---------------------------------------------

### CHC2959T
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_beta_cat.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2959T/scReadCounts_output_beta_cat.tsv

### CHC2960T
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2960T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_beta_cat.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2960T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2960T/scReadCounts_output_beta_cat.tsv

### CHC3133T
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3133T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_beta_cat.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3133T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3133T/scReadCounts_output_beta_cat.tsv

### CHC3377T
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_beta_cat.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3377T/scReadCounts_output_beta_cat.tsv

### CHC3610T
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3610T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_beta_cat.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3610T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3610T/scReadCounts_output_beta_cat.tsv

### CHC3662T
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3662T/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_beta_cat.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3662T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3662T/scReadCounts_output_beta_cat.tsv

### CHC2959N
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_beta_cat.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2959N/scReadCounts_output_beta_cat.tsv

### CHC3377N
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377N/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_beta_cat.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377N/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3377N/scReadCounts_output_beta_cat.tsv



##### 5. Search 3133T mutations in the 2 NT 2959N and 3377N #####
#----------------------------------------------------------------

### CHC2959N
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC3133T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2959N/scReadCounts_output_3133T_mut.tsv

### CHC3377N
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377N/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC3133T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377N/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3377N/scReadCounts_output_3133T_mut.tsv


##### 5. Search 3662T mutations in the 2 NT 2959N and 3377N #####
#----------------------------------------------------------------

### CHC2959N
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC3662T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2959N/scReadCounts_output_3662T_mut.tsv

### CHC3377N
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377N/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC3662T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377N/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3377N/scReadCounts_output_3662T_mut.tsv


##### 6. Search 3377T mutations in the 2 NT 2959N and 3377N #####
#----------------------------------------------------------------

### CHC2959N
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC3377T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2959N/scReadCounts_output_3377T_mut.tsv

### CHC3377N
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377N/outs/gex_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC3377T.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC3377N/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC3377N/scReadCounts_output_3377T_mut.tsv



##### 4. TEST WITH ATAC

### CHC2959N
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/atac_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC2959T_CHC2960T_union.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959N/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2959N/ATAC_test_scReadCounts_output_2959T_2960T.tsv

### CHC2959T
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959T/outs/atac_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC2959T_CHC2960T_union.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2959T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2959T/ATAC_test_scReadCounts_output_2959T_2960T.tsv

### CHC2960T
source /home/amelie/miniconda2/etc/profile.d/conda.sh
conda activate python3.6
cd /home/amelie/SCReadCounts-1.1.8.Python-3.7
export PYTHON3=/home/amelie/miniconda2/envs/python3.6/bin/python3
bin/scReadCounts -r /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2960T/outs/atac_possorted_bam.bam -s /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/0.Input_data/mutation_table_CHC2959T_CHC2960T_union.tsv -G STARsolo -b /mnt/MD1200/NGS_works/Multiome_scRNAseq_scATACseq_works/HB_PJ2105182/Cellranger_outputs/Filtered_gtf_files_5/CHC2960T/outs/barcodes.tsv -o /mnt/MD1200_5/Workdir_amelie/snRNAseq_ATACseq_multiome/2.Analyses/1.snRNAseq/9.Clonal_evolution/Somatic_mutations/CHC2960T/ATAC_test_scReadCounts_output_2959T_2960T.tsv



