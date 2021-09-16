#!/bin/bash

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#         File: pipeline.sh
#               Execute dropletQC analysis
#       Author: Walter Muskovic
#  Modified on: 2020/12/17
#      Version: 0.0.1
#      Example: ./pipeline.sh
#               Executing this script will create the analysis and figures 
#               presented in the manuscript Muskovic et al, journal_name, 2021
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

# Download 10x Genomics scRNA-seq data sets and required annotation files
qsub get_10x.sh

# Download and prepare mouse brain reference datasets to help with the 
# annotation of the mouse brain cell types
qsub get_mouse_brain_reference_data.sh

# Download and prepare mouse kidney scRNA-seq datsets from Denisenko et al
qsub get_mouse_kidney_data.sh

# Process 5' mouse splenocytes dataset
Rscript MS_QC.R

# Process HEK293 dataset
Rscript apo_QC.R

# Run DropletUtils::emptyDrops
qsub run_emptyDrops.sh

# Run velocyto for 10x samples and calculate the nuclear fraction QC metric 
# using the velocyto output
qsub -hold_jid download_10x -t 1-4:1 velocyto_10x.sh

# Perform a quick cell type annotation
qsub cell_type_annotation.sh

# Combine the nuclear fraction scores and emptyDrops info with the metadata of 
# each Seurat object
qsub add_QC.sh

# Assess the nuclear fraction QC metric for each sample
qsub assess_QC.sh

# Create figures
qsub figures.sh

# Create supplementary figures
qsub supplementary_figure_1.sh
bash supplementary_figure_2.sh
qsub supplementary_figure_3.sh
qsub supplementary_figure_4.sh
qsub supplementary_figure_5.sh
qsub supplementary_figure_6.sh
qsub supplementary_figure_7.sh
qsub supplementary_figure_8_excluded.sh

# Speed test
for run in {1..10}; do
  qsub -t 1-4:1 speed_test.sh;
done
