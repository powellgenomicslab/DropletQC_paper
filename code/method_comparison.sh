#!/bin/bash

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#         File: method_comparison.sh
#        Title: Comparison of empty droplet detection methods
#       Author: Walter Muskovic
#  Modified on: 2021/09/07
#      Version: 0.0.1
#      Example: ./method_comparison.sh
#  Description: Run this script to compare the ability of different methods to 
#               exclude empty droplets from the four 10x Genomics datasets; GBM,
#               MB, HL and PBMCs. Methods tested are; DropletUtils::emptyDrops,
#               EmptyNN, CellBender, Cell Ranger 2.2.0 and Cell Ranger 6.1.1.
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



# EmptyDrops -------------------------------------------------------------------

# Run emptyDrops with 'lower' parameter set to 100 (default value)
conda activate Seurat
Rscript comparison_emptyDrops.R 100

# Run emptyDrops with 'lower' parameter set to 500
Rscript comparison_emptyDrops.R 500



# EmptyNN ----------------------------------------------------------------------

# Run emptynn with 'threshold' parameter set to 100 (default value)
Rscript comparison_emptyNN.R 100

# Run emptynn with 'threshold' parameter set to 500
Rscript comparison_emptyNN.R 500



# CellBender -------------------------------------------------------------------

qsub cellbender.sh



# Cell Ranger 6.1.1 ------------------------------------------------------------

qsub cellranger6.sh



# Combine output ---------------------------------------------------------------

Rscript combine_comparison.R


