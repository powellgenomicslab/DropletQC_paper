#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 32
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=6G
#$ -N annotation

# Activate conda environment that contains the Seurat R package
conda activate Seurat
echo "Started running cell type classification script"
Rscript cell_type_annotation.R
echo "Completed running cell type classification script"
