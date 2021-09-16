#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 4
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=4G
#$ -N supp_figures

conda activate Seurat
echo "Started running script to create supplementary figure 5 for the manuscript"
Rscript supplementary_figure_5.R 