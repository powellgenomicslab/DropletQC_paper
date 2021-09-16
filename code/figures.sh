#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 4
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=4G
#$ -N plot_figures

conda activate Seurat
echo "Started running script to create figures for the manuscript"
Rscript figures.R 
