#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 8
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=4G
#$ -N emptyDrops

conda activate Seurat
echo "Started running emptyDrops"
Rscript run_emptyDrops.R
echo "Finished running emptyDrops"
