#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 4
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=4G
#$ -N addQC

conda activate Seurat
echo "Started running script to add QC metrics for sample $SGE_TASK_ID"
Rscript add_QC.R $SGE_TASK_ID
