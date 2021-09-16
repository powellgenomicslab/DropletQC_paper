#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 8
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=2G
#$ -N speed_test

conda activate Seurat
echo "Running a speed test for sample ${SGE_TASK_ID}"
Rscript speed_test.R $SGE_TASK_ID
