#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 8
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=16G
#$ -N get_MB_ref

conda activate Seurat
Rscript get_mouse_brain_reference_data.R
