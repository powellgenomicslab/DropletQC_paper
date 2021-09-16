#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 8
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=2G
#$ -N get_MK_ref

# Prepare to download SRR
conda activate Seurat
Rscript prep_mouse_kidney.R

# Download and align fastqs
while read sample_name; do
  echo "Submitting script for ${sample_name}"
  qsub align_mouse_kidney.sh ${sample_name}
done < sample_names.tsv

# Run QC
Rscript mouse_kidney_data_QC.R


