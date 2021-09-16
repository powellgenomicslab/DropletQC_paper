#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 16
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=4G
#$ -N cellbender

source activate cellbender
base_dir="/directflow/SCCGGroupShare/projects/walmus/repositories/dropletQC_paper/data"

# Hodgkin's lymphoma
cellbender remove-background \
                 --input ${base_dir}/HL/outs/raw_feature_bc_matrix.h5 \
                 --output ${base_dir}/HL/comparison/cb.h5 \
                 --expected-cells 5000 \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150

# Glioblastoma
cellbender remove-background \
                 --input ${base_dir}/GBM/outs/raw_feature_bc_matrix.h5 \
                 --output ${base_dir}/GBM/comparison/cb.h5 \
                 --expected-cells 5000 \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150

# Mouse brain
cellbender remove-background \
                 --input ${base_dir}/MB/outs/raw_feature_bc_matrix.h5 \
                 --output ${base_dir}/MB/comparison/cb.h5 \
                 --expected-cells 10000 \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150             

# PBMC
cellbender remove-background \
                 --input ${base_dir}/PBMC/outs/raw_feature_bc_matrix.h5 \
                 --output ${base_dir}/PBMC/comparison/cb.h5 \
                 --expected-cells 10000 \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150            
