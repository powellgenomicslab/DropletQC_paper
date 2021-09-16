#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 32
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=4G
#$ -N velocyto_10x

## Generate loom file containing spliced and unspliced matrices using velocyto

conda activate velocyto

# Define sample and annotation to use
samples=("GBM" "PBMC" "HL" "MB")
orgs=("human" "human" "human" "mouse")
sample=${samples[$SGE_TASK_ID - 1]}
org=${orgs[$SGE_TASK_ID - 1]}

# Run velocyto
echo "#># Started running velocyto for sample ${sample}, using ${org} annotation file"
velocyto run10x --samtools-threads 32 -vv ../data/${sample} ../data/${org}.gtf
echo "#># Completed running velocyto"

echo "#># Started calculating nuclear fraction"
conda activate scvelo
python --version
#Python 3.6.10
python nuclear_fraction.py ../data/${sample}/velocyto $sample
echo "#># Completed calculating nuclear fraction"
