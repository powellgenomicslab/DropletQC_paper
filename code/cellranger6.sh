#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 16
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=4G
#$ -N cellranger6

# Define paths
base_dir="/directflow/SCCGGroupShare/projects/walmus/repositories/dropletQC_paper/data"
human_transcriptome="/directflow/SCCGGroupShare/projects/walmus/share/cellranger/refdata-gex-GRCh38-2020-A"
mouse_transcriptome="/directflow/SCCGGroupShare/projects/walmus/share/cellranger/refdata-gex-mm10-2020-A"
cellranger_path="/directflow/SCCGGroupShare/projects/walmus/share/cellranger/cellranger-6.1.1"

# Align Hodgkin's lymphoma fastqs
mkdir $base_dir/HL/cellranger6
cd $base_dir/HL/cellranger6
${cellranger_path}/cellranger count --id=HL \
                   --transcriptome=${human_transcriptome} \
                   --fastqs=${base_dir}/HL/Parent_NGSC3_DI_HodgkinsLymphoma_fastqs \
                   --localcores=16 \
                   --localmem=64 \
                   --expect-cells=5000

# Align glioblastoma fastqs
mkdir $base_dir/GBM/cellranger6
cd $base_dir/GBM/cellranger6
${cellranger_path}/cellranger count --id=GBM \
                   --transcriptome=${human_transcriptome} \
                   --fastqs=${base_dir}/GBM/Parent_SC3v3_Human_Glioblastoma_fastqs/Parent_SC3v3_Human_Glioblastoma \
                   --localcores=16 \
                   --localmem=64 \
                   --expect-cells=5000

# Align mouse brain fastqs
mkdir $base_dir/MB/cellranger6
cd $base_dir/MB/cellranger6
${cellranger_path}/cellranger count --id=MB \
                   --transcriptome=${mouse_transcriptome} \
                   --fastqs=${base_dir}/MB/SC3_v3_NextGem_DI_Neuron_10K_fastqs \
                   --localcores=16 \
                   --localmem=64 \
                   --expect-cells=10000

# Align PBMC fastqs
mkdir $base_dir/PBMC/cellranger6
cd $base_dir/PBMC/cellranger6
${cellranger_path}/cellranger count --id=PBMC \
                   --transcriptome=${human_transcriptome} \
                   --fastqs=${base_dir}/PBMC/Parent_NGSC3_DI_PBMC_fastqs \
                   --localcores=16 \
                   --localmem=64 \
                   --expect-cells=10000

