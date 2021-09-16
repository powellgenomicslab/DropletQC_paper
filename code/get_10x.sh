#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 16
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=1G
#$ -N download_10x


## GBM
mkdir -p ../data/GBM ../data/GBM/outs ../data/GBM/outs/filtered_feature_bc_matrix
echo "#># Downloading GBM data files"

# possorted_genome_bam
wget -nc --output-document="../data/GBM/outs/possorted_genome_bam.bam" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_SC3v3_Human_Glioblastoma/Parent_SC3v3_Human_Glioblastoma_possorted_genome_bam.bam
echo "793bad85dbbdcfee12319b3ae67b8ce9  ../data/GBM/outs/possorted_genome_bam.bam" | md5sum -c -
wget -nc --output-document="../data/GBM/outs/possorted_genome_bam.bam.bai" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_SC3v3_Human_Glioblastoma/Parent_SC3v3_Human_Glioblastoma_possorted_genome_bam.bam.bai

# raw_feature_bc_matrix
wget -nc --output-document="../data/GBM/outs/raw_feature_bc_matrix.tar.gz" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_SC3v3_Human_Glioblastoma/Parent_SC3v3_Human_Glioblastoma_raw_feature_bc_matrix.tar.gz
tar -zxvf ../data/GBM/outs/raw_feature_bc_matrix.tar.gz -C ../data/GBM/outs/
rm ../data/GBM/outs/raw_feature_bc_matrix.tar.gz

# raw_feature_bc_matrix.h5
wget -nc --output-document="../data/GBM/outs/raw_feature_bc_matrix.h5" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_SC3v3_Human_Glioblastoma/Parent_SC3v3_Human_Glioblastoma_raw_feature_bc_matrix.h5

# fastq
wget -nc --output-document="../data/GBM/fastqs.tar" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_SC3v3_Human_Glioblastoma/Parent_SC3v3_Human_Glioblastoma_fastqs.tar



## PBMC
mkdir -p ../data/PBMC ../data/PBMC/outs ../data/PBMC/outs/filtered_feature_bc_matrix
echo "#># Downloading PBMC data files"

# possorted_genome_bam
wget -nc --output-document="../data/PBMC/outs/possorted_genome_bam.bam" https://cg.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_possorted_genome_bam.bam
echo "8c36cff42c9d7a73d0e8bfc4c6339787  ../data/PBMC/outs/possorted_genome_bam.bam" | md5sum -c -
wget -nc --output-document="../data/PBMC/outs/possorted_genome_bam.bam.bai" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_possorted_genome_bam.bam.bai

# raw_feature_bc_matrix
wget -nc --output-document="../data/PBMC/outs/raw_feature_bc_matrix.tar.gz" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_raw_feature_bc_matrix.tar.gz
tar -zxvf ../data/PBMC/outs/raw_feature_bc_matrix.tar.gz -C ../data/PBMC/outs/
rm ../data/PBMC/outs/raw_feature_bc_matrix.tar.gz

# raw_feature_bc_matrix.h5
wget -nc --output-document="../data/PBMC/outs/raw_feature_bc_matrix.h5" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_raw_feature_bc_matrix.h5

# fastq
wget -nc --output-document="../data/PBMC/fastqs.tar" https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_fastqs.tar



## Hodgkin's Lymphoma
mkdir -p ../data/HL ../data/HL/outs ../data/HL/outs/filtered_feature_bc_matrix
echo "#># Downloading Hodgkin's Lymphoma data files"

# possorted_genome_bam
wget -nc --output-document="../data/HL/outs/possorted_genome_bam.bam" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_HodgkinsLymphoma/Parent_NGSC3_DI_HodgkinsLymphoma_possorted_genome_bam.bam
echo "75672e0dd546995c74fd2a0780f39c7b  ../data/HL/outs/possorted_genome_bam.bam" | md5sum -c -
wget -nc --output-document="../data/HL/outs/possorted_genome_bam.bam.bai" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_HodgkinsLymphoma/Parent_NGSC3_DI_HodgkinsLymphoma_possorted_genome_bam.bam.bai

# raw_feature_bc_matrix
wget -nc --output-document="../data/HL/outs/raw_feature_bc_matrix.tar.gz" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_HodgkinsLymphoma/Parent_NGSC3_DI_HodgkinsLymphoma_raw_feature_bc_matrix.tar.gz
tar -zxvf ../data/HL/outs/raw_feature_bc_matrix.tar.gz -C ../data/HL/outs/
rm ../data/HL/outs/raw_feature_bc_matrix.tar.gz

# raw_feature_bc_matrix.h5
wget -nc --output-document="../data/HL/outs/raw_feature_bc_matrix.h5" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_HodgkinsLymphoma/Parent_NGSC3_DI_HodgkinsLymphoma_raw_feature_bc_matrix.h5

# fastq
wget -nc --output-document="../data/HL/fastqs.tar" https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/4.0.0/Parent_NGSC3_DI_HodgkinsLymphoma/Parent_NGSC3_DI_HodgkinsLymphoma_fastqs.tar



## Mouse brain
mkdir -p ../data/MB ../data/MB/outs ../data/MB/outs/filtered_feature_bc_matrix
echo "#># Downloading mouse brain data files"

# possorted_genome_bam
wget -nc --output-document="../data/MB/outs/possorted_genome_bam.bam" https://cg.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_DI_Neuron_10K/SC3_v3_NextGem_DI_Neuron_10K_possorted_genome_bam.bam
echo "567d76ec472f96785cc0fbab16b78c3e  ../data/MB/outs/possorted_genome_bam.bam" | md5sum -c -
wget -nc --output-document="../data/MB/outs/possorted_genome_bam.bam.bai" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_DI_Neuron_10K/SC3_v3_NextGem_DI_Neuron_10K_possorted_genome_bam.bam.bai

# raw_feature_bc_matrix
wget -nc --output-document="../data/MB/outs/raw_feature_bc_matrix.tar.gz" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_DI_Neuron_10K/SC3_v3_NextGem_DI_Neuron_10K_raw_feature_bc_matrix.tar.gz
tar -zxvf ../data/MB/outs/raw_feature_bc_matrix.tar.gz -C ../data/MB/outs/
rm ../data/MB/outs/raw_feature_bc_matrix.tar.gz

# raw_feature_bc_matrix.h5
wget -nc --output-document="../data/MB/outs/raw_feature_bc_matrix.h5" https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_DI_Neuron_10K/SC3_v3_NextGem_DI_Neuron_10K_raw_feature_bc_matrix.h5

# fastq
wget -nc --output-document="../data/MB/fastqs.tar" https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/4.0.0/SC3_v3_NextGem_DI_Neuron_10K/SC3_v3_NextGem_DI_Neuron_10K_fastqs.tar



## GBM 5'
mkdir -p ../data/GBM5p ../data/GBM5p/outs ../data/GBM5p/outs/filtered_feature_bc_matrix
echo "#># Downloading GBM data files"

# possorted_genome_bam
wget -nc --output-document="../data/GBM5p/outs/possorted_genome_bam.bam" https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/4.0.0/Parent_SC5v1_Human_Glioblastoma/Parent_SC5v1_Human_Glioblastoma_possorted_genome_bam.bam
echo "d5e712edca1590ba2e80ed0e5d62b239  ../data/GBM5p/outs/possorted_genome_bam.bam" | md5sum -c -
wget -nc --output-document="../data/GBM5p/outs/possorted_genome_bam.bam.bai" https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/Parent_SC5v1_Human_Glioblastoma/Parent_SC5v1_Human_Glioblastoma_possorted_genome_bam.bam.bai

# raw_feature_bc_matrix
wget -nc --output-document="../data/GBM5p/outs/raw_feature_bc_matrix.tar.gz" https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/Parent_SC5v1_Human_Glioblastoma/Parent_SC5v1_Human_Glioblastoma_raw_feature_bc_matrix.tar.gz
tar -zxvf ../data/GBM5p/outs/raw_feature_bc_matrix.tar.gz -C ../data/GBM5p/outs/
rm ../data/GBM5p/outs/raw_feature_bc_matrix.tar.gz



## Mouse splenocytes 5' (MS)
mkdir -p ../data/MS ../data/MS/outs ../data/MS/outs/filtered_feature_bc_matrix
echo "#># Downloading MS data files"

# possorted_genome_bam
wget -nc --output-document="../data/MS/outs/possorted_genome_bam.bam" https://cf.10xgenomics.com/samples/cell-vdj/6.0.1/SC5v2_mouseSplenocytes_10Kcells_Connect_single_channel_SC5v2_mouseSplenocytes_10Kcells_Connect_single_channel/SC5v2_mouseSplenocytes_10Kcells_Connect_single_channel_SC5v2_mouseSplenocytes_10Kcells_Connect_single_channel_count_sample_alignments.bam
echo "53a2f2951a4c6e56f57f296e48d27300  ../data/MS/outs/possorted_genome_bam.bam" | md5sum -c -
wget -nc --output-document="../data/MS/outs/possorted_genome_bam.bam.bai" https://cf.10xgenomics.com/samples/cell-vdj/6.0.1/SC5v2_mouseSplenocytes_10Kcells_Connect_single_channel_SC5v2_mouseSplenocytes_10Kcells_Connect_single_channel/SC5v2_mouseSplenocytes_10Kcells_Connect_single_channel_SC5v2_mouseSplenocytes_10Kcells_Connect_single_channel_count_sample_alignments.bam.bai

# raw_feature_bc_matrix
wget -nc --output-document="../data/MS/outs/sample_feature_bc_matrix.tar.gz" https://cf.10xgenomics.com/samples/cell-vdj/6.0.1/SC5v2_mouseSplenocytes_10Kcells_Connect_single_channel_SC5v2_mouseSplenocytes_10Kcells_Connect_single_channel/SC5v2_mouseSplenocytes_10Kcells_Connect_single_channel_SC5v2_mouseSplenocytes_10Kcells_Connect_single_channel_count_sample_feature_bc_matrix.tar.gz
tar -zxvf ../data/MS/outs/sample_feature_bc_matrix.tar.gz -C ../data/MS/outs/
rm ../data/MS/outs/sample_feature_bc_matrix.tar.gz



## Reference files
echo "#># Downloading reference annotation files"

# refdata-gex-GRCh38-2020-A 
wget -nc --output-document="../data/refdata-gex-GRCh38-2020-A.tar.gz" https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
echo "dfd654de39bff23917471e7fcc7a00cd  ../data/refdata-gex-GRCh38-2020-A.tar.gz" | md5sum -c -
tar -zxvf ../data/refdata-gex-GRCh38-2020-A.tar.gz -C ../data/
rm ../data/refdata-gex-GRCh38-2020-A.tar.gz
cp ../data/refdata-gex-GRCh38-2020-A/genes/genes.gtf ../data/human.gtf
rm -rf ../data/refdata-gex-GRCh38-2020-A

# refdata-gex-mm10-2020-A
wget -nc --output-document="../data/refdata-gex-mm10-2020-A.tar.gz" https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
echo "886eeddde8731ffb58552d0bb81f533d  ../data/refdata-gex-mm10-2020-A.tar.gz" | md5sum -c -
tar -zxvf ../data/refdata-gex-mm10-2020-A.tar.gz -C ../data/
rm ../data/refdata-gex-mm10-2020-A.tar.gz
cp ../data/refdata-gex-mm10-2020-A/genes/genes.gtf ../data/mouse.gtf
rm -rf ../data/refdata-gex-mm10-2020-A
