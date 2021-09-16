#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 4
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=16G
#$ -N align_MK


# Define variables and directories
sample_name=$1
transcriptome_dir="/directflow/SCCGGroupShare/projects/walmus/share/cellranger/refdata-gex-mm10-2020-A"
cd /directflow/SCCGGroupShare/projects/walmus/repositories/dropletQC_paper/data/MK/${sample_name}



# Prefetch SRA files
conda activate sra-tools
while read SRR; do
  echo "Prefetching $SRR"
  prefetch -o ${SRR} $SRR &
done < SRR_list.txt
wait



# Extract fastq files using fastq-dump.2.10.0
while read SRR; do
  echo "Extracting fastq files from $SRR"
  fastq-dump.2.10.0 --gzip --split-files -O ./ ./${SRR} &
done < SRR_list.txt
wait



# Combine fastq files in a separate directory
mkdir fastq
cat *_1.fastq.gz > fastq/${sample_name}_S1_L001_R1_001.fastq.gz
cat *_2.fastq.gz > fastq/${sample_name}_S1_L001_R2_001.fastq.gz



# Align and quantify with cellranger
cellranger count --id=${sample_name} \
                   --transcriptome=${transcriptome_dir} \
                   --fastqs=fastq \
                   --localcores=4 \
                   --localmem=64



# Remove unneeded SRR files
rm SRR*
# Remove unneeded fastq files
rm -rf fastq

# Copy BAM and raw counts dir, delete remaining files
#cp ${sample_name}/outs/XXX ./XXX
#cp ${sample_name}/outs/XXX ./XXX
#cp ${sample_name}/outs/XXX ./XXX
#rm -rf ${sample_name}
