#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 4
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=4G
#$ -N get_data

# Activate sra-tools
conda activate sra-tools

# Prefetch SRA files
while read SRR; do
  echo "Processing $SRR"
  prefetch -o data/${SRR} $SRR &
done <data/SRR_Acc_List.txt
wait

# Extract fastq files using fastq-dump.2.10.0
while read SRR; do
  fastq-dump.2.10.0 --gzip --split-files -O data/ data/${SRR} &
done <data_track/SRR_Acc_List.txt
wait

# Download gene annotation and macaque genome (Mmul_30)
cd data
wget http://ftp.ensembl.org/pub/release-103/gtf/macaca_mulatta/Macaca_mulatta.Mmul_30.103.chr.gtf.gz &
wget http://ftp.ensembl.org/pub/release-103/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_30.dna.toplevel.fa.gz &

# Combine fastq files

# 85T_Fresh
cat SRR13198365_2.fastq.gz SRR13198366_2.fastq.gz SRR13198367_2.fastq.gz SRR13198368_2.fastq.gz SRR13198369_2.fastq.gz SRR13198370_2.fastq.gz SRR13198371_2.fastq.gz SRR13198372_2.fastq.gz > 85T_Fresh_2.fastq.gz &
cat SRR13198365_3.fastq.gz SRR13198366_3.fastq.gz SRR13198367_3.fastq.gz SRR13198368_3.fastq.gz SRR13198369_3.fastq.gz SRR13198370_3.fastq.gz SRR13198371_3.fastq.gz SRR13198372_3.fastq.gz > 85T_Fresh_3.fastq.gz &

# 86T_Fresh
cat SRR13198373_2.fastq.gz SRR13198374_2.fastq.gz SRR13198375_2.fastq.gz SRR13198376_2.fastq.gz SRR13198377_2.fastq.gz SRR13198378_2.fastq.gz SRR13198379_2.fastq.gz SRR13198380_2.fastq.gz > 86T_Fresh_2.fastq.gz &
cat SRR13198373_3.fastq.gz SRR13198374_3.fastq.gz SRR13198375_3.fastq.gz SRR13198376_3.fastq.gz SRR13198377_3.fastq.gz SRR13198378_3.fastq.gz SRR13198379_3.fastq.gz SRR13198380_3.fastq.gz > 86T_Fresh_3.fastq.gz &

# 85T_Cryo
cat SRR13198381_2.fastq.gz SRR13198382_2.fastq.gz SRR13198383_2.fastq.gz SRR13198384_2.fastq.gz > 85T_Cryo_2.fastq.gz &
cat SRR13198381_3.fastq.gz SRR13198382_3.fastq.gz SRR13198383_3.fastq.gz SRR13198384_3.fastq.gz > 85T_Cryo_3.fastq.gz &

# 86T_Cryo
cat SRR13198385_2.fastq.gz SRR13198386_2.fastq.gz SRR13198387_2.fastq.gz SRR13198388_2.fastq.gz > 86T_Cryo_2.fastq.gz &
cat SRR13198385_3.fastq.gz SRR13198386_3.fastq.gz SRR13198387_3.fastq.gz SRR13198388_3.fastq.gz > 86T_Cryo_3.fastq.gz &
wait

# Get count data
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE162663&format=file
