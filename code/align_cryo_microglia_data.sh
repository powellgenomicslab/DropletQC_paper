#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 32
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=10G
#$ -N STAR

# Create genome index
cd /directflow/SCCGGroupShare/projects/walmus/repositories/cryo_microglia
conda activate STAR
#STAR  --runMode genomeGenerate \
#  --runThreadN 32 \
#  --genomeDir data/STAR \
#  --genomeFastaFiles data/Macaca_mulatta.Mmul_10.dna.toplevel.fa \
#  --sjdbGTFfile data/Macaca_mulatta.Mmul_10.103.chr.gtf

# Align with STARsolo
ulimit -n -H 
ulimit -n -S
ulimit -n 4096

for sample in 85T_Fresh 86T_Fresh 85T_Cryo 86T_Cryo; do
echo "aligning ${sample}"

STAR --genomeDir data/STAR \
  --readFilesIn data/${sample}_3.fastq.gz data/${sample}_2.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17  --soloUMIlen 10 \
  --soloCellFilter None \
  --soloCBwhitelist data/whitelist.txt \
  --runThreadN 32 \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
  --outFileNamePrefix data/alignments/${sample} \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 60000000000

echo "done"
done

# Index
for sample in 85T_Fresh 86T_Fresh 85T_Cryo 86T_Cryo; do
echo "creating index for ${sample}"
samtools index data/alignments/${sample}Aligned.sortedByCoord.out.bam
echo "done"
done

