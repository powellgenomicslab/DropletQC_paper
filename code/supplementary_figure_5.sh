#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 4
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=4G
#$ -N supp_figure_2


echo "Started running scripts to create supplementary figure 5 for the manuscript"

qsub get_cryo_microglia_data.R
qsub -hold_jid get_data align_cryo_microglia_data.R
Rscript plot_cryo_microglia_data.R
