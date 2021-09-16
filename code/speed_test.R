# Script information -----------------------------------------------------------

# title: Speed test
# author: Walter Muskovic
# date: 2021-01-27
# description: In this script we will run dropletQC::nuclear_fraction() and see 
# how long it takes to run through the BAM file



# Import libraries -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(dropletQC)
})



# Get sample name --------------------------------------------------------------

# Get the sample on which to run the speed test from the command line argument 
# (integer 1-4)
cmd_args <- commandArgs(trailingOnly=TRUE)
i <- as.integer(cmd_args[1])
sample_name <- c("GBM", "HL", "PBMC", "MB")[i]



# Number of reads --------------------------------------------------------------

# Get the number of reads in the BAM file
bam_info <- Rsamtools::idxstatsBam(
  str_glue('/directflow/SCCGGroupShare/projects/walmus/repositories/dropletQC_paper/data/{sample_name}/outs/possorted_genome_bam.bam'))
num_reads <- sum(bam_info$mapped)



# Run speed tests --------------------------------------------------------------

start_time <- Sys.time()
nf <- nuclear_fraction_tags(outs = str_glue('/directflow/SCCGGroupShare/projects/walmus/repositories/dropletQC_paper/data/{sample_name}/outs'),
                            cores = 8)
end_time <- Sys.time()
time_elapsed1 <- as.numeric(difftime(time1 = end_time,
                                    time2 = start_time,
                                    units = "sec"))


start_time <- Sys.time()
nf <- nuclear_fraction_annotation(annotation_path = ifelse(test = sample_name!="MB",
                                                           yes = "/directflow/SCCGGroupShare/projects/walmus/repositories/dropletQC_paper/data/human.gtf",
                                                           no = "/directflow/SCCGGroupShare/projects/walmus/repositories/dropletQC_paper/data/mouse.gtf"),
                                  bam = str_glue('/directflow/SCCGGroupShare/projects/walmus/repositories/dropletQC_paper/data/{sample_name}/outs/possorted_genome_bam.bam'),
                                  barcodes = str_glue('/directflow/SCCGGroupShare/projects/walmus/repositories/dropletQC_paper/data/{sample_name}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'),
                                  cores = 8)
end_time <- Sys.time()
time_elapsed2 <- as.numeric(difftime(time1 = end_time,
                                     time2 = start_time,
                                     units = "sec"))



# Save time --------------------------------------------------------------------

# Create data frame to hold results
speed_results <- data.frame(sample=sample_name,
                            reads=num_reads,
                            nuclear_fraction_tags=time_elapsed1,
                            nuclear_fraction_annotation=time_elapsed2)

# Save out
if(!file.exists("/directflow/SCCGGroupShare/projects/walmus/repositories/dropletQC_paper/data_track/run_times.tsv")){
  write_tsv(speed_results, "/directflow/SCCGGroupShare/projects/walmus/repositories/dropletQC_paper/data_track/run_times.tsv", 
            col_names = TRUE, append = FALSE)
} else {
  write_tsv(speed_results, "/directflow/SCCGGroupShare/projects/walmus/repositories/dropletQC_paper/data_track/run_times.tsv", 
            col_names = FALSE, append = TRUE)
}

print(str_glue('Finished running speed test for sample {sample_name}'))
