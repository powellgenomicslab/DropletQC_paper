# Script information -----------------------------------------------------------

# title: Run emptyDrops
# author: Walter Muskovic
# date: 2021-09-07
# description: This script will run emptyDrops from the DropletUtils package 



# Import libraries --------------------------------------------------------

library(tidyverse)
library(DropletUtils)
library(Seurat)
library(tictoc)
library(Matrix)



# Parse arguments ---------------------------------------------------------

i <- as.integer(commandArgs(trailingOnly=TRUE)[1])



# Run emptyDrops ----------------------------------------------------------

# Create function to run emptyDrops and save out the results

emptyDropsExtra <- function(sample_name, lower_threshold=i, ed_seed = 123){
  
  # Import UMI counts for *all* barcodes from raw_feature_bc_matrix
  print(str_glue('Importing raw_feature_bc_matrix for sample {sample_name}'))
  tenx.counts <- Read10X(
    str_glue('data/{sample_name}/outs/raw_feature_bc_matrix/')
  )
  
  ## Run emptryDrops()
  print(str_glue('Running empty drops for sample {sample_name} with lower parameter = {lower_threshold}'))
  tic()
  set.seed(ed_seed)
  e.out <- emptyDrops(tenx.counts, lower = lower_threshold)
  is.cell <- e.out$FDR <= 0.01
  print(sum(is.cell, na.rm = TRUE))
  print(table(Limited=e.out$Limited, Significant=is.cell))
  # Restrict to barcodes called as cells by emptyDrops
  tenx.counts.filtered <- tenx.counts[,which(is.cell)]
  toc()
  
  ## Save out emptyDrops results
  print(str_glue('Saving out cell barcodes for sample {sample_name}'))
  if(!dir.exists(str_glue('data/{sample_name}/comparison'))){
    dir.create(str_glue('data/{sample_name}/comparison'))
  }
  write_tsv(x = data.frame(sort(colnames(tenx.counts.filtered))),
            file = str_glue('data/{sample_name}/comparison/barcodes_{lower_threshold}.tsv'),
            col_names = FALSE) 
  
  print(str_glue('Completed running empty drops for sample {sample_name}'))
}

c("HL", "GBM", "PBMC", "MB") %>% walk(emptyDropsExtra)
