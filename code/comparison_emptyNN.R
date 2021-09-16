# Script information -----------------------------------------------------------

# title: Run emptyDrops
# author: Walter Muskovic
# date: 2021-09-07
# description: This script will run emptyDrops from the DropletUtils package 



# Import libraries --------------------------------------------------------

library(tidyverse)
library(EmptyNN)
library(Seurat)
library(tictoc)
library(Matrix)



# Parse arguments ---------------------------------------------------------

i <- as.integer(commandArgs(trailingOnly=TRUE)[1])



# Run emptyDrops ----------------------------------------------------------

# Create function to run EmptyNN and save out the results

EmptyNNExtra <- function(sample_name, threshold_value = i){
  
  # Import UMI counts for *all* barcodes from raw_feature_bc_matrix and
  # transpose so that rows are barcodes and columns are genes
  print(str_glue('Importing raw_feature_bc_matrix for sample {sample_name}'))
  tenx.counts <- Read10X(
    str_glue('data/{sample_name}/outs/raw_feature_bc_matrix/')
  ) %>%
    t()
  
  # Run emptynn()
  tic()
  nn.res <- emptynn(tenx.counts,
                    threshold = threshold_value,
                    k = 10,
                    iteration = 10,
                    verbose = TRUE)
  toc()
  
  ## Save out emptynn results
  print(str_glue('Saving out emptynn results for sample {sample_name}'))
  saveRDS(nn.res, str_glue('data/{sample_name}/comparison/emptynn_{threshold_value}.rds'))
  
  print(str_glue('Completed running emptynn for sample {sample_name}'))
}

c("HL", "GBM", "PBMC", "MB") %>% walk(EmptyNNExtra)
