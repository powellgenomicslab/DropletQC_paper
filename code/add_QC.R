# Script information -----------------------------------------------------------

# title: Add QC metrics
# author: Walter Muskovic
# date: 2020-12-16
# description: In this script we aim to combine each Seurat object's metadata
# with additional QC metrics and other info 



# Import libraries -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(DropletQC)
  library(DropletUtils)
})



# Define directories & import Seurat object ------------------------------------

# Get the sample on which to gather QC metrics from the command line
# argument (integer 1-4)
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- c("GBM", "HL", "PBMC", "MB", "GBM5p")[i]

# Import Seurat
input.seurat <- readRDS(str_glue('../data/{sample_name}/seurat_SCT_anno.rds'))



# Add emptyDrops info ----------------------------------------------------------

# Import results from running emptyDrops
ed <- readRDS(str_glue('../data/{sample_name}/ED.rds')) %>% 
  as.data.frame() %>%
  rownames_to_column(var="CB") %>%
  drop_na() %>%
  filter(CB%in%colnames(input.seurat))
# should match exactly
all(ed$CB==colnames(input.seurat))

# Add info to the Seurat object metadata
input.seurat$emptyDrops_LogProb <- ed$LogProb
input.seurat$emptyDrops_PValue <- ed$PValue
input.seurat$emptyDrops_FDR <- ed$FDR



# Add cell barcodes and UMAP coords --------------------------------------------

# UMAP coords
input.seurat <- AddMetaData(input.seurat,
                            as.data.frame(input.seurat@reductions$umap@cell.embeddings))

# Cell barocdes
input.seurat$cell_barcode <- colnames(input.seurat)

# log10 UMI count
input.seurat$UMI_count <- input.seurat$nCount_RNA
input.seurat$log10_UMI_count <- log10(input.seurat$nCount_RNA)



# Add nuclear fraction calculated with dropletQC -------------------------------

nf <- nuclear_fraction_tags(outs = str_glue('../data/{sample_name}/outs'), cores = 16)
nf <-rename(nf, nuclear_fraction = "nuclear_fraction_dropletQC")
input.seurat <- AddMetaData(input.seurat, nf)



# Import nuclear fraction calculated with velocyto output ----------------------

# Import nuclear fraction and match order with Seurat object
nf <- read_csv(
  file = str_glue('../data/{sample_name}/velocyto/nuclear_fraction.csv'),
  col_names = c("CB", "nf"), skip = 1)
nf$CB <- str_split_fixed(nf$CB, ":", 2)[,2]
nf$CB <- substr(nf$CB,1,nchar(nf$CB)-1)
nf$CB <- str_glue('{nf$CB}-1')
nf <- nf[match(colnames(input.seurat), nf$CB),]
all(nf$CB==colnames(input.seurat))

# Add to the Seurat object
input.seurat$nuclear_fraction_velocyto <- nf$nf



# Save out ---------------------------------------------------------------------

saveRDS(input.seurat, str_glue('../data/{sample_name}/seurat_SCT_anno_QC.rds'))

print(str_glue('Finished adding QC metrics for sample: {sample_name}'))
