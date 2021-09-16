# Script information -------------------------------------------------------

# title: Combine output from different methods
# author: Walter Muskovic
# date: 2021-09-08
# description: In this script we are combining the output produced by; 
# emptyDrops, EmptyNN, CellBender and Cell Ranger 6.1.1.



# Import libraries  -------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(DropletQC)
})

samples <- c("HL", "GBM", "PBMC", "MB")
umi <- c(100,500)



# EmptyDrops --------------------------------------------------------------

process_ED <- function(sample_name, umi_cutoff){
  
  print(str_glue('Processing emptyDrops output for sample {sample_name} with umi cutoff {umi_cutoff}'))
  
  # Import empty drops results
  ed <- read_tsv(str_glue('data/{sample_name}/comparison/barcodes_{umi_cutoff}.tsv'), col_names = FALSE)
  
  # Import raw counts
  seurat.object <- Read10X(str_glue('data/{sample_name}/outs/raw_feature_bc_matrix'))
  
  # Subset to cells called by empty drops
  seurat.object <- seurat.object[,ed$X1]
  seurat.object <- CreateSeuratObject(seurat.object, project = sample_name)
  
  # Add mitochondrial gene content
  if(sample_name=="MB"){
    seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, "^mt-")
  } else {
    seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, "^MT-")
  }
  
  # Add nuclear fraction
  seurat.object$nuclear_fraction <- nuclear_fraction_tags(bam = str_glue('data/{sample_name}/outs/possorted_genome_bam.bam'),
                                                          barcodes = colnames(seurat.object),
                                                          cores=16)$nuclear_fraction
  
  # Identify any empty droplets
  ed <- identify_empty_drops(select(seurat.object@meta.data, nuclear_fraction, nCount_RNA), umi_rescue = 0) %>%
    mutate(percent.mt = seurat.object$percent.mt)
  
  # Add other relevant info
  ed <- mutate(ed,
               method="EmptyDrops",
               sample = sample_name,
               param = umi_cutoff)
  
  return(ed)
  }

# Save out
c(lapply(samples, function(i) process_ED(i, umi_cutoff=umi[1])),
lapply(samples, function(i) process_ED(i, umi_cutoff=umi[2]))) %>%
  do.call(rbind,.) %>%
  write_tsv("data/ed.tsv")



# EmptyNN -----------------------------------------------------------------

process_ENN <- function(sample_name, umi_cutoff){
  print(str_glue('Processing EmptyNN output for sample {sample_name} with umi cutoff {umi_cutoff}'))
  
  # Import EmptyNN results
  ed <- readRDS(str_glue('data/{sample_name}/comparison/emptynn_{umi_cutoff}.rds'))
  
  # Import raw counts
  seurat.object <- Read10X(str_glue('data/{sample_name}/outs/raw_feature_bc_matrix'))
  
  # Check that all cells were not excluded
  if(!any(ed$nn.keep)){ return(NULL) }
  
  # Subset to cells called by emptyNN
  seurat.object <- seurat.object[,ed$nn.keep]
  seurat.object <- CreateSeuratObject(seurat.object, project = sample_name)
  
  # Add mitochondrial gene content
  if(sample_name=="MB"){
    seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, "^mt-")
  } else {
    seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, "^MT-")
  }
  
  # Add nuclear fraction
  seurat.object$nuclear_fraction <- nuclear_fraction_tags(bam = str_glue('data/{sample_name}/outs/possorted_genome_bam.bam'),
                                                          barcodes = colnames(seurat.object),
                                                          cores=16)$nuclear_fraction
  
  # Identify any empty droplets
  ed <- identify_empty_drops(select(seurat.object@meta.data, nuclear_fraction, nCount_RNA), umi_rescue = 0) %>%
    mutate(percent.mt = seurat.object$percent.mt)
  
  # Add other relevant info
  ed <- mutate(ed,
               method="EmptyNN",
               sample = sample_name,
               param = umi_cutoff)
  return(ed)
}

# Save out
c(lapply(samples, function(i) process_ENN(i, umi_cutoff=umi[1])),
  lapply(samples, function(i) process_ENN(i, umi_cutoff=umi[2]))) %>%
  do.call(rbind,.) %>%
  write_tsv("data/enn.tsv")



# CellBender --------------------------------------------------------------

process_CB <- function(sample_name){
  print(str_glue('Processing CellBender output for sample {sample_name}'))
  
  # Import CellBender results
  ed <- Read10X_h5(str_glue('data/{sample_name}/comparison/cb_filtered.h5'))
  
  # Import raw counts
  seurat.object <- Read10X(str_glue('data/{sample_name}/outs/raw_feature_bc_matrix'))
  
  # Subset to cells called by CellBender
  seurat.object <- seurat.object[,colnames(ed)]
  seurat.object <- CreateSeuratObject(seurat.object, project = sample_name)
  
  # Add mitochondrial gene content
  if(sample_name=="MB"){
    seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, "^mt-")
  } else {
    seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, "^MT-")
  }

  # Add nuclear fraction
  seurat.object$nuclear_fraction <- nuclear_fraction_tags(bam = str_glue('data/{sample_name}/outs/possorted_genome_bam.bam'),
                                                          barcodes = colnames(seurat.object),
                                                          cores=16)$nuclear_fraction
  
  # Identify any empty droplets
  ed <- identify_empty_drops(select(seurat.object@meta.data, nuclear_fraction, nCount_RNA), umi_rescue = 0) %>%
    mutate(percent.mt = seurat.object$percent.mt)
  
  # Add other relevant info
  ed <- mutate(ed,
               method="CellBender",
               sample = sample_name,
               param = "default")
  return(ed)
}

# Save out
lapply(samples, process_CB) %>%
  do.call(rbind,.) %>%
  write_tsv("data/cb.tsv")



# Cell Ranger 6.1.1 -------------------------------------------------------

process_CR <- function(sample_name){
  print(str_glue('Processing Cell Ranger 6.1.1 output for sample {sample_name}'))
  
  # Import Cell Ranger 6.1.1 results
  ed <- Read10X_h5(str_glue('data/{sample_name}/cellranger6/{sample_name}/outs/filtered_feature_bc_matrix.h5'))
  
  # Import raw counts
  seurat.object <- Read10X(str_glue('data/{sample_name}/outs/raw_feature_bc_matrix'))
  
  # Subset to cells called by CellBender
  seurat.object <- seurat.object[,colnames(ed)]
  seurat.object <- CreateSeuratObject(seurat.object, project = sample_name)
  
  # Add mitochondrial gene content
  if(sample_name=="MB"){
    seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, "^mt-")
  } else {
    seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, "^MT-")
  }
  
  # Add nuclear fraction
  seurat.object$nuclear_fraction <- nuclear_fraction_tags(bam = str_glue('data/{sample_name}/outs/possorted_genome_bam.bam'),
                                                          barcodes = colnames(seurat.object),
                                                          cores=16)$nuclear_fraction
  
  # Identify any empty droplets
  ed <- identify_empty_drops(select(seurat.object@meta.data, nuclear_fraction, nCount_RNA), umi_rescue = 0) %>%
    mutate(percent.mt = seurat.object$percent.mt)
  
  # Add other relevant info
  ed <- mutate(ed,
               method="Cell_Ranger_6.1.1",
               sample = sample_name,
               param = "default")
  return(ed)
}

# Save out
lapply(samples, process_CR) %>%
  do.call(rbind,.) %>%
  write_tsv("data/cr.tsv")



# Code graveyard ----------------------------------------------------------

nf <- lapply(c("ed", "enn", "cb", "cr"),
       function(i) read_tsv(str_glue('data/{i}.tsv'))) %>%
  do.call(rbind,.)

nf.summary <- nf %>%
  filter(percent.mt<15) %>%
  group_by(sample, method, param, cell_status) %>%
  summarise(count=n()) %>%
  pivot_wider(names_from = cell_status, values_from = count) %>%
  mutate(cell_plus_empty = cell + empty_droplet,
         percent_empty = empty_droplet/cell_plus_empty*100) %>%
  replace_na(list(empty_droplet=0, cell_plus_empty=0, percent_empty=0))

# EmptyNN failed for samples HL and GBM with param = 500
nf.summary <- rbind(nf.summary,
data.frame(sample=c("GBM","HL"),
           method="EmptyNN",
           param="500",
           cell=0,
           empty_droplet=0,
           cell_plus_empty=0,
           percent_empty=0))

# Adjust for plotting
nf.summary  <- nf.summary %>%
  mutate(method_param = case_when(
    method == "CellBender" ~ "CellBender",
    method == "Cell_Ranger_6.1.1" ~ "Cell Ranger 6.1.1",
    method == "EmptyDrops" & param == "100" ~"EmptyDrops lower=100",
    method == "EmptyDrops" & param == "500" ~"EmptyDrops lower=500",
    method == "EmptyNN" & param == "100" ~"EmptyNN threshold=100",
    method == "EmptyNN" & param == "500" ~"EmptyNN threshold=500",
  )) %>%
  mutate(method_param = str_wrap(method_param, width = 10)) %>%
  mutate(sample_long = case_when(
    sample == "GBM" ~ "Glioblastoma",
    sample == "MB" ~ "Mouse brain",
    sample == "HL" ~ "Hodgkin's lymphoma",
    sample == "PBMC" ~ "PBMCs"
  )) %>%
  mutate(sample_long = factor(sample_long,
                              levels=c("Mouse brain",
                                       "Hodgkin's lymphoma",
                                       "Glioblastoma",
                                       "PBMCs")))
# Save out
saveRDS(nf.summary, "data/nf_summary.rds")
