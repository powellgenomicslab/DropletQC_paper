# Load packages -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(Seurat)
library(DropletQC)
library(DropletUtils)
library(tictoc)
library(Matrix)
library(UCell)
library(patchwork)



# Define samples and directories ------------------------------------------

samples <- c("healthy","proapoptotic","apoptotic","notsorted")



# Download data -----------------------------------------------------------

data/APO
download.file(url = "ftp.sra.ebi.ac.uk/vol1/run/ERR337/ERR3379692/apoptotic_genome.bam",
              destfile = "data/APO/apoptotic_genome.bam")
download.file(url = "ftp.sra.ebi.ac.uk/vol1/run/ERR337/ERR3379693/healthy_genome.bam",
              destfile = "data/APO/healthy_genome.bam")
download.file(url = "ftp.sra.ebi.ac.uk/vol1/run/ERR337/ERR3379694/notsorted_genome.bam",
              destfile = "data/APO/notsorted_genome.bam")
download.file(url = "ftp.sra.ebi.ac.uk/vol1/run/ERR337/ERR3379695/proapoptotic_genome.bam",
              destfile = "data/APO/proapoptotic_genome.bam")



# emptyDrops --------------------------------------------------------------

emptyDropsExtra <- function(sample_name, ed_seed = 123){
  
  # Import UMI counts for *all* barcodes from raw_feature_bc_matrix
  print(str_glue('Importing raw_feature_bc_matrix for sample {sample_name}'))
  tenx.counts <- Read10X(
    str_glue('data/APO/{sample_name}_outs/outs/raw_feature_bc_matrix/')
  )
  
  ## Run emptryDrops()
  print(str_glue('Running empty drops for sample {sample_name}'))
  tic()
  set.seed(ed_seed)
  e.out <- emptyDrops(tenx.counts, lower = 700)
  is.cell <- e.out$FDR <= 0.01
  print(sum(is.cell, na.rm = TRUE))
  print(table(Limited=e.out$Limited, Significant=is.cell))
  # Restrict to barcodes called as cells by emptyDrops
  tenx.counts.filtered <- tenx.counts[,which(is.cell)]
  toc()
  
  ## Save out emptyDrops results
  saveRDS(e.out, str_glue('data/APO/{sample_name}/ED.rds'))
  
  ## Write filtered barcodes out
  print(str_glue('Writing files to filtered_feature_bc_matrix directory for sample {sample_name}'))
  dir.create(str_glue('data/APO/{sample_name}/filtered_feature_bc_matrix'))
  write_tsv(x = data.frame(sort(colnames(tenx.counts.filtered))),
            file = str_glue('data/APO/{sample_name}/filtered_feature_bc_matrix/barcodes.tsv.gz'),
            col_names = FALSE) 
  
  ## Features are unchanged, so just copy
  file.copy(from = str_glue('data/APO/{sample_name}_outs/outs/raw_feature_bc_matrix/features.tsv.gz'),
            to = str_glue('data/APO/{sample_name}/filtered_feature_bc_matrix/features.tsv.gz'))
  
  ## Write filtered matrix
  # This is a bit convoluted - it's just to keep the two header lines
  writeMM(obj = tenx.counts.filtered, file = str_glue('data/APO/{sample_name}/filtered_feature_bc_matrix/matrix.mtx'), colnames= FALSE)
  mtx <- read_tsv(str_glue('data/APO/{sample_name}/filtered_feature_bc_matrix/matrix.mtx'), skip=1, col_names = FALSE)
  top_lines <- read_lines(str_glue('data/APO/{sample_name}_outs/outs/raw_feature_bc_matrix/matrix.mtx.gz'), n_max = 2)
  write_lines(top_lines, file = str_glue('data/APO/{sample_name}/filtered_feature_bc_matrix/matrix.mtx'))
  write.table(mtx, file = str_glue('data/APO/{sample_name}/filtered_feature_bc_matrix/matrix.mtx'), col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
  R.utils::gzip(str_glue('data/APO/{sample_name}/filtered_feature_bc_matrix/matrix.mtx'))
  
  print(str_glue('Completed running empty drops for sample {sample_name}'))
}

samples %>% walk(emptyDropsExtra)



# Nuclear fraction --------------------------------------------------------

for (sample_name in samples){
  print(str_glue('Calculating nuclear fraction for sample {sample_name}'))
  nf <- nuclear_fraction_tags(bam = str_glue('data/APO/{sample_name}_outs/outs/possorted_genome_bam.bam'),
                              barcodes = str_glue('data/APO/{sample_name}/filtered_feature_bc_matrix/barcodes.tsv.gz'),
                              cores = 12)
  write_tsv(nf, str_glue('data/APO/{sample_name}/nf.tsv'))
  }



# Combine QC info ---------------------------------------------------------

# Define function to combine QC info
extract_metadata<- function(sample_name){
  
  print(str_glue('Processing {sample_name}'))
  
  # Create Seurat object
  seurat.object <- Read10X(str_glue('data/APO/{sample_name}/filtered_feature_bc_matrix'))
  seurat.object <- CreateSeuratObject(seurat.object, project = sample_name)
  
  # Add the percentage of reads that map to the mitochondrial genome
  seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object,
                                                       pattern = "^MT-")  
  # Add emptyDrops info
  ed <- data.frame(readRDS(str_glue('data/APO/{sample_name}/ED.rds')))
  seurat.object <- AddMetaData(seurat.object, ed)

  # Add nuclear fraction
  nf <- read_tsv(str_glue('data/APO/{sample_name}/nf.tsv'))
  seurat.object$nuclear_fraction <- nf$nuclear_fraction
  
  # Return metadata
  return(seurat.object@meta.data)
}

# Combine QC info from all samples
combined.metadata <- lapply(samples, extract_metadata)
combined.metadata <- do.call(rbind, combined.metadata)
saveRDS(combined.metadata, str_glue('data/APO/combined.metadata.rds'))



# Summarise QC info -------------------------------------------------------

combined.metadata <- readRDS("data/APO/combined.metadata.rds")

# Detect empty droplets and damaged cells
ed <- lapply(samples, function(i){
  
  print(i)
  
  # Get subset
  combined.metadata.subset <- filter(combined.metadata, orig.ident==i)
  
  # Identify empty droplets
  ed <- identify_empty_drops(nf_umi = select(combined.metadata.subset, nuclear_fraction, nCount_RNA))
  
  # Add cell type info
  ed$cell_type <- "HEK293"
  
  # Identify damaged cells
  ed <- identify_damaged_cells(nf_umi_ed_ct = ed, verbose = FALSE)$df
  ed <- select(ed, cell_status)
  
  # Add to other metadata
  de <- merge(combined.metadata.subset, ed, by=0)
  
  return(de)
  
}) %>%
  do.call(rbind,.) %>%
  mutate(
    cell_status = case_when(
      cell_status == "damaged_cell" ~ "damaged cell",
      cell_status == "empty_droplet" ~ "empty droplet",
      TRUE                      ~ "cell"
    )
  ) %>%
  dplyr::rename('MT%' = 'percent.mt')


# Cannot identify damaged cells in APO sample, as all cells appear to be damaged
ed$cell_status[ed$cell_status=="cell"&ed$orig.ident=="apoptotic"] <- "damaged cell"

# Save out this summary
saveRDS(ed, str_glue('data/APO/combined.metadata.subset.rds'))

