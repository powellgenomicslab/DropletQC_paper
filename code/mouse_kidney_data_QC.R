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


# Get mouse kidney sample metadata ----------------------------------------

# sample names
samples <- list.files("data/MK/GSE141115_RAW/") %>%
  str_remove_all(".txt.gz") %>%
  str_split_fixed("_",2) %>%
  data.frame() %>%
  `colnames<-`(c("accession", "library_short"))

# cell type info
metadata <- read.delim("data/MK/metadata.tsv", stringsAsFactors = FALSE)
metadata$library_short <- str_split_fixed(metadata$Library,"_",2)[,1]

# Add accession to cell type info
metadata <- merge(samples, metadata, by="library_short")
rm(samples)

# Add cell barcode short version with "-1" at the end
metadata$cell_barcode_short <- str_glue('{str_sub(metadata$cell_barcode, -16, -1)}-1')



# emptyDrops --------------------------------------------------------------

emptyDropsExtra <- function(sample_name, ed_seed = 123){
  
  # Import UMI counts for *all* barcodes from raw_feature_bc_matrix
  print(str_glue('Importing raw_feature_bc_matrix for sample {sample_name}'))
  tenx.counts <- Read10X(
    str_glue('data/MK/{sample_name}/raw_feature_bc_matrix/')
  )
  
  ## Diagnostic plot 1 of 2: barcode rank plot
  print(str_glue('Creating barcode rank plot for sample {sample_name}'))
  bcrank <- barcodeRanks(tenx.counts, lower = 100)
  # Only plot unique points (for plotting speed)
  uniq <- !duplicated(bcrank$rank)
  png(filename = str_glue('figures/barcode_rank_{sample_name}.png'),
      width = 15, height = 15, units = "cm", res=300)
  plot(bcrank$rank[uniq],
       bcrank$total[uniq],
       log="xy",
       xlab="Rank",
       ylab="Total UMI count",
       cex.lab=1.2)
  abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
  abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
  legend("bottomleft", legend=c("Inflection", "Knee"), 
         col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
  dev.off()
  
  ## Run emptryDrops()
  print(str_glue('Running empty drops for sample {sample_name}'))
  tic()
  set.seed(ed_seed)
  e.out <- emptyDrops(tenx.counts, lower = 100)
  is.cell <- e.out$FDR <= 0.01
  print(sum(is.cell, na.rm = TRUE))
  print(table(Limited=e.out$Limited, Significant=is.cell))
  # Restrict to barcodes called as cells by emptyDrops
  tenx.counts.filtered <- tenx.counts[,which(is.cell)]
  toc()
  
  ## Save out emptyDrops results
  saveRDS(e.out, str_glue('data/MK/{sample_name}/ED.rds'))
  
  ## Diagnostic plot 2 of 2: total count against the negative log-probability. 
  png(filename = str_glue('figures/UMI_logProb_QC_{sample_name}.png'),
      width = 15, height = 15, units = "cm", res=300)
  plot(e.out$Total,
       -e.out$LogProb,
       col=ifelse(is.cell, "red", "black"),
       xlab="Total UMI count", ylab="-Log Probability")
  dev.off()
  
  ## Write filtered barcodes out
  print(str_glue('Writing files to filtered_feature_bc_matrix directory for sample {sample_name}'))
  dir.create(str_glue('data/MK/{sample_name}/filtered_feature_bc_matrix'))
  write_tsv(x = data.frame(sort(colnames(tenx.counts.filtered))),
            file = str_glue('data/MK/{sample_name}/filtered_feature_bc_matrix/barcodes.tsv.gz'),
            col_names = FALSE) 
  
  ## Features are unchanged, so just copy
  file.copy(from = str_glue('data/MK/{sample_name}/raw_feature_bc_matrix/features.tsv.gz'),
            to = str_glue('data/MK/{sample_name}/filtered_feature_bc_matrix/features.tsv.gz'))
  
  ## Write filtered matrix
  # This is a bit convoluted - it's just to keep the two header lines
  writeMM(obj = tenx.counts.filtered, file = str_glue('data/MK/{sample_name}/filtered_feature_bc_matrix/matrix.mtx'), colnames= FALSE)
  mtx <- read_tsv(str_glue('data/MK/{sample_name}/filtered_feature_bc_matrix/matrix.mtx'), skip=1, col_names = FALSE)
  top_lines <- read_lines(str_glue('data/MK/{sample_name}/raw_feature_bc_matrix/matrix.mtx.gz'), n_max = 2)
  write_lines(top_lines, file = str_glue('data/MK/{sample_name}/filtered_feature_bc_matrix/matrix.mtx'))
  write.table(mtx, file = str_glue('data/MK/{sample_name}/filtered_feature_bc_matrix/matrix.mtx'), col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
  R.utils::gzip(str_glue('data/MK/{sample_name}/filtered_feature_bc_matrix/matrix.mtx'))
  
  print(str_glue('Completed running empty drops for sample {sample_name}'))
}

unique(metadata$accession) %>% walk(emptyDropsExtra)



# Nuclear fraction --------------------------------------------------------

for (sample_name in unique(metadata$accession)){
  print(str_glue('Calculating nuclear fraction for sample {sample_name}'))
  nf <- nuclear_fraction_tags(bam = str_glue('data/MK/{sample_name}/possorted_genome_bam.bam'),
                              barcodes = str_glue('data/MK/{sample_name}/filtered_feature_bc_matrix/barcodes.tsv.gz'),
                              cores = 4)
  write_tsv(nf, str_glue('data/MK/{sample_name}/nf.tsv'))
  }



# Combine QC info ---------------------------------------------------------

# Define function to combine QC info
extract_metadata<- function(sample_name){
  
  print(str_glue('Processing {sample_name}'))
  
  # Create Seurat object
  seurat.object <- Read10X(str_glue('data/MK/{sample_name}/filtered_feature_bc_matrix'))
  seurat.object <- CreateSeuratObject(seurat.object, project = sample_name)
  
  # Add the percentage of reads that map to the mitochondrial genome
  seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object,
                                                       pattern = "^mt-")  
  # Add emptyDrops info
  ed <- data.frame(readRDS(str_glue('data/MK/{sample_name}/ED.rds')))
  seurat.object <- AddMetaData(seurat.object, ed)

  # Add nuclear fraction
  nf <- read_tsv(str_glue('data/MK/{sample_name}/nf.tsv'))
  seurat.object$nuclear_fraction <- nf$nuclear_fraction
  
  # Add cell type and sample info
  metadata.filtered <- filter(metadata, accession==sample_name) %>% data.frame()
  row.names(metadata.filtered) <- metadata.filtered$cell_barcode_short
  seurat.object <- AddMetaData(seurat.object, metadata.filtered)
  
  # Return metadata
  return(seurat.object@meta.data)
}

# Combine QC info from all samples
combined.metadata <- lapply(unique(metadata$accession),extract_metadata)
combined.metadata <- do.call(rbind,combined.metadata)
saveRDS(combined.metadata, str_glue('data/MK/combined.metadata.rds'))



# Summarise QC info -------------------------------------------------------

combined.metadata <- readRDS("data/MK/combined.metadata.rds")

# Restrict to barcodes called as cells in the study
combined.metadata <- filter(combined.metadata, complete.cases(cell_type))

# Detect empty droplets and damaged cells
ed <- lapply(unique(combined.metadata$accession), function(i){
  
  print(i)
  
  # Get subset
  combined.metadata.subset <- filter(combined.metadata, accession==i)
  
  # Identify empty droplets
  ed <- identify_empty_drops(nf_umi = select(combined.metadata.subset, nuclear_fraction, nCount_RNA))
  
  # Add cell type info
  ed$cell_type <- combined.metadata.subset$cell_type
  
  # Remove "Unknown" or cell types with only 1 cell (that's not marked as an empty droplet)
  ed <- ed %>%
    rownames_to_column(var="cell_barcode") %>%
    group_by(cell_type,cell_status) %>%
    filter(n()>1) %>%
    filter(cell_type!="Unknown") %>%
    ungroup() %>%
    data.frame()
  ed <- data.frame(ed)
  row.names(ed) <- ed$cell_barcode
  ed <- select(ed, -cell_barcode)
  
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
  )

# Restrict to the experiment comparing two tissue dissociation protocols
# (warm/cold) combined with fresh/cryopreserved/methanol-fixed presrvation
table(select(ed, dissociation_protocol, Preservation), useNA = "always")
ed <- ed %>%
  filter(complete.cases(dissociation_protocol)) %>%
  replace_na(list(Preservation = "fresh"))
table(select(ed, dissociation_protocol, Preservation), useNA = "always")

# Save out this summary
saveRDS(ed, str_glue('data/MK/combined.metadata.subset.rds'))

