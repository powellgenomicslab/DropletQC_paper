# Script information -----------------------------------------------------------

# title: Run emptyDrops
# author: Walter Muskovic
# date: 2021-01-20
# description: This script will run emptyDrops from the DropletUtils
# package for each of the 10x samples, creating some QC figures along the way.
# We also write out the filtered cells to a filtered_feature_bc_matrix, with the
# same strucure that would be output by Cell Ranger.



# Import libraries -------------------------------------------------------------

library(tidyverse)
library(DropletUtils)
library(Seurat)
library(tictoc)
library(Matrix)



# Run emptyDrops ---------------------------------------------------------------

# Create function to:
  # 1. Run emptyDrops 
  # 2. Save out some QC plots
  # 3. Save out the emptyDrops results
  # 4. Save out the cells that pass the emptyDrops threshold in a 
  # filtered_feature_bc_matrix directory, with the same structure to what we'd 
  # get from Cell Ranger


# Note that a lot of this code is identical to that from the DropletUtils
# vignette. Our only modification from the default parameters when using the
# emptyDrops() function is to increase the default of the `lower` argument from
# 100 to 500. Droplets with fewer UMIs than this are assumed to be empty
# droplets - the sequencing depth for all these libraries is pretty high.

emptyDropsExtra <- function(sample_name, ed_seed = 123){

  # Import UMI counts for *all* barcodes from raw_feature_bc_matrix
  print(str_glue('Importing raw_feature_bc_matrix for sample {sample_name}'))
  tenx.counts <- Read10X(
    str_glue('../data/{sample_name}/outs/raw_feature_bc_matrix/')
    )
  
  ## Diagnostic plot 1 of 2: barcode rank plot
  print(str_glue('Creating barcode rank plot for sample {sample_name}'))
  bcrank <- barcodeRanks(tenx.counts, lower = 500)
  # Only plot unique points (for plotting speed)
  uniq <- !duplicated(bcrank$rank)
  png(filename = str_glue('../figures/barcode_rank_{sample_name}.png'),
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
  e.out <- emptyDrops(tenx.counts, lower = 500)
  is.cell <- e.out$FDR <= 0.01
  print(sum(is.cell, na.rm = TRUE))
  print(table(Limited=e.out$Limited, Significant=is.cell))
  # Restrict to barcodes called as cells by emptyDrops
  tenx.counts.filtered <- tenx.counts[,which(is.cell)]
  toc()
  
  ## Save out emptyDrops results
  saveRDS(e.out, str_glue('../data/{sample_name}/ED.rds'))
  
  ## Diagnostic plot 2 of 2: total count against the negative log-probability. 
  png(filename = str_glue('../figures/UMI_logProb_QC_{sample_name}.png'),
      width = 15, height = 15, units = "cm", res=300)
  plot(e.out$Total,
       -e.out$LogProb,
       col=ifelse(is.cell, "red", "black"),
       xlab="Total UMI count", ylab="-Log Probability")
  dev.off()
  
  ## Write filtered barcodes out
  print(str_glue('Writing files to filtered_feature_bc_matrix directory for sample {sample_name}'))
  write_tsv(x = data.frame(sort(colnames(tenx.counts.filtered))),
            file = str_glue('../data/{sample_name}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'),
            col_names = FALSE) 
  
  ## Features are unchanged, so just copy
  file.copy(from = str_glue('../data/{sample_name}/outs/raw_feature_bc_matrix/features.tsv.gz'),
            to = str_glue('../data/{sample_name}/outs/filtered_feature_bc_matrix/features.tsv.gz'))
  
  ## Write filtered matrix
  # This is a bit convoluted - it's just to keep the two header lines
  writeMM(obj = tenx.counts.filtered, file = str_glue('../data/{sample_name}/outs/filtered_feature_bc_matrix/matrix.mtx'), colnames= FALSE)
  mtx <- read_tsv(str_glue('../data/{sample_name}/outs/filtered_feature_bc_matrix/matrix.mtx'), skip=1, col_names = FALSE)
  top_lines <- read_lines(str_glue('../data/{sample_name}/outs/raw_feature_bc_matrix/matrix.mtx.gz'), n_max = 2)
  write_lines(top_lines, file = str_glue('../data/{sample_name}/outs/filtered_feature_bc_matrix/matrix.mtx'))
  write.table(mtx, file = str_glue('../data/{sample_name}/outs/filtered_feature_bc_matrix/matrix.mtx'), col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
  R.utils::gzip(str_glue('../data/{sample_name}/outs/filtered_feature_bc_matrix/matrix.mtx'))
  
  print(str_glue('Completed running empty drops for sample {sample_name}'))
}

c("HL", "GBM", "PBMC", "MB", "GBM5p") %>% walk(emptyDropsExtra)
