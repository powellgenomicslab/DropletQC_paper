# Script information -----------------------------------------------------------

# title: Plot figures
# author: Walter Muskovic
# date: 2020-12-18
# description: In this script we are going to create the figures presented
# in the manuscript

# We also create a data frame with:
# sample name
# cell barcodes
# UMAP_1 and UMAP_2 coords
# Seurat cluster info
# Cell type 
# the number of UMIs, add log10
# MT% - percentage of reads mapping to the mitochondrial genome
# emptyDrops info
# nuclear fraction score calculated with the dropletQC package
# nuclear fraction score calculated from the velocyto output
# flagged cells (one of three; "cell", "damaged_cell", "empty_droplet)


# Import libraries -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(janitor)
  library(ggExtra)
  library(ggpubr)
})



# QC data frame ----------------------------------------------------------------

# Create a data frame with just the QC info
getQC <- function(sample_name){
  print(str_glue("Reading in sample {sample_name}"))
  input.seurat  <- readRDS(
    str_glue('../data/{sample_name}/seurat_SCT_anno_QC_filter.rds'))
  input.seurat$sample <- sample_name
  seurat.df <- dplyr::select(input.seurat[[]], sample, cell_barcode,
                             UMAP_1, UMAP_2, seurat_clusters, cell_type,
                             UMI_count, log10_UMI_count, 
                             percent.mt, emptyDrops_LogProb, emptyDrops_PValue,
                             emptyDrops_FDR, nuclear_fraction_dropletQC,
                             nuclear_fraction_velocyto, flag)
  return(seurat.df)
}

qc_examples <- lapply(c("GBM","HL","PBMC","MB"), getQC) %>% 
  do.call(rbind,.) %>%
  clean_names()

# Save out
save(qc_examples, file = "../data_track/qc_examples.Rdata", compress = "xz")



# Figure 2 ----------------------------------------------------------------

# Create plots showing damaged cells and empty droplets for four different
# 10x Genomics datasets

# red - #b2182b
# purple - #88419d
# blue - #2166ac

# Change names for plotting
qc_examples <-
  mutate(
    qc_examples,
    flag = case_when(
      flag == "damaged_cell" ~ "damaged cell",
      flag == "empty_droplet" ~ "empty droplet",
      TRUE                      ~ "cell"
    )
  )

# QC plot split by cell type
cell_type_plot <- function(i) {
  
  plot_title <-
    c("Mouse brain", "Hodgkin's Lymphoma","Glioblastoma", "PBMCs")[i]
  sample_name <- c("MB", "HL", "GBM", "PBMC")[i]
  
  ctplot <- ggplot(
    filter(qc_examples, sample == sample_name),
    aes(x = nuclear_fraction_droplet_qc,
        y = log10_umi_count,
        col = flag)
  ) +
    geom_point() +
    scale_color_manual(values=c("#88419d", "#2166ac", "#b2182b"))+
    theme(legend.position = "left") +
    ggtitle(plot_title) + 
    labs(y="Log10(UMIs detected)", x = "Nuclear fraction") + 
    theme(legend.title=element_blank())
  
  ctplot <- ggMarginal(
    ctplot,
    type = c("density"),
    margins = c("both"),
    size = 5
  )
}

glist <- lapply(1:4, cell_type_plot)
combined_plots <- ggarrange(plotlist = glist, labels =  letters[1:4], ncol = 2, nrow = 2)
ggsave(filename = str_glue('../figures/Figure_2.png'), plot = combined_plots, device = "png", width = 30, height = 20, units = "cm")
ggsave(filename = str_glue('../figures/Figure_2.pdf'), plot = combined_plots, device = "pdf", width = 30, height = 20, units = "cm")
