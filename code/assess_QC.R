# Script information -----------------------------------------------------------

# title: Assess QC metrics
# author: Walter Muskovic
# date: 2020-12-16
# description: In this script we are starting with a dataset that has already 
# been filtered by:
  # 1. excluding cells with >15% mitochondrial gene content
  # 2. removed empty droplets with DropletUtils::emptyDrops
# Now we use the nuclear fraction score to demonstrate how damaged cells or 
# droplets containing mostly ambient RNA can be identified



# Import libraries -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(ggExtra)
  library(ggpubr)
  library(DropletQC)
})



# Define directories & import data ---------------------------------------------

# Get the sample from the command line argument (integer 1-4)
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- c("GBM", "HL", "PBMC", "MB", "GBM5p")[i]

# Import Seurat
input.seurat <- readRDS(str_glue('../data/{sample_name}/seurat_SCT_anno_QC.rds'))



# Identify empty droplets -------------------------------------------------

if(sample_name!="GBM5p"){
  ed <- identify_empty_drops(nf_umi = data.frame(nf=input.seurat$nuclear_fraction_dropletQC,
                                                 umi=input.seurat$nCount_RNA))
}

# There are significantly fewer reads mapped to intronic regions in the 5p
# dataset: 5.6% compared to 24.9% for the 3' dataset. So the nuclear fraction
# scores are lower and there's insufficient separation between cells and empty
# droplets, so the automatic cut-off is inappropriate. Use the rescue parameters
# to deal with this.
if(sample_name=="GBM5p"){
  ed <- identify_empty_drops(nf_umi = data.frame(nf=input.seurat$nuclear_fraction_dropletQC,
                                                 umi=input.seurat$nCount_RNA),
                             nf_rescue = 0.015, umi_rescue = 100)
  }



# Identify damaged cells --------------------------------------------------

ed$ct <- input.seurat$cell_type
ed <- identify_damaged_cells(nf_umi_ed_ct = ed)
ed <- ed$df

# Add to Seurat metadata
input.seurat$flag <- ed$cell_status
# Get metadata
seurat_metadata <- input.seurat[[]]



# Plot nuclear fraction filtering ----------------------------------------------

# Plot the NMF ambient pattern vs the nuclear fraction for all of the cells in the sample
g1 <- ggplot(seurat_metadata,
             aes(x=nuclear_fraction_dropletQC, y=log10_UMI_count)) +
  geom_point() +
  theme(legend.position="left") +
  ggtitle(sample_name) +
  xlim(c(0,1))
g1 <- ggMarginal(g1, type = c("density"), margins = c("both"), size = 5)

# Now the same plot, this time split by cell type
g2 <- ggplot(seurat_metadata,
             aes(x=nuclear_fraction_dropletQC,
                 y=log10_UMI_count, col=cell_type)) +
  geom_point() +
  theme(legend.position="left") +
  ggtitle("Split by cell type") +
  xlim(c(0,1))
g2 <- ggMarginal(g2, type = c("density"),
                 margins = c("both"),
                 size = 5, groupColour = TRUE)

# Now individual plots for each cell type
cell_type_plot <- function(ct){
  ctplot <- ggplot(filter(seurat_metadata, cell_type==ct),
                   aes(x=nuclear_fraction_dropletQC,
                       y=log10_UMI_count, col=flag)) +
    geom_point() +
    theme(legend.position="left") +
    ggtitle(ct) +
    ylim(range(seurat_metadata$log10_UMI_count)) +
    xlim(c(0,1))
  ctplot <- ggMarginal(ctplot,
                       type = c("density"), margins = c("both"), size = 5)
}

# Get a UMAP as well so that we can visualise which cells have been excluded
g8 <- DimPlot(input.seurat, group.by = c("flag"))

# We use ggpubr::ggarrange() instead of patchwork here to arrange the plots,
# because patchwork doesn't play nicely with ggpubr::ggMarginal() which we
# used to add the density plots in the margins of the scatter plots
combined_plots <- append(append(list(g1, g2),
                                lapply(unique(seurat_metadata$cell_type),
                                       cell_type_plot)),
                         list(g8))
combined_plots <- ggarrange(plotlist = combined_plots,
                            labels =  letters[1:(length(unique(seurat_metadata$cell_type))+3)],
                            ncol = 4,
                            nrow = ceiling((length(unique(seurat_metadata$cell_type))+3)/4))
ggsave(filename = str_glue('../figures/{sample_name}_filtering.pdf'),
       plot = combined_plots, device = "pdf",
       width = 80,
       height = ceiling((length(unique(seurat_metadata$cell_type))+3)/4)*15, units = "cm")



# Save object ------------------------------------------------------------------

# We added to the Seurat object metadata a vector (flag), defining what we 
# suspect are damaged cells and empty droplets

# Save out the updated object
saveRDS(input.seurat,
        str_glue('../data/{sample_name}/seurat_SCT_anno_QC_filter.rds'))

print(str_glue('Completed assessing QC metrics for sample: {sample_name}'))
