# Script information -----------------------------------------------------------

# title: Plot supplementary figures
# author: Walter Muskovic
# date: 2021-08-26
# description: In this script we are going to create an additional supplementary
# figure



# Import libraries -------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggpubr)
})



# Load data ---------------------------------------------------------------

# Load expression data
load("data_track/qc_examples.Rdata")
pbmc <- Read10X("data/PBMC/outs/raw_feature_bc_matrix/")[,qc_examples$cell_barcode[qc_examples$sample=="PBMC"]]
gbm <- Read10X("data/GBM/outs/raw_feature_bc_matrix/")[,qc_examples$cell_barcode[qc_examples$sample=="GBM"]]
hl <- Read10X("data/HL/outs/raw_feature_bc_matrix/")[,qc_examples$cell_barcode[qc_examples$sample=="HL"]]
mb <- Read10X("data/MB/outs/raw_feature_bc_matrix/")[,qc_examples$cell_barcode[qc_examples$sample=="MB"]]

# Get proportion of UMIs mapped to lncRNAs
lncRNA <- data.frame(sample=c(rep("GBM", sum(qc_examples$sample=="GBM")),
                              rep("HL", sum(qc_examples$sample=="HL")),
                              rep("PBMC", sum(qc_examples$sample=="PBMC")),
                              rep("MB", sum(qc_examples$sample=="MB"))),
                     cell_barcode = c(colnames(gbm), colnames(hl), colnames(pbmc), colnames(mb)),
                     MALAT1_NEAT1 = c(colSums(gbm[c("MALAT1","NEAT1"),]),
                                colSums(hl[c("MALAT1","NEAT1"),]),
                                colSums(pbmc[c("MALAT1","NEAT1"),]),
                                colSums(mb[c("Malat1","Neat1"),])))
                      

# Add to qc_examples data frame
qc_examples <- merge(qc_examples, lncRNA, by = c("cell_barcode", "sample")) %>%
  mutate(MALAT1_NEAT1_percentage = 100*(MALAT1_NEAT1/umi_count))



# Plot --------------------------------------------------------------------

# Change names for plotting
qc_examples <-
  mutate(
    qc_examples,
    flag = case_when(
      flag == "damaged_cell" ~ "damaged cell",
      flag == "empty_droplet" ~ "empty droplet",
      TRUE                      ~ "cell"
    ),
    sample = case_when(
      sample == "GBM" ~ "Glioblastoma",
      sample == "HL" ~ "Hodgkin's Lymphoma",
      sample == "PBMC" ~ "PBMCs",
      sample == "MB" ~ "Mouse brain"
    )
  )

# Plot
p1 <- lapply(c("Mouse brain", "Hodgkin's Lymphoma", "Glioblastoma","PBMCs"), function(i){
  filter(qc_examples, sample==i) %>%
    ggplot(aes(x=flag, y=MALAT1_NEAT1_percentage, fill=flag)) +
    geom_boxplot() +
    facet_wrap(~sample) +
    scale_fill_manual(values=c("damaged cell"="#88419d",
                               "empty droplet"="#2166ac",
                               "cell"="#b2182b")) +
    theme(legend.position = "none") +
    xlab("") +
    ylab("MALAT1/NEAT1 % of UMIs")
}) %>%
ggarrange(plotlist = ., labels =letters[1:4], ncol = 2, nrow = 2)

# Save out figure
ggsave(filename = str_glue('figures/Figure_S3.png'),
       plot = p1,
       device = "png",
       width = 30, height = 20, units = "cm")

ggsave(filename = str_glue('figures/Figure_S3.pdf'),
       plot = p1,
       device = "pdf",
       width = 30, height = 20, units = "cm")

