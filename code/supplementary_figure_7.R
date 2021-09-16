# Script information -----------------------------------------------------------

# title: Plot supplementary figures
# author: Walter Muskovic
# date: 2021-08-26
# description: In this script we are going to create an additional supplementary
# figure



# Import libraries -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(dropletQC)
})



# Load data --------------------------------------------------------------------

gbm5p <- readRDS(str_glue('data/GBM5p/seurat_SCT_anno_QC_filter.rds'))
gbm5p <- gbm5p[[]] %>% arrange(desc(flag))

ms <- readRDS(str_glue('data/MS/seurat_SCT_anno.rds'))
ms <- ms[[]] %>% arrange(cell_status)



# Percentages -------------------------------------------------------------

# Get percentages for supplementary table

# GBM
table(gbm5p$flag)/nrow(gbm5p) * 100
#         cell  damaged_cell empty_droplet 
#92.448513      4.424104      3.127384 

# MS
table(ms$cell_status)/nrow(ms) * 100
#cell  damaged_cell empty_droplet 
#89.939252      7.641479      2.419269 


# Figure S7 --------------------------------------------------------------------

# Create plots showing damaged cells and empty droplets for the GBM 5' dataset
# and mouse spleen 5' dataset

p1 <- mutate(
  gbm5p,
  flag = case_when(
    flag == "damaged_cell" ~ "damaged cell",
    flag == "empty_droplet" ~ "empty droplet",
    TRUE                      ~ "cell"
  )
) %>%
  ggplot(aes(x = nuclear_fraction_dropletQC,
                                  y = log10(nCount_RNA), colour = flag)) +
  geom_point() +
  ggtitle("Glioblastoma 5'v1") +
  scale_color_manual(values=c("#88419d", "#2166ac", "#b2182b")) +
  labs(y="Log10(UMIs detected)", x = "Nuclear fraction") +
  theme(legend.position = "left") +
  theme(legend.title=element_blank())

p2 <- mutate(
  ms,
  cell_status = case_when(
    cell_status == "damaged_cell" ~ "damaged cell",
    cell_status == "empty_droplet" ~ "empty droplet",
    TRUE                      ~ "cell"
  )
) %>%
  ggplot(aes(x = nuclear_fraction,
                        y = log10(nCount_RNA), colour = cell_status)) +
  geom_point() +
  ggtitle("Mouse Splenocytes 5'v2") +
  scale_color_manual(values=c("#88419d", "#2166ac", "#b2182b")) +
  labs(y="Log10(UMIs detected)", x = "Nuclear fraction") +
  theme(legend.position = "left") +
  theme(legend.title=element_blank())

ggsave(filename = str_glue('figures/Figure_S3.png'), plot = wrap_plots(p1,p2, ncol = 2, guides = "collect") + plot_annotation(tag_levels = 'a'), device = "png", width = 25, height = 10, units = "cm")
ggsave(filename = str_glue('figures/Figure_S3.pdf'), plot = wrap_plots(p1,p2, ncol = 2, guides = "collect") + plot_annotation(tag_levels = 'a'), device = "pdf", width = 25, height = 10, units = "cm")
