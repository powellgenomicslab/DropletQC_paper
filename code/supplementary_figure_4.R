# Script information -----------------------------------------------------------

# title: Plot supplementary figures
# author: Walter Muskovic
# date: 2021-01-28
# description: In this script we are going to create the supplementary figures 
# presented in the manuscript



# Import libraries -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(DropletQC)
  #library(ggExtra)
  #library(ggpubr)
})



# Load data --------------------------------------------------------------------

load("data_track/qc_examples.Rdata")

# Simplify unresolved clusters for plotting
qc_examples <- qc_examples %>%
  mutate(
    cell_type = replace(
      cell_type,
      cell_type == "astrocyte_unresolved_1",
      "astrocyte_unresolved"
    ),
    cell_type = replace(
      cell_type,
      cell_type == "astrocyte_unresolved_2",
      "astrocyte_unresolved"
    ),
    cell_type = replace(
      cell_type,
      cell_type == "neuron_unresolved_1",
      "neuron_unresolved"
    ),
    cell_type = replace(
      cell_type,
      cell_type == "neuron_unresolved_2",
      "neuron_unresolved"
    ),
    cell_type = replace(
      cell_type,
      cell_type == "neuron_unresolved_3",
      "neuron_unresolved"
    ),
    cell_type = replace(
      cell_type,
      cell_type == "neuron_unresolved_4",
      "neuron_unresolved"
    ),
    cell_type = replace(
      cell_type,
      cell_type == "neuron_unresolved_5",
      "neuron_unresolved"
    ),
    cell_type = replace(
      cell_type,
      cell_type == "neuron_unresolved_6",
      "neuron_unresolved"
    )
  )

# Wrap cell type names for plotting
qc_examples <- mutate(qc_examples,
                      cell_type = str_replace_all(cell_type, "_", " ")) %>%
  mutate(cell_type = str_wrap(cell_type, width = 15))



# Figure S3 --------------------------------------------------------------------

# Create plots showing how both the distributions of the nuclear fraction and 
# UMI counts can be quite different for different cell types in the same sample

supp3 <- function(sample_name, plot_title){
  p1 <- ggplot(filter(qc_examples, sample==sample_name),
              aes(x=cell_type, y=nuclear_fraction_droplet_qc)) + 
    geom_violin()+
    labs(title=plot_title,
         x="Cell type", 
         y="Nuclear fraction") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  
  p2 <- ggplot(filter(qc_examples, sample==sample_name),
               aes(x=cell_type, y=log10_umi_count)) + 
    geom_violin()+
    labs(title=plot_title,
         x="Cell type", 
         y="log10(UMI count)") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  return(p1 + p2)
}

# Save out plot
sample_names <- c("MB", "GBM", "PBMC", "HL")
plot_titles <- c("Mouse brain", "Glioblastoma", "PBMCs", "Hodgkin's lymphoma")
s1 <- map(1:4, function(x) supp3(sample_names[x], plot_titles[x]))
ggsave(filename = str_glue('figures/Figure_S4.png'),
       plot = s1[[1]] / s1[[2]] / s1[[3]] / s1[[4]],
       device = "png", width = 35, height = 30, units = "cm")
ggsave(filename = str_glue('figures/Figure_S4.pdf'),
       plot = s1[[1]] / s1[[2]] / s1[[3]] / s1[[4]],
       device = "pdf", width = 35, height = 30, units = "cm")
print("Finished creating supplementary figures")
