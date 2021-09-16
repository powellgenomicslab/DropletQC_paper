# Script information -----------------------------------------------------------

# title: Plot supplementary figures
# author: Walter Muskovic
# date: 2021-08-26
# description: In this script we are going to create an additional supplementary
# figure



# Import libraries -------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(palettetown)
  library(patchwork)
})



# Figure S2 --------------------------------------------------------------------
# Plot the number of cells called by each method
nf.summary <- readRDS("data/nf_summary.rds")

p1 <- ggplot(nf.summary, 
             aes(x=method_param,
                 y=cell_plus_empty,
                 fill=sample_long)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~sample_long, scale="free_y") +
  theme(legend.position="none") +
  ylab("Number of cells called") +
  xlab("Method") +
  scale_fill_poke(pokemon = 'Quilava', spread = 4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

# Plot the number of percentage of empty droplets remaining
p2 <- p1 +
  aes(y=percent_empty) +
  ylab("Percentage of empty droplets")

ggsave("figures/Figure_S2.pdf",
       plot = wrap_plots(p1, p2, nrow = 2) + plot_annotation(tag_levels = 'a'),
       device = "pdf", width = 20,height = 25,units = "cm")
ggsave("figures/Figure_S2.png",
       plot = wrap_plots(p1, p2, nrow = 2) + plot_annotation(tag_levels = 'a'),
       device = "png", width = 20,height = 25,units = "cm")
