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
})



# Figure S6 --------------------------------------------------------------------

nf <- lapply(c("ed", "enn", "cb", "cr"),
             function(i) read_tsv(str_glue('data/{i}.tsv'))) %>%
  do.call(rbind,.) %>%
  filter(method=="EmptyDrops", param=="500") %>%
  mutate(`MT%`=percent.mt) %>%
  arrange(percent.mt)%>%
  mutate(
    sample_long = case_when(
      sample == "MB" ~ "Mouse brain",
      sample == "GBM" ~ "Glioblastoma",
      sample == "HL" ~ "Hodgkin's lymphoma",
      sample == "PBMC" ~ "PBMCs")
  )

p1 <- lapply(c("MB", "HL", "GBM", "PBMC"), function(i){
  filter(nf, sample == i) %>%
    ggplot(aes(x=nuclear_fraction, y=log10(nCount_RNA), colour=`MT%`)) +
    geom_point() + 
    labs(y="Log10(UMIs detected)", x = "Nuclear fraction") +
    xlim(range(nf$nuclear_fraction)) +
    facet_wrap(~ sample_long, labeller = label_wrap_gen(multi_line=FALSE), ncol = 2, scales = "free_y") +
    scale_colour_viridis_c(limits=c(0, 50), oob = scales::squish)
}) %>%
wrap_plots(ncol=2)+ plot_annotation(tag_levels = 'a')

ggsave("figures/Figure_S6.png",
       plot = p1,
       device = "png", width = 25,height = 20,units = "cm")

ggsave("figures/Figure_S6.pdf",
       plot = p1,
       device = "pdf", width = 25,height = 20,units = "cm")
