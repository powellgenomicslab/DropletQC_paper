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



# Figure S5 --------------------------------------------------------------------

# Create plots for the HEK293 datasets

ed <- readRDS("data/APO/combined.metadata.subset.rds")
samples <- c("healthy","proapoptotic","apoptotic","notsorted")

# Plot samples
p1 <- lapply(samples[1:3], function(i){
  filter(ed, orig.ident == i) %>%
    ggplot(aes(x=nuclear_fraction, y=log10(nCount_RNA), colour=`MT%`)) +
    geom_point() + 
    labs(y="Log10(UMIs detected)", x = "Nuclear fraction") +
    labs(fill = "Dose (mg)") + 
    xlim(range(ed$nuclear_fraction)) +
    facet_wrap(~ orig.ident, labeller = label_wrap_gen(multi_line=FALSE), ncol = 2, scales = "free_y") +
    scale_colour_viridis_c()
})

# Barchart of cell status split by sample
ed.prop <- ed %>% 
  filter(orig.ident %in% samples[1:3]) %>%
  select(orig.ident, cell_status) %>%
  group_by(orig.ident, cell_status) %>%
  summarise(count=n()) %>%
  group_by(orig.ident) %>%
  mutate(Percentage = count / sum(count)*100)
# Plot
p1[[4]] <- ggplot(ed.prop, aes(fill=cell_status, y=Percentage, x=orig.ident)) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Sample") +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("cell"="#88419d", "damaged cell"="#2166ac", "empty droplet"="#b2182b"))

# MT% for each sample, split by cell status
p1[[5]] <-  ed %>% 
  filter(orig.ident!="notsorted") %>% 
  ggplot(aes(x=factor(orig.ident), y=`MT%`, fill=as.factor(cell_status))) +
  geom_violin() +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("cell"="#88419d", "damaged cell"="#2166ac", "empty droplet"="#b2182b")) +
  xlab("Sample")

# Combine plots
combined_plots <- (p1[[1]] | p1[[2]]) /
(p1[[3]] | p1[[4]]) /
p1[[5]] + plot_annotation(tag_levels = 'a')
combined_plots + plot_annotation()

# Save out plots
ggsave("figures/Figure_S5.png",
       plot = combined_plots,
       device = "png", width = 25,height = 25,units = "cm")

ggsave("figures/Figure_S5.pdf",
       plot = combined_plots,
       device = "pdf", width = 25,height = 25,units = "cm")
