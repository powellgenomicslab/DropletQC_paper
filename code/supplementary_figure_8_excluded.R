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

# Create plots for the mouse kidney data sets
ed <- readRDS("data/MK/combined.metadata.subset.rds")

# Barchart of cell status split by dissociation/preservation
ed.prop <- ed %>% 
  select(Preservation, dissociation_protocol, cell_status, accession) %>%
  group_by(Preservation, dissociation_protocol, cell_status, accession) %>%
  summarise(count=n()) %>%
  group_by(accession) %>%
  mutate(percentage = count / sum(count)*100)

# Plot
p1 <- ggplot(ed.prop, aes(fill=cell_status, y=percentage, x=accession)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("cell"="#88419d", "damaged cell"="#2166ac", "empty droplet"="#b2182b")) +
  facet_wrap(dissociation_protocol ~ Preservation, scales='free', labeller = label_wrap_gen(multi_line=FALSE))

# Save
ggsave(filename = "figures/Figure_S8_excluded.png",plot = p1,device = "png",width = 35, height = 15,units = "cm")
ggsave(filename = "figures/Figure_S8_excluded.pdf",plot = p1,device = "pdf",width = 35, height = 15,units = "cm")

