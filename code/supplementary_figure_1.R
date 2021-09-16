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
  library(DropletQC)
  library(janitor)
  library(ggExtra)
  library(ggpubr)
  library(patchwork)
})

samples <- c("MB", "HL", "GBM", "PBMC")


# Load data ---------------------------------------------------------------

# Get 10,000 empty droplets excluded by EmptyDrops if not already done
if(!file.exists("data/spike.rds")){
  
  spike <- lapply(samples,
         function(sample_name){

  # import raw matrix
  counts.raw <- Read10X_h5(str_glue('data/{sample_name}/outs/raw_feature_bc_matrix.h5'))
  
  # Import barcodes called as cells
  ed <- read_tsv(str_glue('data/{sample_name}/comparison/barcodes_500.tsv'), col_names = "cb")
           
  # Take top 10,000 cell barcodes with; >= 100 counts and not called by
  # EmptyDrops as cells
  empty.droplets <- colSums(counts.raw) %>%
    sort(decreasing = TRUE) %>%
    .[.>=100] %>%
    names() %>%
    .[!.%in%ed$cb] %>%
    head(10000)
  
  # Get nuclear fraction for these cell barcodes
  nf <- nuclear_fraction_tags(bam = str_glue('data/{sample_name}/outs/possorted_genome_bam.bam'),
                              barcodes = empty.droplets,
                              cores = 16)$nuclear_fraction
  
  return(
    data.frame(sample=sample_name,
             cell_barcode=empty.droplets,
             log10_umi_count=log10(colSums(counts.raw[,empty.droplets])),
             nuclear_fraction_droplet_qc=nf)
    )
  }) %>%
    do.call(rbind,.)
  saveRDS(spike, "data/spike.rds")

} else {
  spike <- readRDS("data/spike.rds")
}

load("data_track/qc_examples.Rdata")
qc_examples <- bind_rows(qc_examples, spike)


# Figure S1  --------------------------------------------------------------

# Create density plots showing populations of;
  # cells  
  # damaged cells 
  # empty droplets called by DropletQC
  # empty droplets called by EmptyDrops
# for four different 10x Genomics datasets;
  # Mouse brain
  # Hodgkin's Lymphoma
  # Glioblastoma
  # PBMCs

# red - #b2182b
# purple - #88419d
# blue - #2166ac

# Change names for plotting
qc_examples <-
  mutate(
    qc_examples,
    flag = case_when(
      flag == "damaged_cell" ~ "damaged cell",
      flag == "empty_droplet" ~ "empty droplet - DropletQC",
      is.na(flag) ~ "empty droplet - EmptyDrops",
      TRUE                      ~ "cell"
    )
  ) %>%
  mutate(
      sample_long = case_when(
        sample == "MB" ~ "Mouse brain",
        sample == "GBM" ~ "Glioblastoma",
        sample == "HL" ~ "Hodgkin's lymphoma")
    )

# Plot
p1 <- lapply(c("MB", "GBM", "HL"), function(sample_name){
ggplot(filter(qc_examples, sample == sample_name, flag=="cell"),
  aes(x = nuclear_fraction_droplet_qc,
      y = log10_umi_count,
      fill = flag,
      colour =flag)) +
  geom_density_2d(bins=10) +
  geom_density_2d(data = filter(qc_examples, sample == sample_name, flag=="damaged cell"),
                  aes(x = nuclear_fraction_droplet_qc,
                      y = log10_umi_count,
                      colour =flag), bins=10)+
  geom_density_2d(data = filter(qc_examples, sample == sample_name, flag=="empty droplet - DropletQC"),
                  aes(x = nuclear_fraction_droplet_qc,
                      y = log10_umi_count,
                      colour =flag), bins=10) +
  geom_density_2d(data = filter(qc_examples, sample == sample_name, flag=="empty droplet - EmptyDrops"),
                  aes(x = nuclear_fraction_droplet_qc,
                      y = log10_umi_count,
                      colour =flag), bins=10) +
  scale_color_manual(values=c("#88419d", "#2166ac", "#b2182b", "#636363")) +
  facet_wrap(~sample_long) +
    theme(legend.title=element_blank()) +
    labs(y="Log10(UMIs detected)", x = "Nuclear fraction")
}) %>% wrap_plots(ncol=3, guides = "collect")

ggsave(filename = str_glue('figures/Figure_S1.png'), plot = p1, device = "png", width = 30, height = 10, units = "cm")
ggsave(filename = str_glue('figures/Figure_S1.pdf'), plot = p1, device = "pdf", width = 30, height = 10, units = "cm")
