
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(Seurat)
library(dropletQC)
library(rhdf5)
library(ggpubr)



# Load data ---------------------------------------------------------------

# Get carcodes from study
Fresh_85T <- h5read("data/GSM4956405_85TFresh_combined.h5", "matrix")$barcodes
Fresh_86T <- h5read("data/GSM4956406_86TFresh_combined.h5", "matrix")$barcodes
Cryo_85T <- h5read("data/GSM4956407_C85T.h5", "matrix")$barcodes
Cryo_86T <- h5read("data/GSM4956408_C86T.h5", "matrix")$barcodes



# Calculate nuclear fraction ----------------------------------------------

Fresh_85T_nf <- nuclear_fraction_annotation(annotation_path = "data/Macaca_mulatta.Mmul_10.103.chr.gtf",
                                            bam = "data/alignments/85T_FreshAligned.sortedByCoord.out.bam",
                                            barcodes = str_remove_all(Fresh_85T,"-1"),
                                            cores = 16,
                                            verbose = TRUE)

Fresh_86T_nf <- nuclear_fraction_annotation(annotation_path = "data/Macaca_mulatta.Mmul_10.103.chr.gtf",
                                            bam = "data/alignments/86T_FreshAligned.sortedByCoord.out.bam",
                                            barcodes = str_remove_all(Fresh_86T,"-1"),
                                            cores = 16,
                                            verbose = TRUE)

Cryo_85T_nf <- nuclear_fraction_annotation(annotation_path = "data/Macaca_mulatta.Mmul_10.103.chr.gtf",
                                           bam = "data/alignments/85T_CryoAligned.sortedByCoord.out.bam",
                                           barcodes = str_remove_all(Cryo_85T,"-1"),
                                           cores = 16,
                                           verbose = TRUE)

Cryo_86T_nf <- nuclear_fraction_annotation(annotation_path = "data/Macaca_mulatta.Mmul_10.103.chr.gtf",
                                           bam = "data/alignments/86T_CryoAligned.sortedByCoord.out.bam",
                                           barcodes = str_remove_all(Cryo_86T,"-1"),
                                           cores = 16,
                                           verbose = TRUE)

# Get Seurat object
Fresh_85T_seurat <- Read10X("data/alignments/85T_FreshSolo.out/Gene/raw") %>% CreateSeuratObject()
Fresh_86T_seurat <- Read10X("data/alignments/86T_FreshSolo.out/Gene/raw") %>% CreateSeuratObject()
Cryo_85T_seurat <- Read10X("data/alignments/85T_CryoSolo.out/Gene/raw") %>% CreateSeuratObject()
Cryo_86T_seurat <- Read10X("data/alignments/86T_CryoSolo.out/Gene/raw") %>% CreateSeuratObject()

# Restrict to barcodes
Fresh_85T_seurat <- subset(Fresh_85T_seurat, cells = str_remove_all(Fresh_85T, "-1"))
Fresh_86T_seurat <- subset(Fresh_86T_seurat, cells = str_remove_all(Fresh_86T, "-1"))
Cryo_85T_seurat <- subset(Cryo_85T_seurat, cells = str_remove_all(Cryo_85T, "-1"))
Cryo_86T_seurat <- subset(Cryo_86T_seurat, cells = str_remove_all(Cryo_86T, "-1"))

# Add nuclear fraction
Fresh_85T_seurat <- AddMetaData(object = Fresh_85T_seurat, metadata = Fresh_85T_nf, col.name = "nf")
Fresh_86T_seurat <- AddMetaData(object = Fresh_86T_seurat, metadata = Fresh_86T_nf, col.name = "nf")
Cryo_85T_seurat <- AddMetaData(object = Cryo_85T_seurat, metadata = Cryo_85T_nf, col.name = "nf")
Cryo_86T_seurat <- AddMetaData(object = Cryo_86T_seurat, metadata = Cryo_86T_nf, col.name = "nf")



# Save out ----------------------------------------------------------------

saveRDS(object = Fresh_85T_seurat, file = "data/Fresh_85T_seurat.rds")
saveRDS(object = Fresh_86T_seurat, file = "data/Fresh_86T_seurat.rds")
saveRDS(object = Cryo_85T_seurat, file = "data/Cryo_85T_seurat.rds")
saveRDS(object = Cryo_86T_seurat, file = "data/Cryo_86T_seurat.rds")

Fresh_85T_seurat <- readRDS("data/Fresh_85T_seurat.rds")
Fresh_86T_seurat <- readRDS("data/Fresh_86T_seurat.rds")
Cryo_85T_seurat <- readRDS("data/Cryo_85T_seurat.rds")
Cryo_86T_seurat <- readRDS("data/Cryo_86T_seurat.rds")



# Identify empty drops ----------------------------------------------------

Fresh_85T_nf_umi <- select(Fresh_85T_seurat[[]], nf, nCount_RNA) %>%
  identify_empty_drops(include_plot = TRUE, umi_rescue = max(Fresh_85T_seurat$nCount_RNA)) %>%
  mutate(ct = "microglia")

Fresh_86T_nf_umi <- select(Fresh_86T_seurat[[]], nf, nCount_RNA) %>%
  identify_empty_drops(include_plot = TRUE, umi_rescue = max(Fresh_86T_seurat$nCount_RNA)) %>%
  mutate(ct = "microglia")

Cryo_85T_nf_umi <- select(Cryo_85T_seurat[[]], nf, nCount_RNA) %>%
  identify_empty_drops(include_plot = TRUE, umi_rescue = max(Cryo_85T_seurat$nCount_RNA)) %>%
  mutate(ct = "microglia")

Cryo_86T_nf_umi <- select(Cryo_86T_seurat[[]], nf, nCount_RNA) %>%
  identify_empty_drops(include_plot = TRUE, umi_rescue = max(Cryo_86T_seurat$nCount_RNA)) %>%
  mutate(ct = "microglia")



# Identify damaged cells --------------------------------------------------

Fresh_85T_nf_umi <- identify_damaged_cells(Fresh_85T_nf_umi, output_plots = TRUE)
Fresh_85T_nf_umi$plots

Fresh_86T_nf_umi <- identify_damaged_cells(Fresh_86T_nf_umi, output_plots = TRUE)
Fresh_86T_nf_umi$plots

Cryo_85T_nf_umi <- identify_damaged_cells(Cryo_85T_nf_umi, output_plots = TRUE)
Cryo_85T_nf_umi$plots

Cryo_86T_nf_umi <- identify_damaged_cells(Cryo_86T_nf_umi, output_plots = TRUE)
Cryo_86T_nf_umi$plots

# Combine results
cell_status <- list(
  data.frame(Sample = "Fresh_85T", cell_status = Fresh_85T_nf_umi$df$cell_status),
  data.frame(Sample = "Fresh_86T", cell_status = Fresh_86T_nf_umi$df$cell_status),
  data.frame(Sample = "Cryo_85T", cell_status = Cryo_85T_nf_umi$df$cell_status),
  data.frame(Sample = "Cryo_86T", cell_status = Cryo_86T_nf_umi$df$cell_status)
) %>%
  do.call(rbind,.)%>%
  group_by(Sample, cell_status) %>%
  summarise(Proportion=n())

# Print for future reference
cell_status %>%
  group_by(Sample) %>%
  mutate(percentage = Proportion / sum(Proportion) *100)

#   Sample    cell_status   Proportion percentage
#1 Cryo_85T  cell                4157     85.9  
#2 Cryo_85T  damaged_cell         650     13.4  
#3 Cryo_85T  empty_droplet         31      0.641
#4 Cryo_86T  cell                4200     85.3  
#5 Cryo_86T  damaged_cell         693     14.1  
#6 Cryo_86T  empty_droplet         28      0.569
#7 Fresh_85T cell                7953     93.8  
#8 Fresh_85T damaged_cell         458      5.40 
#9 Fresh_85T empty_droplet         68      0.802
#10 Fresh_86T cell                6085     96.0  
#11 Fresh_86T damaged_cell         180      2.84 
#12 Fresh_86T empty_droplet         73      1.15



# Plot --------------------------------------------------------------------

# Percent stacked barchart
p1 <- ggplot(cell_status, aes(fill=cell_status, y=Proportion, x=Sample)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("#88419d", "#2166ac", "#b2182b")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Nuclear fraction vs UMI count plots
p2 <- ggplot(Fresh_85T_nf_umi$df, aes(x = nf,y = log10(nCount_RNA), colour = cell_status)) +
  geom_point() +
  ggtitle("Fresh_85T") +
  scale_color_manual(values=c("#88419d", "#2166ac", "#b2182b")) +
  theme(legend.position = "none") +
  labs(y="Log10(UMIs detected)", x = "Nuclear fraction")

p3 <- ggplot(Fresh_86T_nf_umi$df, aes(x = nf,y = log10(nCount_RNA), colour = cell_status)) +
  geom_point() +
  ggtitle("Fresh_86T") +
  scale_color_manual(values=c("#88419d", "#2166ac", "#b2182b")) +
  theme(legend.position = "none")  +
  labs(y="Log10(UMIs detected)", x = "Nuclear fraction")

p4 <- ggplot(Cryo_85T_nf_umi$df, aes(x = nf,y = log10(nCount_RNA), colour = cell_status)) +
  geom_point() +
  ggtitle("Cryo_85T") +
  scale_color_manual(values=c("#88419d", "#2166ac", "#b2182b")) +
  theme(legend.position = "none")  +
  labs(y="Log10(UMIs detected)", x = "Nuclear fraction")

p5 <- ggplot(Cryo_86T_nf_umi$df, aes(x = nf,y = log10(nCount_RNA), colour = cell_status)) +
  geom_point() +
  ggtitle("Cryo_86T") +
  scale_color_manual(values=c("#88419d", "#2166ac", "#b2182b")) +
  theme(legend.position = "none")  +
  labs(y="Log10(UMIs detected)", x = "Nuclear fraction")

p_combined <- ggarrange(
  ggarrange(p2, p3, p4, p5, labels =  letters[1:4], ncol = 2, nrow = 2),
  ggarrange(p1, labels="e"), ncol=2, widths = c(1,0.5)
  )

# Save out plots
ggsave(filename = str_glue('figures/Figure_S2.pdf'),
       plot = p_combined,
       device = "pdf", width = 30, height = 20, units = "cm")

ggsave(filename = str_glue('figures/Figure_S2.png'),
       plot = p_combined,
       device = "png", width = 30, height = 20, units = "cm")
