# Script information ------------------------------------------------------

# title: Prepare Allen Brain Atlas mouse and developing mouse brain 
# reference dataset 
# author: Walter Muskovic
# date: 2020-12-16
# description: Prepare Allen Brain Atlas Mouse reference.From the atlas:
# This data set includes single-cell transcriptomes from multiple cortical
# areas and the hippocampal formation, including 1,093,785 total cells.
# Samples were collected from dissections of brain regions from 
# ~8 week-old male and female mice, from pan-neuronal transgenic lines. 
# In addition to the Allen Brain Atlas, we also download and prepare a 
# reference of developing mouse brain from Rosenberg et al, Science, 2018


# Import libraries --------------------------------------------------------

# Load required R packages
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(scPred)
  library(R.matlab)
  })


# Download ----------------------------------------------------------------

# Define directories
data_dir <- c("../data/reference")

print('Downloading Allen Brain Map mouse scRNA-seq dataset')

if(!file.exists(str_glue('{data_dir}/matrix.csv'))){
  metadata <- "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/metadata.csv"
  gene_expression_matrix <- "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/matrix.csv"
  download.file(metadata, destfile = str_glue('{data_dir}/metadata.csv'))
  download.file(gene_expression_matrix, destfile = str_glue('{data_dir}/matrix.csv'))
}



# Preprocessing -----------------------------------------------------------

# Import metadata
metadata <- read_csv(str_glue('{data_dir}/metadata.csv'))

# Peek
subclass_summary <- data.frame(table(metadata$subclass_label)) %>% arrange(Freq)
subclass_summary

# Remove NA and Meis2, which has only one cell
metadata <- filter(metadata, subclass_label!="Meis2")

# Downsample each subclass to have no more than 200 cells
no_downsampling <- filter(metadata, subclass_label%in%subclass_summary$Var1[subclass_summary$Freq<200])
downsampling <- filter(metadata, !subclass_label%in%subclass_summary$Var1[subclass_summary$Freq<200]) %>%
  group_by(subclass_label) %>%
  sample_n(200, replace = FALSE)
metadata_downsampled <- rbind(no_downsampling, downsampling)

# Import expression data and restrict to the downsampled cells
aba.data <- data.table::fread(str_glue('{data_dir}/matrix.csv'))
aba.data <- aba.data[aba.data$sample_name%in%metadata_downsampled$sample_name,]

# Rotate
aba.data.t <- t(aba.data[,2:ncol(aba.data)])
colnames(aba.data.t) <- aba.data$sample_name

# Create Seurat object
aba <- CreateSeuratObject(counts = aba.data.t,
                          project = "Allen_brain_atlas_mouse")

# Add subclass info as metadata
aba$subclass_label <- metadata_downsampled$subclass_label[match(colnames(aba),metadata_downsampled$sample_name)]
aba$subclass_color <- metadata_downsampled$subclass_color[match(colnames(aba),metadata_downsampled$sample_name)]

# Normalise with SCTransform and carry out standard steps for visualisation and
# clustering
aba <- SCTransform(aba, return.only.var.genes = FALSE) # 26 minutes
aba <- RunPCA(aba, verbose = TRUE)
aba <- RunUMAP(aba, dims = 1:50, verbose = TRUE)
aba <- FindNeighbors(aba, dims = 1:50, verbose = TRUE)
aba <- FindClusters(aba, verbose = TRUE)
saveRDS(aba, str_glue('{data_dir}/Allen_Brain_mouse_seurat.rds'))

# Create overview plot if not already created
ggsave(file = "../figures/Allen_Brain_mouse.pdf",
         plot = DimPlot(aba, group.by = "subclass_label", label = TRUE, repel = TRUE),
         device = "pdf",
         width = 35, height = 25, units = "cm")

# Prepare model to use with scPred for cell type classification
# Get feature space
reference <- getFeatureSpace(aba, "subclass_label")
# Train classifier
#reference <- trainModel(reference)
reference <- trainModel(reference, allowParallel = TRUE)
# save out
saveRDS(reference, str_glue('{data_dir}/Allen_Brain_mouse_scPred_model.rds'))



# Devloping mouse brain ---------------------------------------------------

# Here we use a reference dataset produced by Rosenberg et al in their 2018
# science paper "Single-cell profiling of the developing mouse brain and
# spinal cord with split-pool barcoding" - DOI: 10.1126/science.aam8999

# Download GSM3017261_150000_CNS_nuclei.mat from GSE110823
dev <- readMat("../reference/GSM3017261_150000_CNS_nuclei.mat")

# Downsample to no more than 200 cells for each subtype
dev.metadata <- tibble(sample_type = str_trim(dev$sample.type),
                       cell_type = str_replace_all(str_trim(str_split_fixed(dev$cluster.assignment, " ",2)[,2]), " ", "_")) %>%
  mutate(id = row_number()) %>%
  filter(sample_type %in% c("p2_brain", "p11_brain"))
# Get number of each cell type
cell_type_summary <- data.frame(table(dev.metadata$cell_type))  %>% arrange(Freq)
# Remove SC_Glut_Gna14 & SC_Glut_Gna14 (spinal cord), then downsample
dev.metadata <- filter(dev.metadata, cell_type!="SC_Glut_Gna14")
dev.metadata <- filter(dev.metadata, cell_type!="SC_Glut_Hmga2")

no_downsampling <- filter(dev.metadata, cell_type%in%cell_type_summary$Var1[cell_type_summary$Freq<200])
downsampling <- filter(dev.metadata, !cell_type%in%cell_type_summary$Var1[cell_type_summary$Freq<200]) %>%
  group_by(cell_type) %>%
  sample_n(200, replace = FALSE)
metadata_downsampled <- rbind(no_downsampling, downsampling)

# Get expression matrix for the downsampled cells
dev.data <- dev$DGE[metadata_downsampled$id,] %>% t() %>% as.matrix()
# Add cell names and gene names
colnames(dev.data) <- dev$barcodes[1,metadata_downsampled$id]
row.names(dev.data) <- str_trim(dev$genes)

# Create Seurat object
dev.seurat <- CreateSeuratObject(counts = dev.data, project = "Rosenberg_et_al")

# Simplify cell_type classification and add to the Seurat object as metadata
cell_types <- metadata_downsampled$cell_type
cell_types[str_detect(cell_types, "Astro_")] <- "Astrocyte"
cell_types[str_detect(cell_types, "CB_")] <- "Neuron_cerebellum"
cell_types[str_detect(cell_types, "CTX_")] <- "Neuron_cortex"
cell_types[str_detect(cell_types, "HIPP_")] <- "Neuron_hippocampus"
cell_types[str_detect(cell_types, "MD_")] <- "Neuron_medulla"
cell_types[str_detect(cell_types, "Migrating_")] <- "Migrating_interneurons"
cell_types[str_detect(cell_types, "SVZ_Stem")] <- "Migrating_interneurons"
cell_types[str_detect(cell_types, "OB_")] <- "Neuron_olfactory_bulb"
cell_types[str_detect(cell_types, "Oligo_")] <- "Oligodendrocyte"
cell_types[str_detect(cell_types, "Purkinje_")] <- "Neuron_cerebellum"
cell_types[str_detect(cell_types, "SUB_Pyr")] <- "Neuron_hippocampus"
cell_types[str_detect(cell_types, "THAL_")] <- "Neuron_thalamus"
cell_types[str_detect(cell_types, "Unresolved")] <- "Unresolved"
cell_types[str_detect(cell_types, "VLMC_")] <- "VLMC"
cell_types[str_detect(cell_types, "CLAU_Pyr")] <- "Neuron_cortex"
cell_types[str_detect(cell_types, "Medium_Spiny_Neurons")] <- "Neuron_striatum"
cell_types[str_detect(cell_types, "MTt_Glut")] <- "Neuron_rostral_midbrain"
cell_types[str_detect(cell_types, "Nigral_Dopaminergic")] <- "Neuron_basal_ganglia"

dev.seurat$cell_type <- cell_types

# Normalise with SCTransform and carry out standard steps for visualisation and
# clustering
dev.seurat <- SCTransform(dev.seurat, return.only.var.genes = FALSE)
dev.seurat <- RunPCA(dev.seurat, verbose = TRUE)
dev.seurat <- RunUMAP(dev.seurat, dims = 1:50, verbose = TRUE)
dev.seurat <- FindNeighbors(dev.seurat, dims = 1:50, verbose = TRUE)
dev.seurat <- FindClusters(dev.seurat, verbose = TRUE)
saveRDS(dev.seurat, str_glue('{data_dir}/developing_mouse_brain_seurat.rds'))

# Create overview plot if not already created
ggsave(file = "../figures/developing_mouse_brain.pdf",
       plot = DimPlot(dev.seurat, group.by = "cell_type",
                      label = TRUE, repel = TRUE) + NoLegend(),
       device = "pdf",
       width = 35, height = 35, units = "cm")

# Prepare model to use with scPred for cell type classification
# Get feature space
reference <- getFeatureSpace(dev.seurat, "cell_type")
# Train classifier
reference <- trainModel(reference)
# save out
saveRDS(reference, str_glue('{data_dir}/developing_mouse_brain_scPred_model.rds'))
