# Key findings ------------------------------------------------------------

# This study:
# compares two tissue dissociation protocols
# compares fresh cells to cryopreserved and methanol-fixed cells

# Findings:
# digestion on ice avoids the stress response observed with 37C dissociation

# cell types more abundant either in the cold or warm dissociations that may
# represent populations that require gentler or harsher conditions to be
# released intact

# cryopreservation of dissociated cells results in a major loss of epithelial
# cell types

# methanol fixation maintains the cellular composition but suffers from ambient
# RNA leakage.

# cell type composition differences are observed between single-cell and
# single-nucleus RNA sequencing libraries



# Load packages -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(Seurat)



# Samples -----------------------------------------------------------------


# 77,656 single-cell, 98,303 single-nucleus, and 15 bulk RNA-seq profiles were
# generated and made publicly available (GSE141115).

# Get type labels from the paper
dir.create("data/MK")
download.file(url = "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-020-02048-6/MediaObjects/13059_2020_2048_MOESM3_ESM.xlsx",
              destfile = "data/MK/mk_cell_types.xlsx")

# Read cell type info in and tidy
mk <- lapply(excel_sheets("data/MK/mk_cell_types.xlsx"), function(i)
  read_excel("data/MK/mk_cell_types.xlsx",
             sheet = i)) %>%
  bind_rows() %>%
  rename(cell_barcode = `...1`,
         dissociation_protocol =`Dissociation protocol`,
         cell_type = `Cell type`)

# Check the number of nuclei
str_detect(mk$cell_barcode, "SN") %>% table()
# Great. This corresponds to the "77,656 single-cell, 98,303 single-nucleus"

# Remove nuclei
mk <- mk[!str_detect(mk$cell_barcode, "SN"),]

# Write out 
write_tsv(mk, "data/MK/metadata.tsv")



# SRA ---------------------------------------------------------------------

sra <- read_csv("data/MK/SraRunTable.txt")

# Remove nuclei
sra <- filter(sra, !str_detect(workflow, "Single-nuclei"))

# Remove bulk
sra <- filter(sra, !str_detect(workflow, "Bulk"))

# Use `Sample Name` as directory names and add file with the 8 SRR numbers to each
table(sra$`Sample Name`)
for(i in unique(sra$`GEO_Accession (exp)`)){
  dir.create(str_glue('data/MK/{i}'))
  filter(sra, `GEO_Accession (exp)`==i) %>%
    select(Run) %>%
    write_tsv(str_glue('data/MK/{i}/SRR_list.txt'), col_names = FALSE)
}

# Also write out file with sample names
data.frame(unique(sra$`Sample Name`)) %>%
  write_tsv("data/MK/sample_names.tsv", col_names = FALSE)



# Expression data ---------------------------------------------------------

# Download data
if(!dir.exists("data/MK/GSE141115_RAW")){
  download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141115&format=file",
              destfile = "data/MK/GSE141115_RAW")
  
  # Remove uneeded nuclei files
  files <- list.files("data/MK/GSE141115_RAW")
  files_to_remove <- files[!str_remove_all(str_split_fixed(files, "_",2)[,2], ".txt.gz")%in%unique(str_split_fixed(mk$cell_barcode,"_",2)[,1])]
  unlink(str_glue("data/MK/GSE141115_RAW/{files_to_remove}"))
  }

# Visualise a sample to make sure the cell type labels match
mk.example <- read.table("data/MK/GSE141115_RAW/GSM4195102_BG1.txt.gz") %>%
  CreateSeuratObject(project="GSM4195102_BG1") %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims=1:15) %>%
  FindNeighbors() %>%
  FindClusters()
mk.example$cell_type <- mk$cell_type[match(colnames(mk.example), mk$cell_barcode)]

ggsave(filename = "figures/MK_check.png",
       plot = DimPlot(mk.example, group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend(),
       device = "png",width = 25,height = 20,units = "cm")

# Yep, the cell type labels match
