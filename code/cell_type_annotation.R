# Script information ------------------------------------------------------

# title: Annotate cell types
# author: Walter Muskovic
# date: 2020-12-15
# description: Quick classification of cell types in each dataset using 
# marker genes or scPred. This is done so that when we assess the QC 
# metrics later we can split the cells by type



# Import libraries --------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(biomaRt)
  library(Matrix)
  library(scPred)
})



# Seurat pre-processing ---------------------------------------------------

get_Seurat <- function(input_name, save_seurat = TRUE){
  # Load the data
  input.data <- Read10X(data.dir = str_glue('../data/{input_name}/outs/filtered_feature_bc_matrix/'))
  # Initialise the Seurat object with the raw (non-normalised data).
  input.seurat <- CreateSeuratObject(counts = input.data, project = input_name)
  # Add the percentage of reads that map to the mitochondrial genome
  if(input_name!="MB"){
    input.seurat[["percent.mt"]] <- PercentageFeatureSet(input.seurat,
                                                         pattern = "^MT-")
  } else {
    input.seurat[["percent.mt"]] <- PercentageFeatureSet(input.seurat,
                                                         pattern = "^mt-")
  }
  # Visualize QC metrics as a violin plot
  ggsave(file = str_glue('../figures/{input_name}_seurat_QC.pdf'),
         plot = VlnPlot(input.seurat,
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        ncol = 3),
         device = "pdf", width=30, height=20, units = "cm")
  # Filter out cells with MT gene content > 15%
  input.seurat <- subset(input.seurat, subset = percent.mt < 15)
  # Run sctransform
  input.seurat <- SCTransform(input.seurat, return.only.var.genes = FALSE)
  # Visualization and clustering
  input.seurat <- RunPCA(input.seurat)
  input.seurat <- RunUMAP(input.seurat, dims = 1:30)
  input.seurat <- FindNeighbors(input.seurat, dims = 1:30)
  input.seurat <- FindClusters(input.seurat)
  # Cell cycle scoring
  if(input_name!="MB"){
    input.seurat <- CellCycleScoring(input.seurat,
                                     s.features = cc.genes$s.genes,
                                     g2m.features = cc.genes$g2m.genes,
                                     set.ident = FALSE)
  } else {
    # Following code snippet borrowed from user "radiaj" post on 
    #https://rjbioinformatics.com/ "Converting mouse to human gene names with 
    #biomaRt package"
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    s.genes.mouse = getLDS(attributes = c("mgi_symbol"),
                           filters = "mgi_symbol",
                           values = cc.genes$s.genes,
                           mart = mouse,
                           attributesL = c("hgnc_symbol"),
                           martL = human,
                           uniqueRows=T)
    s.genes.mouse <- unique(s.genes.mouse[, 1])
    g2m.genes.mouse = getLDS(attributes = c("mgi_symbol"),
                             filters = "mgi_symbol",
                             values = cc.genes$g2m.genes,
                             mart = mouse,
                             attributesL = c("hgnc_symbol"),
                             martL = human,
                             uniqueRows=T)
    g2m.genes.mouse <- unique(g2m.genes.mouse[, 1])
    input.seurat <- CellCycleScoring(input.seurat,
                                     s.features = s.genes.mouse,
                                     g2m.features = g2m.genes.mouse,
                                     set.ident = FALSE)
  }
  # UMAP plot with identified clusters
  ggsave(file = str_glue('../figures/{input_name}_UMAP_clusters.pdf'),
         plot = DimPlot(input.seurat, label = TRUE) + NoLegend() +
           FeaturePlot(input.seurat,
                       c("nCount_RNA", "S.Score","G2M.Score","percent.mt"),
                       min.cutoff = "q5", max.cutoff = "q99"),
         device = "pdf", width=40, height=15, units = "cm")
  
  # Save Seurat object if requested
  if(save_seurat){
    saveRDS(input.seurat, str_glue('../data/{input_name}/seurat_SCT.rds'))
  }
  
  return(input.seurat)
}

gbm <- get_Seurat("GBM")
pbmc <- get_Seurat("PBMC")
hl <- get_Seurat("HL")
mb <- get_Seurat("MB")
gbm <- get_Seurat("GBM5p")



# GBM marker genes --------------------------------------------------------

gbm <- readRDS(str_glue('../data/GBM/seurat_SCT.rds'))

# Define marker genes

# reference 1 : Dusart et al, Cell Reports, 2019 - "A Systems-Based Map of Human
# Brain Cell-Type Enriched Genes and Malignancy-Associated Endothelial Changes",
# https://doi.org/10.1016/j.celrep.2019.09.088)

# reference 2 : Neftel et al, Cell, 2019 - "An Integrative Model of Cellular
# States, Plasticity, and Genetics for Glioblastoma",
# https://doi.org/10.1016/j.cell.2019.06.024)

# reference 3 : Couturier et al, Nature Communications, 2020 - "Single-cell
# RNA-seq reveals that glioblastoma recapitulates a normal neurodevelopmental 
# hierarchy", https://doi.org/10.1038/s41467-020-17186-5)

# reference 4 : Wang et al, Cancer Discovery, 2019 - "The Phenotypes of 
# Proliferating Glioblastoma Cells Reside on a Single Axis of Variation", 
# https://doi.org/10.1158/2159-8290.cd-19-0329)

oligodendrocyte <- c("MOG", "MAG", "CNP") # ref1
oligodendrocyte <- c(oligodendrocyte, "MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11") %>% unique() # ref2
microglia_macrophage <- c("C1QA", "AIF1", "LAPTM5") # ref1
microglia_macrophage <- c(microglia_macrophage, "CD14", "AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R") %>% unique() # ref2
t_cell <- c("CD2", "CD3D", "CD3E", "CD3G") # ref2
endothelial <- c("CD34", "CLEC14A", "VWF") # ref1
endothelial <- c(endothelial, "ESAM") # ref3
endothelial <- c(endothelial, "APOLD1", "CLDN5") # ref4

# Define and apply function to add scores for cell types using the defined
# marker genes
score_marker_genes <- function(seurat_object, nbins=24){
  AddModuleScore(seurat_object, features = list(oligodendrocyte, microglia_macrophage, t_cell, endothelial),
                 name=c("oligodendrocyte", "microglia_macrophage", "t_cell", "endothelial"), nbin=nbins)
}

gbm <- score_marker_genes(gbm)

# Define and apply function to plot scores
plot_marker_genes <- function(seurat_object, image_name){
  ggsave(file = str_glue('../figures/{image_name}.pdf'),
         plot = FeaturePlot(seurat_object,
                            features = c("oligodendrocyte1", "microglia_macrophage2", "t_cell3", "endothelial4"),
                            min.cutoff = "q10", max.cutoff = "q90",
                            ncol=2, label=TRUE, order = TRUE),
         device = "pdf", width = 40, height = 30, units = "cm")
}
plot_marker_genes(gbm, "gbm_markers")

# Add annotations to the object, plot and save out again
gbm$cell_type <- "malignant"
gbm$cell_type[gbm$seurat_clusters%in%c(4,5)] <- "oligodendrocyte"
gbm$cell_type[gbm$seurat_clusters%in%c(1,7,6,12)] <- "myeloid"
gbm$cell_type[gbm$seurat_clusters==10] <- "t_cell"
gbm$cell_type[gbm$seurat_clusters==14] <- "endothelial"
ggsave(file = str_glue('../figures/gbm_cell_types.pdf'),
       plot = DimPlot(gbm, group.by = "cell_type", label=TRUE) + NoLegend(),
       device = "pdf", width = 20, height = 15, units = "cm")
saveRDS(gbm, str_glue('../data/GBM/seurat_SCT_anno.rds'))



# HL marker genes ---------------------------------------------------------

hl <- readRDS(str_glue('../data/HL/seurat_SCT.rds'))

# reference 5 :  Aoki et al, Cancer Discovery, 2020 - "Single-Cell Transcriptome
# Analysis Reveals Disease-Defining T-cell Subsets in the Tumor Microenvironment
# of Classic Hodgkin Lymphoma", https://doi.org/10.1158/2159-8290.cd-19-0680)

# reference 6 :  Schafflick et al, Nature Communications, 2020 - "Integrated 
# single cell analysis of blood and cerebrospinal fluid leukocytes in multiple 
# sclerosis", https://doi.org/10.1038/s41467-019-14118-w)

# reference 7 :  Ohneda et al, Int J Mol Sci., 2019 - "Mouse Tryptase Gene 
# Expression is Coordinately Regulated by GATA1 and GATA2 in Bone Marrow-Derived
# Mast Cells", https://dx.doi.org/10.3390%2Fijms20184603)

b_cell <- c("MS4A1") # ref5
macrophage <- c("CD68", "IDO1") # ref5
plasmacytoid_dendritic_cell <- c("CLEC4C","NRP1") # ref5
erythrocyte <- c("HBB", "HBA1", "HBA2") # ref6
cytotoxic_t_cell <- c("GZMA", "GZMK", "IFNG") # ref5
regulatory_t_cell <- c("FOXP3", "IL2RA", "IKZF2") # ref5
t_helper <- c("CXCL13","PDCD1","FABP5") # ref5
naive_t_cell <- c("CCR7", "IL7R", "LEF1") # ref5
progenitor <- c("CD34") # ref5
mast_cell <- c("TPSAB1","TPSB2", "KIT", "GATA1", "GATA2") # ref7

score_marker_genes <- function(seurat_object, nbins=24){
  AddModuleScore(seurat_object, 
                 features = list(b_cell, macrophage, plasmacytoid_dendritic_cell, erythrocyte, cytotoxic_t_cell, regulatory_t_cell, t_helper, naive_t_cell, progenitor, mast_cell),
                 name=c("b_cell", "macrophage", "plasmacytoid_dendritic_cell", "erythrocyte", "cytotoxic_t_cell", "regulatory_t_cell", "t_helper", "naive_t_cell", "progenitor", "mast_cell"),
                 nbin=nbins)
}
hl <- score_marker_genes(hl)

plot_marker_genes <- function(seurat_object, image_name){
  ggsave(file = str_glue('../figures/{image_name}.pdf'),
         plot = FeaturePlot(seurat_object,
                            features = c("b_cell1", "macrophage2", "plasmacytoid_dendritic_cell3", "erythrocyte4", "cytotoxic_t_cell5", "regulatory_t_cell6", "t_helper7", "naive_t_cell8", "progenitor9", "mast_cell10"),
                            min.cutoff = "q10", max.cutoff = "q90",
                            ncol=4, label=TRUE, order = TRUE),
         device = "pdf", width = 50, height = 40, units = "cm")
}
plot_marker_genes(hl, "hl_markers")

hl$cell_type <- "NA"
hl$cell_type[hl$seurat_clusters%in%c(8,2)] <- "b_cell"
hl$cell_type[hl$seurat_clusters==4] <- "macrophage"
hl$cell_type[hl$seurat_clusters%in%c(9,10)] <- "malignant"
hl$cell_type[hl$seurat_clusters==11] <- "plasmacytoid_dendritic_cell"
hl$cell_type[hl$seurat_clusters==5] <- "cytotoxic_t_cell"
hl$cell_type[hl$seurat_clusters%in%c(0,7,6)] <- "regulatory_t_cell"
hl$cell_type[hl$seurat_clusters==1] <- "t_helper"
hl$cell_type[hl$seurat_clusters==3] <- "naive_t_cell"
hl$cell_type[hl$seurat_clusters==12] <- "progenitor"
hl$cell_type[hl$mast_cell10>2] <- "mast_cell"

ggsave(file = str_glue('../figures/hl_cell_types.pdf'),
       plot = DimPlot(hl, group.by = "cell_type", label=TRUE) + NoLegend(),
       device = "pdf", width = 20, height = 15, units = "cm")
saveRDS(hl, str_glue('../data/HL/seurat_SCT_anno.rds'))



# PBMC scPred annotation --------------------------------------------------

# Get reference and query 
query <- readRDS(str_glue('../data/PBMC/seurat_SCT.rds'))
reference <- merge(scPred::pbmc_1, scPred::pbmc_2) %>%
  SCTransform() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

# Find more finely resolved clusters
query <- FindClusters(query, resolution = 1.2)

# Get feature space
reference <- getFeatureSpace(reference, "cell_type")
# Train classifier
reference <- trainModel(reference)
# Classify
query <- scPredict(query, reference)

# Look at how clusters and annotations compare
as_tibble(query[[]]) %>% 
  dplyr::select(scpred_prediction, seurat_clusters) %>%
  table

# Assign cell types
query$cell_type <- "NA"
query$cell_type[query$seurat_clusters%in%c(12,13,15,17,22)] <- "B.cell"
query$cell_type[query$seurat_clusters%in%c(1,5,6,14)] <- "CD4.T.cell"
query$cell_type[query$seurat_clusters%in%c(0,8,9,19)] <- "CD8.T.cell"
query$cell_type[query$seurat_clusters==16] <- "cDC"
query$cell_type[query$seurat_clusters%in%c(2,3,4,10)] <- "cMono" 
query$cell_type[query$seurat_clusters%in%c(7)] <- "ncMono"
query$cell_type[query$seurat_clusters==11] <- "NK.cell"
query$cell_type[query$seurat_clusters==18] <- "pDC"
query$cell_type[query$seurat_clusters==21] <- "Plasma.cell"
query$cell_type[query$seurat_clusters==20] <- "unresolved"

ggsave(file = str_glue('../figures/pbmc_cell_types.pdf'),
       plot = DimPlot(query, group.by = "cell_type", label=TRUE) + NoLegend(),
       device = "pdf", width = 20, height = 15, units = "cm")
saveRDS(query, str_glue('../data/PBMC/seurat_SCT_anno.rds'))



# Mouse brain scPred annotation -------------------------------------------

mb <- readRDS(str_glue('../data/MB/seurat_SCT.rds'))

# Find more finely resolved clusters
mb <- FindClusters(mb, resolution = 1.2)

# Get trained scPred models
ref1 <- readRDS("../data/reference/Allen_Brain_mouse_scPred_model.rds")
ref2 <- readRDS("../data/reference/developing_mouse_brain_scPred_model.rds")

# Classify
query <- scPredict(mb, ref1)
ggsave(file = str_glue('../figures/mb_allen_classification.pdf'),
       plot = DimPlot(query,
                      group.by = "scpred_prediction", label=TRUE, repel = TRUE),
       device = "pdf", width = 20, height = 15, units = "cm")

# Look at how clusters and annotations compare
as_tibble(query[[]]) %>% 
  dplyr::select(scpred_prediction, seurat_clusters) %>%
  table

# Most of the cells are unassigned by scPred. 
# We can however discover that:
# Seurat cluster 12, 14 & 25 = astrocyte
# Seurat cluster 20 = endothelial cell
# Seurat cluster 21 = cajal-retzius
# Seurat cluster 22 = smooth muscle cell/pericyte
# Seurat cluster 24 = microglia/perivascular macrophage
# Seurat cluster 26 = oligodendrocyte

# We will likely have more luck with the training data produced from the 
# developing mouse brain reference from Rosenberg et al, Science, 2018
query <- scPredict(mb, ref2)
ggsave(file = str_glue('../figures/mb_developing_brain_classification.pdf'),
       plot = DimPlot(query,
                      group.by = "scpred_prediction", label=TRUE, repel = TRUE),
       device = "pdf", width = 35, height = 25, units = "cm")

# Look at how clusters and annotations compare
as_tibble(query[[]]) %>% 
  dplyr::select(scpred_prediction, seurat_clusters) %>%
  table

# From this classification we discover:
# Seurat cluster 0 & 19 = Neuron_hippocampus
# Seurat cluster 1, 7 & 11 = Neuron_cortex
# Seurat cluster 2 = Neuron_unresolved_1
# Seurat clusters 3, 13 & 17 = Migrating_interneuron
# Seurat cluster 4 = Neuron_unresolved_2
# Seurat clusters 5, 8 & 19 = Neuron_hippocampus
# Seurat cluster 6 = Neuron_unresolved_3
# Seurat cluster 9 = Neuron_unresolved_4
# Seurat cluster 10 = Astrocyte_unresolved_1
# Seurat cluster 12 = Astrocyte
# Seurat cluster 14 = Astrocyte_unresolved_2
# Seurat cluster 15 = Neuron_unresolved_5
# Seurat cluster 16 = Neuron_cerebellum
# Seurat cluster 18 = erythrocyte
# Seurat cluster 20 = endothelial cell
# Seurat cluster 21 = cajal-retzius
# Seurat cluster 22 = smooth muscle cell/pericyte
# Seurat cluster 23 = Neuron_unresolved_6
# Seurat cluster 24 = microglia/perivascular macrophage
# Seurat cluster 25 = ependymal cell
# Seurat cluster 26 = oligodendrocyte/OPC

# Assign cell types
mb$cell_type <- NA
mb$cell_type[mb$seurat_clusters%in%c(0,5,8,19)] <- "neuron_hippocampus"
mb$cell_type[mb$seurat_clusters%in%c(1,7,11)] <- "neuron_cortex"
mb$cell_type[mb$seurat_clusters%in%c(2)] <- "neuron_unresolved_1"
mb$cell_type[mb$seurat_clusters%in%c(3,13,17)] <- "migrating_interneuron"
mb$cell_type[mb$seurat_clusters%in%c(4)] <- "neuron_unresolved_2"
mb$cell_type[mb$seurat_clusters%in%c(6)] <- "neuron_unresolved_3"
mb$cell_type[mb$seurat_clusters%in%c(9)] <- "neuron_unresolved_4"
mb$cell_type[mb$seurat_clusters%in%c(10)] <- "astrocyte_unresolved_1"
mb$cell_type[mb$seurat_clusters%in%c(12)] <- "astrocyte"
mb$cell_type[mb$seurat_clusters%in%c(14)] <- "astrocyte_unresolved_2"
mb$cell_type[mb$seurat_clusters%in%c(15)] <- "neuron_unresolved_5"
mb$cell_type[mb$seurat_clusters%in%c(16)] <- "neuron_cerebellum"
mb$cell_type[mb$seurat_clusters%in%c(18)] <- "erythrocyte"
mb$cell_type[mb$seurat_clusters%in%c(20)] <- "endothelial_cell"
mb$cell_type[mb$seurat_clusters%in%c(21)] <- "cajal-retzius"
mb$cell_type[mb$seurat_clusters%in%c(22)] <- "smooth_muscle_cell/pericyte"
mb$cell_type[mb$seurat_clusters%in%c(23)] <- "neuron_unresolved_6"
mb$cell_type[mb$seurat_clusters%in%c(24)] <- "microglia/perivascular_macrophage"
mb$cell_type[mb$seurat_clusters%in%c(25)] <- "ependymal_cell"
mb$cell_type[mb$seurat_clusters%in%c(26)] <- "oligodendrocyte/OPC"

ggsave(file = str_glue('../figures/mb_cell_types.pdf'), 
       plot = DimPlot(mb, group.by = "cell_type"),
       device = "pdf", width = 25, height = 15, units = "cm")
saveRDS(mb, str_glue('../data/MB/seurat_SCT_anno.rds'))



# GBM5p marker genes ------------------------------------------------------

gbm5p <- readRDS(str_glue('../data/GBM5p/seurat_SCT.rds'))

# Apply function to add scores for cell types
gbm5p <- score_marker_genes(gbm5p)

# Apply function to plot scores
plot_marker_genes(gbm5p, "gbm5p_markers")

# Add annotations to the object, plot and save out again
gbm5p$cell_type <- "malignant"
gbm5p$cell_type[gbm5p$seurat_clusters%in%c(3,7)] <- "oligodendrocyte"
gbm5p$cell_type[gbm5p$seurat_clusters%in%c(2,4,8,9,12)] <- "myeloid"
gbm5p$cell_type[gbm5p$seurat_clusters==10] <- "t_cell"
gbm5p$cell_type[gbm5p$seurat_clusters==13] <- "endothelial"

ggsave(file = str_glue('../figures/gbm5p_cell_types.pdf'),
       plot = DimPlot(gbm5p, group.by = "cell_type", label=TRUE) + NoLegend(),
       device = "pdf", width = 20, height = 15, units = "cm")
saveRDS(gbm5p, str_glue('../data/GBM5p/seurat_SCT_anno.rds'))
