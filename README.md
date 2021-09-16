# DropletQC manuscript code

Droplet-based single cell RNA-sequencing methods utilise microfluidics to encapsulate individual cells in water-in-oil droplets, dramatically increasing throughput compared to plate-based protocols. While encapsulating cells, droplets also capture cell-free "ambient" RNA, a complex mixture of transcripts released from damaged, stressed or dying cells. High concentrations of ambient RNA hinder identification of cell-containing droplets and degrade downstream biological interpretation, by violating the assumption that a droplet contains RNA from a single cell. The concentration of ambient RNA depends on the nature of the input cell suspension, being more abundant in freshly dissociated solid tissue and samples containing fragile cell types such as brain tissue.

We developed the [dropletQC R package](https://powellgenomicslab.github.io/DropletQC/ "DropletQC R package") which utilises a novel QC metric, the nuclear fraction, to  identify droplets containing ambient RNA or damaged cells. The nuclear fraction quantifies the amount of RNA in a droplet that originated from unspliced pre-mRNA. Ambient RNA consists mostly of mature cytoplasmic mRNA and is relatively depleted of unspliced nuclear precursor mRNA. Droplets containing only ambient RNA will have a low nuclear fraction score while damaged cells, due to loss of cytoplasmic RNA, will tend to have a higher score than intact cells.

The schematic below illustrates how the nuclear fraction score, in combination with total UMI counts, can be used to identify empty droplets and damaged cells.

![Figure_1](https://github.com/powellgenomicslab/DropletQC_paper/raw/main/figures/Figure_1.png)

To demonstrate the method on real-world data, we  applied it to four single cell RNA-seq datasets produced and made publicly available by 10x Genomics:

1. Mouse E18 Combined Cortex, Hippocampus and Subventricular Zone Cells, ~10k cells [dataset link](https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/SC3_v3_NextGem_DI_Neuron_10K "dataset link")
2. Hodgkin's Lymphoma, Dissociated Tumour, ~5k cells [dataset link](https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/Parent_NGSC3_DI_HodgkinsLymphoma "dataset link")
3. Human Glioblastoma Multiforme, ~5k cells [dataset link](https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/Parent_SC3v3_Human_Glioblastoma "dataset link")
4. PBMCs from a Healthy Donor, ~10k cells [dataset link](https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/Parent_NGSC3_DI_PBMC "dataset link")

The figure below illustrates how the method can be applied to identify both empty droplets and damaged cells.

![Figure_2](https://github.com/powellgenomicslab/DropletQC_paper/raw/main/figures/Figure_2.png)

This repository contains the code used to produced all of the analysis and figures published in:
[DropletQC: improved identification of empty droplets and damaged cells in single-cell RNA-seq data](https://doi.org/10.1101/2021.08.02.454717)
