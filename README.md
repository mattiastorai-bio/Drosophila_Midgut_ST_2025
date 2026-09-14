# Redhai, Storai *et al.*, 2026

This repository contains scripts and files supporting the Master's Thesis: </br>

Storai M., **A Spatial Transcriptomics Atlas of Drosophila Midgut** (2026). Ruprecht-Karls-Universität Heidelberg, supervised by Prof. Michael Boutros (DKFZ) and Dr Siamak Redhai (DKFZ).

## Summary

The *Drosophila melanogaster* midgut is an established model for studying stem cell maintenance and tissue homeostasis. Recently, single cell RNA sequencing has provided extended characterization of midgut cell types, however, this technique inherently loses spatial information. In this thesis, using the 10x Genomics Visium HD platform, we generated the first spatial transcriptomics atlas of the adult *Drosophila* midgut. This resource was integrated with established scRNA-seq atlases to describe the spatial gene expression profile and the distribution of intestinal cell types at high-resolution. My findings have led to the identification of over 100 new cell-type-specific and region-specific markers, as well as the characterization of two previously undetected cell types within single cell datasets. Moreover, characterising gene programs across spatial domains revealed that ISCs are transcriptionally similar with region-specific gene expression emerging in immediate daughter cells and fully established in differentiated cells. Using hierarchical clustering and transcriptome comparison tools we validated previous boundary regions and further identified a new boundary. To understand gene expression programs in three dimensions, we utilized micro-CT scans of the gut and integrated spatial gene expression, revealing topologically associated gene programs. To better characterize the mechanisms involved in region-specific gene expression and cell type formation, I combined temperature-inducible CRISPR knock-out of the *Apc1* tumour suppressor in progenitor cells with spatial transcriptomics. *Apc1* mutant midguts displayed localized disruptions in cell type formation and an anteriorized shift in gene expression profiles, suggesting that Wnt signalling is involved in the spatial pattern of gene expression and cell type formation. Overall, this thesis builds a spatially informed resource for the *Drosophila* midgut and uses it to refine and expand previous knowledge, linking cell identity and tissue architecture, and offering new biological insights into regionalized features, cell populations and intestinal function.

## Repository Structure

The code is divided into two primary sections matching the experimental conditions:

```
Drosophila\_Midgut\_ST\_2025/
├── WT/              # All scripts \& analysis for the Wild Type dataset (111125\_ST)
└── Mutant/          # All scripts \& comparative analysis for the Apc KO mutant dataset (03slide4)
```

### WT/

Contains the core spatial transcriptomics processing pipeline, integration with single-cell datasets, trajectory inference, topological association analyses, and 3D visualizations developed for the wild-type condition.

### Mutant/

Contains the modified pipeline scripts for processing the mutant slide (03slide4) along with novel comparative analyses:

* `Gene\_Spatial\_Comparison\_WT\_vs\_Mutant.Rmd`: Direct comparison of target gene spatial profiles
* `Global\_Difference\_Comparison.Rmd`: Genome-wide transcriptional divergence between WT and mutant
* `Hierarchical\_Clustering\_Of\_Regions.Rmd`: Gap statistic boundary detection in the mutant
* `concordance.Rmd`: Gene expression concordance analysis

## Requirements

The analyses were performed in R. Key packages include:

* **Spatial transcriptomics**: Seurat, SpatialExperiment, spacexr (RCTD), ggspavis
* **Single-cell analysis**: Seurat, anndataR, SeuratDisk
* **Trajectory inference**: monocle3, slingshot
* **Deconvolution**: spacexr, SPOTlight, DOTr
* **Functional enrichment**: clusterProfiler, enrichplot, ReactomePA, org.Dm.eg.db
* **Visualization**: ggplot2, plotly, pheatmap, ComplexHeatmap, bslib (Shiny)
* **Data handling**: dplyr, data.table, Matrix, hdf5r, arrow, rhdf5

See [`packages.csv`](packages.csv) for a complete list of packages used across all scripts.

