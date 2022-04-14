# Bioinformatics analysis of the transcriptome from tubuloid kidney organoids

> The Kramann laboratory has established a protocol to develop tubuloid kidney organoids derived from adult kidney cells.
These organoids recapitulate better the characteristics of adult tissue. These can be also used for disease modelling by
gene editing (e.g. polycystic kidney disease - PKD). Here the transcriptome of several normal and PKD organoids are
characterized at a single-cell resolution. Together with ADPKD and control kidney tissue.

Code repository for the manuscript: Xu\*, Kuppe\* et al. 2022.   
For full description about samples, methods and data availability, please refer to the corresponding sections from the manuscript.

## Structure 

The repository is divided in folders, each represents a set of scripts for a given analysis. Scripts are enumerated in increasing order following their rational order within the analysis,
with a descriptive filename and title. Certain scripts in a subfolder might present dependencies over the output of other analysis. Thus these [analyses](#Analyses) are presented in order.

You might take a look at the `markdown` reports (`.md` files) and jupyter notebooks (`ipynb` files) for a quick understanding of data processing/modeling, inputs and outputs. 
In general, each code file has an output folder with figures and files. Only embedded figures are included in the code repository to build the reports.
Please refer to the [data availability](#data-availability) for intermediate output files from the analyses.

## Analyses (in order)

### 1. Data ([`data/`](data))

It contains the initial input data (processed data) and external data (prior knowledge manually curated) for the analysis. The [data availability](#Data-availability) is described below.

### 2. Individual analyses (`Individual_analysis_*/`)
These are several independent analyses on each sample (sequencing library) to characterize their heterogeneity in an unsupervised manner 
(QC filtering, unsupervised clustering, and reiteration). Doing so, these individual analyses aim to processed the data and annotate distinct clusters in each sample.

* Individual analysis of sorted kidney cells: [sorted CD13+ cells](Individual_analysis_CK120_CD13), [sorted CD24+ cells](Individual_analysis_CK121_CD24)
* Individual analysis of several kidney organoids: [early control organoid](Individual_analysis_CK5_early_organoid), 
[late control organoid](Individual_analysis_CK119_late_organoid), [Clevers' control organoid](Individual_analysis_JX1_HC_organoid), [PKD1-KO PKD organoid](Individual_analysis_JX2_PKD1KO_organoid),
[PKD2-KO PKD organoid](Individual_analysis_JX3_PKD2KO_organoid).
* Individual analysis of 3x ADPKD patient kidney tissue: [ADPKD3](Individual_analysis_ADPKD3), [CK224](Individual_analysis_CK224_PT25_PKD2-), [CK225](Individual_analysis_CK225_PT8_PKD1-).
* Individual analysis of 2x healthy patient kidney tissue: [Control1](Individual_analysis_Control1), [Control2](Individual_analysis_Control2)

### 3. Integrated atlas of human kidney tissue for differential expression analysis in ADPKD ([`Cystic_signature`](Cystic_signature))

The integrated analysis harmonize several samples for a statistical test for differential gene expression in ADPKD condition (tissue) as compared to control: [Cystic signature](Cystic_signature)

### 4. Merged analysis to make figures ([`Merged_analysis/`](Merged_analysis))

The merged analysis simply combine output data from several samples to find patterns in gene expression across them.
These include: UMAP plots with standarized color scheme, cell proportions and marker expression on normalized data.

### 5. Tubuloid cell mapping into human kidney tissue ([`RefMap/`](RefMap))

[Cell type annotation transfer](RefMap/01_symphony_mapping.R) for tubuloid cells via Symphony compression from the [Integrated atlas](Cystic_signature/01_Harmony_integration.md) of human kidney tissue, 
using the [curated cell type annotations](https://github.com/saezlab/Xu_tubuloid/blob/master/Cystic_signature/04_clusters.md) from the previous analysis. 

### 6. Cell-to-cell communication analysis ([`Cell-cell_communication`](Cell-cell_communication))

Cell-to-cell communication analysis via LIANA (CellPhoneDB) for each kidney condition, and differential communication using CrossTalkeR.


## Samples

Samples were profiled by single-cell (organoid and sorted cells) or single-nuclei (human tissue samples) RNA-seq 10x technology.

| Sample | Type | Condition | Description |
| :--- | :--- | :--- | :--- |
| CK5_early_organoid | Organoid | Healthy/Control | Organoid generated from CD24+ sorted cells from human adult kidney tissue at an early stage |
| CK119_late_organoid | Organoid | Healthy/Control | Organoid generated from CD24+ sorted cells from human adult kidney tissue at a late stage |
| JX1_HC_organoid | Organoid | Healthy/Control | Organoid generated following Hans Clever's protocol for kidney organoids at a late stage development |
| JX2_PKD1KO_organoid | Organoid | PKD/disease | Organoid generated from CD24+ sorted cells from human adult kidney tissue, for which PKD1 was gene-edited to reproduce PKD phenotype, developed at a late stage |
| JX3_PKD2KO_organoid | Organoid | PKD/disease | Organoid generated from CD24+ sorted cells from human adult kidney tissue, for which PKD2 was gene-edited to reproduce PKD phenotype, developed at a late stage |
| CK120_CD13 | Sorted cells | Healthy/Control | CD13+ sorted cells from human adult kidney tissue |
| CK121_CD24 | Sorted cells | Healthy/Control | CD24+ sorted cells from human adult kidney tissue |
| CK224_PT25_PKD2- | Kidney tissue | ADPKD/disease | human specimen with ADPKD (PKD2- genotype) | 
| CK225_PT8_PKD1- | Kidney tissue | ADPKD/disease | human specimen with ADPKD (PKD1- genotype) | 
| ADPKD3 | Kidney tissue | ADPKD/disease | human specimen with ADPKD (Not described genotype) | 
| Control1 | Kidney tissue | Control/Healthy | human specimen with healthy kidney function | 
| Control2 | Kidney tissue | Control/Healthy | human specimen with healthy kidney function | 

## Data availability
Processed UMI gene expression values from the single-cell RNA-seq and intermediate output files from the analyses are deposited at https://doi.org/10.6084/m9.figshare.11786238 

## Environment

**Main required packages**

| Package name | source | version |
| :--- | :---: | ---: |
| Seurat | https://satijalab.org/seurat/install.html | 3.1.0 |
| ggplot2 | CRAN / https://ggplot2.tidyverse.org | 3.2.1 |
| clustree | Bioconductor | 0.4.1 |
| ComplexHeatmap | github (@jokergoo/ComplexHeatmap) | 2.0.0 |
| cowplot | CRAN | 1.0.0 |
| genesorteR | github (@mahmoudibrahim/genesorteR) | 0.3.1 |
| scanpy | https://scanpy.readthedocs.io/en/stable/ | 1.8.1 |


