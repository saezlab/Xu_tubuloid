# Bioinformatics analysis of the transcriptome from tubuloid-derived kidney organoids

> The Kramann laboratory has established a protocol to develop tubuloid kidney organoids derived from adult kidney cells.
These organoids recapitulate better the characteristics of adult tissue. These can be also used for disease modelling by
gene editing (e.g. polycystic kidney disease - PKD). Here the transcriptome of several normal and PKD organoids are
characterized at a gingle-cell resolution. Together with ADPKD and control kidney tissue.

## Index 
You might take a look at the `markdown` reports (`.md` files) included within this repository. 

* Individual analysis of sorted kidney cells: [sorted CD13+ cells](Individual_analysis_CK120_CD13), [sorted CD24+ cells](Individual_analysis_CK121_CD24)
* Individual analysis of several kidney organoids: [early control organoid](Individual_analysis_CK5_early_organoid), 
[late control organoid](Individual_analysis_CK119_late_organoid), [Clevers' control organoid](Individual_analysis_JX1_HC_organoid), [PKD1-KO PKD organoid](Individual_analysis_JX2_PKD1KO_organoid),
[PKD2-KO PKD organoid](Individual_analysis_JX3_PKD2KO_organoid).
* Individual analysis of 3x ADPKD patient kidney tissue: [ADPKD3](Individual_analysis_ADPKD3), [CK224](Individual_analysis_CK224_PT25_PKD2-), [CK225](Individual_analysis_CK225_PT8_PKD1-).
* Individual analysis of 2x healthy patient kidney tissue: [Control1](Individual_analysis_Control1), [Control2](Individual_analysis_Control2)
* Integrated analysis for differential expression in ADPKD tissue: [Cystic signature](Cystic_signature)
* [Merged analysis](Merged_analysis)  to make figures for the manuscript: UMAP plots with standarized color scheme and marker expression on normalized data.

> These markdowns consists on individual, integrated and merged analyses. 
The individual analyses aim to processed the data and annotate distinct clusters in each sample.
The integrated analysis harmonize several samples for a statistical test for differential gene expression in ADPKD condition (tissue).
The merged analysis simply combine processed data from several samples to find patterns in gene expression across them.


## Data description

Samples were profiled by single-cell (organoid and sorted cells) or single-nuclei RNA-seq 10x technology.

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
Processed UMI gene expression values from the single-cell RNA-seq are deposited at https://doi.org/10.6084/m9.figshare.11786238 - and will be published (with DOI) upon manuscript acceptance.


## Reproducibility
Reproducible results are available following these steps:

1. [Optionally] Request access to the raw sequencing 10x data (FastQ files) to Rafael Kramann following the manuscript details. Then preprocess the data following your preferences.
2. Download the data matrices from CellRanger. The processed data with gene expression values are available at https://doi.org/10.6084/m9.figshare.11786238. These data must be located in `./data/sc/`.
3. Run the make script to re-build the html rmarkdowns reports

```bash
bash make.sh
```

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


