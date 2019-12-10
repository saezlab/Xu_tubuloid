# Bioinformatics analysis on Tubuloid samples of single-cell RNA-seq
The following samples were profiled at the single-cell resolution:
* CK5 early organoid. Organoid generated from CD24+ sorted cells from human adult kidney tissue at an early stage.
* CK119 late organoid. Organoid generated from CD24+ sorted cells from human adult kidney tissue at a late stage.
* CK120 CD13. CD13+ sorted cells from human adult kidney tissue.
* CK120 CD24. CD24+ sorted cells from human adult kidney tissue.

## Material and methods
> Bioinformatics analysis of single-cell RNAseq. FastQ files were aligned to the human reference genome (GRCh38 assembly), 
and the expression was quantified using CellRanger (10x Genomics, version 3.1). Then the dataset was analyzed using the Seurat 
R package (version 3.0.2) [1]. Quality control diagnostics were performed by the examination of the total number of detected genes 
per cell, mitochondrial gene expression and library sizes for each sample. Cells with less than 200 or greater than 6,000 expressed 
genes were discarded to avoid bad quality cells and cell doublets, respectively. Taking into account that proximal tubule cells 
potentially express mitochondrial genes in high magnitude, a lenient threshold of greater than 80% mitochondrial expression was 
defined to drop low quality cells at first instance, but revisited after cell annotation to discard non-proximal tubule cells with 
high mitochondrial expression content. Raw gene expression was normalized by the library size of each cell, and multiplied by a scaling 
factor of 10,000 afterwards, as it is proposed in the Seurat workflow. A feature selection prior high dimensionality reduction was 
performed selecting the 2,000 most highly variable genes using the ‘vst’ method implemented in Seurat. Then, principal component analysis 
was applied to select the first 25 principal components for graph-based clustering using Shared-Nearest Neighbour algorithm for cell clustering, 
and Louvain to find clusters (resolution of 0.5). The cell type assignment among the resulting cell clusters was manually curated based on 
cluster-specific markers. For this, gene specificity scores and gene conditional probabilities for each cluster were obtained using 
genesorteR (version 0.3.1) [2]. Additionally, differentially over-expressed genes from each cluster were obtained using transcriptome-wide 
wilcoxon rank sum tests implemented in Seurat, but limiting the individual gene tests to only those genes that were expressed in 10% 
of the cells within the cluster with an average over-expression of 0.25 log2-fold-change. When it is stated in the results, PTPRC and 
PECAM1 gene expression was used to identify contaminant hematopoietic and endothelial cell populations in sorted cells, respectively. 
Similarly, the remaining unrelated cell clusters to the biological context under study with high content of ribosomal and mitochondrial 
genes were also discarded. Proliferating cell clusters were identified by scoring the average expression of cell cycle genes [3], as 
compared to a control set as proposed in the CellCycleScoring function from Seurat with default parameters. UMAP was used to reduce 
the 25 principal components into 2-dimensional space, showing the cell cluster labels from the cell assignment. Raw sequencing data 
and processed matrices from scRNA-seq were deposited in the Gene Expression Omnibus database with the accession number _GEOXXXXX_.

## Usage
You can reproduce the results following these instructions:
1. [Optionally] Download the raw 10x data (FastQ files). Then preprocess as stated in the material and methods section.
2. Download the data matrices from CellRanger, that are included in the supplementary files from the GEO entry.
3. Run the shuttle to build the html rmarkdowns reports

```bash
bash run_rmd.sh
```
## References
1. Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM, Hao Y, Stoeckius M, Smibert P, Satija R. Comprehensive Integration of Single-Cell Data. Cell. 2019 Jun 13;177(7):1888-1902.e21
2. MM Ibrahim and  Kramann R. genesorteR: Feature Ranking in Clustered Single Cell Data. bioRxiv
3. Tirosh I, Izar B, Prakadan SM, Wadsworth MH, Treacy D, Trombetta JJ, Rotem A, Rodman C, Lian C, Murphy G, Fallahi-Sichani M, Dutton-Regester K, Lin JR, Cohen O, Shah P, Lu D, Genshaft AS, Hughes TK, Ziegler CG, Kazer SW, Gaillard A, Kolb KE, Villani AC, Johannessen CM, Andreev AY, Van Allen EM, Bertagnolli M, Sorger PK, Sullivan RJ, Flaherty KT, Frederick DT, Jané-Valbuena J, Yoon CH, Rozenblatt-Rosen O, Shalek AK, Regev A, Garraway LA. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science. 2016 Apr 8;352(6282):189-96
