# Bioinformatics analysis on Tubuloid-derived single-cell RNA-seq samples
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
the 25 principal components into 2-dimensional space, showing the cell cluster labels from the cell assignment.

> *References*:
> 1. Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM, Hao Y, Stoeckius M, Smibert P, Satija R. Comprehensive Integration of Single-Cell Data. Cell. 2019 Jun 13;177(7):1888-1902.e21
> 2. MM Ibrahim and  Kramann R. genesorteR: Feature Ranking in Clustered Single Cell Data. bioRxiv
> 3. Tirosh I, Izar B, Prakadan SM, Wadsworth MH, Treacy D, Trombetta JJ, Rotem A, Rodman C, Lian C, Murphy G, Fallahi-Sichani M, Dutton-Regester K, Lin JR, Cohen O, Shah P, Lu D, Genshaft AS, Hughes TK, Ziegler CG, Kazer SW, Gaillard A, Kolb KE, Villani AC, Johannessen CM, Andreev AY, Van Allen EM, Bertagnolli M, Sorger PK, Sullivan RJ, Flaherty KT, Frederick DT, Jané-Valbuena J, Yoon CH, Rozenblatt-Rosen O, Shalek AK, Regev A, Garraway LA. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science. 2016 Apr 8;352(6282):189-96

## Usage
You might take a look at the `markdown` reports (`.md` files) included within this repository. This consists on the individual single-cell analysis of four samples, each one goes through a bioinformatics analysis to identify cell type populations from the kidney in the samples. Finaly, a merged analysis with the four samples is performed to create the final figures for the manuscript.

* [Individual analysis of sorted CD13+ cells](Individual_analysis_CK120_CD13)
* [Individual analysis of sorted CD24+ cells](Individual_analysis_CK121_CD24)
* [Individual analysis of early organoid](Individual_analysis_CK5_early_organoid)
* [Individual analysis of late organoid](Individual_analysis_CK119_late_organoid)
* [Merged analysis to make figures for the manuscript](Merged_analysis)


Alternatively, you can reproduce the results following these instructions. Please note that this is going to over-write the markdowns from the repository:

1. [Optionally] Request access to the raw sequencing 10x data (FastQ files) to Rafael Kramann following the manuscript details. Then preprocess the data following your preferences.
2. Download the data matrices from CellRanger. The processed data with gene expression values are available at https://doi.org/10.6084/m9.figshare.11786238. These data must be located in `./data/sc/`.
3. Run the make script to re-build the html rmarkdowns reports

```bash
bash make.sh
```

## Environment

**Required packages**

| Package name | source | version |
| :--- | :---: | ---: |
| Seurat | https://satijalab.org/seurat/install.html | 3.1.0 |
| ggplot2 | CRAN / https://ggplot2.tidyverse.org | 3.2.1 |
| clustree | Bioconductor | 0.4.1 |
| ComplexHeatmap | github (@jokergoo/ComplexHeatmap) | 2.0.0 |
| cowplot | CRAN | 1.0.0 |
| genesorteR | github (@mahmoudibrahim/genesorteR) | 0.3.1 |


**Detailed environment for reproducibility**

Herein it is described the `sessionInfo()` from `R`:
```
R version 3.6.1 (2019-07-05)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] cowplot_1.0.0        clustree_0.4.1       ggraph_2.0.0.9000   
 [4] ggplot2_3.2.1        ComplexHeatmap_2.0.0 genesorteR_0.3.1    
 [7] Matrix_1.2-17        dplyr_0.8.3          GSEABase_1.46.0     
[10] graph_1.62.0         annotate_1.62.0      XML_3.98-1.20       
[13] AnnotationDbi_1.46.1 IRanges_2.18.2       S4Vectors_0.22.1    
[16] Biobase_2.44.0       BiocGenerics_0.30.0  Seurat_3.1.0        

loaded via a namespace (and not attached):
  [1] Rtsne_0.15          colorspace_1.4-1    rjson_0.2.20       
  [4] ggridges_0.5.1      mclust_5.4.5        circlize_0.4.7     
  [7] GlobalOptions_0.1.0 clue_0.3-57         farver_1.1.0       
 [10] leiden_0.3.1        listenv_0.7.0       npsurv_0.4-0       
 [13] graphlayouts_0.5.0  ggrepel_0.8.1       bit64_0.9-7        
 [16] codetools_0.2-16    splines_3.6.1       R.methodsS3_1.7.1  
 [19] lsei_1.2-0          knitr_1.24          polyclip_1.10-0    
 [22] zeallot_0.1.0       jsonlite_1.6        ica_1.0-2          
 [25] cluster_2.1.0       png_0.1-7           R.oo_1.22.0        
 [28] pheatmap_1.0.12     uwot_0.1.4          ggforce_0.3.1      
 [31] sctransform_0.2.0   compiler_3.6.1      httr_1.4.1         
 [34] backports_1.1.4     assertthat_0.2.1    lazyeval_0.2.2     
 [37] tweenr_1.0.1        htmltools_0.3.6     tools_3.6.1        
 [40] rsvd_1.0.2          igraph_1.2.4.1      gtable_0.3.0       
 [43] glue_1.3.1          RANN_2.6.1          reshape2_1.4.3     
 [46] Rcpp_1.0.2          vctrs_0.2.0         gdata_2.18.0       
 [49] ape_5.3             nlme_3.1-141        gbRd_0.4-11        
 [52] lmtest_0.9-37       xfun_0.9            stringr_1.4.0      
 [55] globals_0.12.4      lifecycle_0.1.0     irlba_2.3.3        
 [58] gtools_3.8.1        future_1.14.0       MASS_7.3-51.4      
 [61] zoo_1.8-6           scales_1.0.0        tidygraph_1.1.2    
 [64] RColorBrewer_1.1-2  yaml_2.2.0          memoise_1.1.0      
 [67] reticulate_1.13     pbapply_1.4-2       gridExtra_2.3      
 [70] stringi_1.4.3       RSQLite_2.1.2       checkmate_1.9.4    
 [73] caTools_1.17.1.2    bibtex_0.4.2        shape_1.4.4        
 [76] Rdpack_0.11-0       SDMTools_1.1-221.1  rlang_0.4.0        
 [79] pkgconfig_2.0.3     bitops_1.0-6        evaluate_0.14      
 [82] lattice_0.20-38     ROCR_1.0-7          purrr_0.3.2        
 [85] labeling_0.3        htmlwidgets_1.3     bit_1.1-14         
 [88] tidyselect_0.2.5    RcppAnnoy_0.0.13    plyr_1.8.4         
 [91] magrittr_1.5        R6_2.4.0            gplots_3.0.1.1     
 [94] DBI_1.0.0           withr_2.1.2         pillar_1.4.2       
 [97] fitdistrplus_1.0-14 survival_2.44-1.1   RCurl_1.95-4.12    
[100] tibble_2.1.3        future.apply_1.3.0  tsne_0.1-3         
[103] crayon_1.3.4        KernSmooth_2.23-16  plotly_4.9.0       
[106] rmarkdown_1.15      viridis_0.5.1       GetoptLong_0.1.7   
[109] data.table_1.12.2   blob_1.2.0          metap_1.1          
[112] digest_0.6.21       xtable_1.8-4        tidyr_1.0.0        
[115] R.utils_2.9.0       RcppParallel_4.4.3  munsell_0.5.0      
[118] viridisLite_0.3.0  
```



