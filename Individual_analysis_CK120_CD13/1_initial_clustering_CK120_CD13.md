CK120 sorted CD13+ cells: initial clustering
================
Javier Perales-Paton - <javier.perales@bioquant.uni-heidelberg.de>

## Load libraries and auxiliar functions

``` r
set.seed(1234)
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(clustree))
suppressPackageStartupMessages(require(cowplot))
source("../src/seurat_fx.R")
```

## Read data into a Seurat Object

``` r
SeuratObject <- getSeuratObject(path = "../data/sc/CK120_CD13/filtered_feature_bc_matrix/",
                     project_name = "CK120_CD13", 
                     mt.pattern = "^MT-", min.cells = 5, 
                     min.features = 200)
```

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

## Define output directory

``` r
# Define output directory
OUTDIR <- paste0("./output/1_initial_clustering/")
if(! dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
```

## Quality control

``` r
print(VlnPlot(SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
```

![](1_initial_clustering_CK120_CD13_files/figure-gfm/QC-1.png)<!-- -->

``` r
print(CombinePlots(plots = list(FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "percent.mt"),
                                FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "nFeature_RNA")))
)
```

![](1_initial_clustering_CK120_CD13_files/figure-gfm/QC-2.png)<!-- -->

## Data cleasing: gene and cell filtering

To avoid unnecessary sparsity and noise in downstream analysis,  
\* We discard cells with low (\<200) and high (\>6000) number of genes
detected to avoid bad quality and doublets, repectively.  
\* We discard genes not detected across cells.  
For this we use cutoffs in concordance with kidney tissue, where
Proximal tubule cells might have a high content of mitochondrial genes
(\<80%). Later we will apply diagnostics to see if this is related.

``` r
nSamples_before <- ncol(SeuratObject)
SeuratObject <- subset(SeuratObject, nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 80)
nSamples_after <- ncol(SeuratObject)
```

Doing so, 7233 out of an initial 7279 total were retrieve for the
analysis, thus discarding a 46 cells.

## Pre-processing the data for cell clustering and cell-type assignment

First, we normalize the data using default parameters

``` r
SeuratObject <- NormalizeData(SeuratObject)
```

Second, we perform a feature selection based on `vst` method to select
high variable genes (top 2000) for cell
clustering.

``` r
SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 2000)

print(LabelPoints(VariableFeaturePlot(SeuratObject),
                  points=VariableFeatures(SeuratObject)[1:10],repel = TRUE))
```

    ## Warning: Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
    ## Please use `as_label()` or `as_name()` instead.
    ## This warning is displayed once per session.

    ## When using repel, set xnudge and ynudge to 0 for optimal results

![](1_initial_clustering_CK120_CD13_files/figure-gfm/sel-1.png)<!-- -->

Third, we center and scale the data prior PCA:

``` r
SeuratObject <- Seurat_scaledata(SeuratObject)
```

    ## Centering and scaling data matrix

Finally we perform PCA, selecting a number of PCs that are relevant
based on the elbow plot. We tested that our results are not sensitive to
this parameter, and the arbitrary selection of number of PCs do not
change final
conclusions.

``` r
SeuratObject <- RunPCA(SeuratObject, features = VariableFeatures(SeuratObject), npcs = 50)
```

    ## PC_ 1 
    ## Positive:  TMSB10, S100A11, TMSB4X, VIM, RPL10, IL32, PTMA, RPL13, RPS12, RPL3 
    ##     KRT18, RPL12, S100A10, S100A6, RPS8, MYL6, RPLP1, RPS18, RPS23, RPS2 
    ##     RPL15, RPL8, RPL18A, ACTB, GSTP1, RPS24, RPL7A, RPL32, RPL28, RPS3A 
    ## Negative:  MT-ND2, GPX3, GSTA1, UGT2B7, GSTA2, MTRNR2L8, SLC22A8, FABP1, APOM, MIOX 
    ##     PLG, PRAP1, AOC1, NAT8, CTSB, PCK1, LGMN, SLC7A8, AZGP1, ALB 
    ##     SLC27A2, SLC5A12, SLC4A4, AQP1, SLC16A9, CA12, ATP1B1, MCCD1, CYP4A11, G6PC 
    ## PC_ 2 
    ## Positive:  GAPDH, MT1G, IL32, CRYAB, KRT18, RPS6, RPL3, MT1F, RPL7, BBOX1 
    ##     RPL13, KRT8, RPL12, C12orf75, ENO1, ALDH1A1, GNB2L1, RPL14, RPL8, PPA1 
    ##     ANXA4, RPL29, RPS4X, RPL26, RPL10, AKR7A3, RPL10A, RPL7A, RPL15, MT1H 
    ## Negative:  NEAT1, MALAT1, SRGN, HLA-DPA1, GMFG, KLF2, PLEKHO1, SLC2A3, ARHGDIB, HLA-DQA1 
    ##     ANXA1, FXYD5, LAPTM5, MEF2C, PECAM1, HLA-DPB1, HLA-DQB1, RNASE1, RSRP1, CTD-3252C9.4 
    ##     NR4A2, BST2, A2M, TIMP1, FOSB, IFI16, FCER1G, CTSS, THBD, PTPRC 
    ## PC_ 3 
    ## Positive:  IGFBP7, GSTP1, SPARC, NUPR1, THY1, S100A6, TUBA1A, GNG11, TMSB4X, ANK2 
    ##     TMSB10, PRAP1, ADGRG1, SPON2, S100A10, RP11-149I23.3, TMEM54, SLC7A8, TAGLN2, TNFRSF12A 
    ##     TXNIP, PPFIBP1, IFITM3, APP, HSPB8, LGMN, PFKP, CLDN4, RARRES2, B2M 
    ## Negative:  PCK1, RBP4, MIOX, DCXR, GLYATL1, BHMT, GPT, HAO2, FMO5, BBOX1 
    ##     SLC6A18, ALDH6A1, AGXT, MT1G, CYP17A1, PTER, RP11-536I6.2, SLC7A13, MYH8, AOX1 
    ##     MT1H, UPP2, MSRB1, AKR7A3, MPC1, SLC22A7, ACY3, ABHD14A, SLC5A10, GAPDH 
    ## PC_ 4 
    ## Positive:  IGFBP5, SLC9A3R2, IFI27, EGFL7, RAMP2, ID1, MGP, TM4SF1, C8orf4, MEIS2 
    ##     TIMP3, NPDC1, RAMP3, EMCN, SOX18, KANK3, TMEM204, CLEC14A, PTPRB, TGFBR2 
    ##     KDR, GIMAP7, PLPP3, CRHBP, PTRF, PLPP1, ADGRF5, UACA, PLAT, SDPR 
    ## Negative:  AIF1, TYROBP, FCER1G, PLEK, CTSS, LYZ, SPI1, HCST, MS4A6A, PTPRC 
    ##     CD37, MNDA, LCP1, CXCR4, LAPTM5, PLAUR, ITGB2, PYCARD, CORO1A, RGS10 
    ##     DUSP2, GPR183, CELF2, PHACTR1, LST1, HLA-DQB1, CD4, CSF1R, C1QA, RASSF5 
    ## PC_ 5 
    ## Positive:  GPX3, UGT2B7, SLC22A8, CTSB, SLC7A8, SLC16A9, SLC4A4, ATP1B1, AQP1, NAT8 
    ##     LGMN, THY1, PRAP1, IGFBP7, CA12, PSAP, ADGRG1, SPP1, SLC5A12, ANK2 
    ##     CTSD, FABP1, MAF, GSTA1, CST3, AOC1, TSPAN1, G0S2, GSTA2, GPAT3 
    ## Negative:  TMSB4X, VIM, S100A6, RNASE1, MYH8, CRNDE, ID3, MARCKS, RBP4, MARCKSL1 
    ##     GYPC, RP5-1021I20.2, AKAP12, S100A11, ANXA2, IGFBP5, IFI27, RPL10, MGP, KLF2 
    ##     ID1, RPS18, RPL29, SLC9A3R2, RPS12, GIMAP7, NEU4, RP11-536I6.2, RPS24, ANXA3

``` r
print(ElbowPlot(SeuratObject,ndims = 50) + geom_vline(xintercept = 25, col="red"))
```

![](1_initial_clustering_CK120_CD13_files/figure-gfm/pca-1.png)<!-- -->

We decided to use 25 PCs for the clustering.

## Cell clustering

We are going to use Shared-Nearest Neighbour with Graph partitioning for
cell clustering. We run the algorithm with multiple resolutions. Later
this is read-out to give an idea how these are consistent from a lower
to high number of expected cell populations in the
sample.

``` r
SeuratObject <- FindNeighbors(SeuratObject, dims = 1:25)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
SeuratObject <- FindClusters(SeuratObject, resolution = seq(from=0.1, to=1, by=0.1))
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7233
    ## Number of edges: 269546
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9105
    ## Number of communities: 3
    ## Elapsed time: 1 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7233
    ## Number of edges: 269546
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8619
    ## Number of communities: 4
    ## Elapsed time: 1 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7233
    ## Number of edges: 269546
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8287
    ## Number of communities: 5
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7233
    ## Number of edges: 269546
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8062
    ## Number of communities: 6
    ## Elapsed time: 1 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7233
    ## Number of edges: 269546
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7861
    ## Number of communities: 7
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7233
    ## Number of edges: 269546
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7686
    ## Number of communities: 7
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7233
    ## Number of edges: 269546
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7504
    ## Number of communities: 7
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7233
    ## Number of edges: 269546
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7350
    ## Number of communities: 9
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7233
    ## Number of edges: 269546
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7235
    ## Number of communities: 11
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7233
    ## Number of edges: 269546
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7112
    ## Number of communities: 11
    ## Elapsed time: 0 seconds

``` r
# We redefine the final partitioning with resolution 0.5
SeuratObject$seurat_clusters <- SeuratObject$RNA_snn_res.0.5 
Idents(SeuratObject) <- SeuratObject$RNA_snn_res.0.5 
```

We will investigate how the selection of multiple resolutions affects
the partition into individual cell
clusters.

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.")
```

![](1_initial_clustering_CK120_CD13_files/figure-gfm/clustree-1.png)<!-- -->

We decided to chose resolution 0.5 for the initial clustering. We will
save this output to archive the outcome. Herein, first column is the
cell barcode, and rest of columns are clusters where each cell belong to
across multiple resolutions. The ones that are active for the study are
those related to resolution 0.5 (`RNA_snn_res.0.5` and
`seurat_clusters`).

``` r
# 1 Clustering outcome
write.table(SeuratObject@meta.data[,c(grep("^RNA_snn_res",
                                           colnames(SeuratObject@meta.data),
                                           value=TRUE),
                                      "seurat_clusters"),],
            file=paste0(OUTDIR,"/init_clustering.tsv"),
            sep="\t", col.names = NA, row.names=TRUE, quote=TRUE)

# 2 Initial idents (same as seurat_clusters)
write.table(data.frame("Ident"=SeuratObject@active.ident),
            file=paste0(OUTDIR,"/active_idents.tsv"),
            sep="\t", col.names = NA, row.names = TRUE, quote=TRUE)
```

## Non-linear dim reduction (umap)

``` r
SeuratObject <- RunUMAP(SeuratObject, dims = 1:25)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 17:17:37 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 17:17:37 Read 7233 rows and found 25 numeric columns

    ## 17:17:37 Using Annoy for neighbor search, n_neighbors = 30

    ## 17:17:37 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 17:17:38 Writing NN index file to temp file /tmp/RtmprUC1lT/file6a8e79060fdb
    ## 17:17:38 Searching Annoy index using 1 thread, search_k = 3000
    ## 17:17:40 Annoy recall = 100%
    ## 17:17:41 Commencing smooth kNN distance calibration using 1 thread
    ## 17:17:41 Initializing from normalized Laplacian + noise
    ## 17:17:41 Commencing optimization for 500 epochs, with 324432 positive edges
    ## 17:18:00 Optimization finished

## Archive processed data for downstream analysis

``` r
DATA_DIR <- paste0(OUTDIR,"/data")
if(!dir.exists(DATA_DIR)) dir.create(DATA_DIR)
```

``` r
saveRDS(SeuratObject, paste0(DATA_DIR,"/SeuratObject.rds"))
```

## Session info

``` r
sessionInfo()
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] cowplot_1.0.0     clustree_0.4.1    ggraph_2.0.0.9000 ggplot2_3.2.1    
    ## [5] Seurat_3.1.0     
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] tsne_0.1-3          nlme_3.1-141        bitops_1.0-6       
    ##   [4] RcppAnnoy_0.0.13    RColorBrewer_1.1-2  httr_1.4.1         
    ##   [7] sctransform_0.2.0   tools_3.6.1         backports_1.1.4    
    ##  [10] R6_2.4.0            irlba_2.3.3         KernSmooth_2.23-16 
    ##  [13] uwot_0.1.4          lazyeval_0.2.2      colorspace_1.4-1   
    ##  [16] withr_2.1.2         npsurv_0.4-0        gridExtra_2.3      
    ##  [19] tidyselect_0.2.5    compiler_3.6.1      plotly_4.9.0       
    ##  [22] labeling_0.3        checkmate_1.9.4     caTools_1.17.1.2   
    ##  [25] scales_1.0.0        lmtest_0.9-37       ggridges_0.5.1     
    ##  [28] pbapply_1.4-2       stringr_1.4.0       digest_0.6.21      
    ##  [31] rmarkdown_1.15      R.utils_2.9.0       pkgconfig_2.0.3    
    ##  [34] htmltools_0.3.6     bibtex_0.4.2        htmlwidgets_1.3    
    ##  [37] rlang_0.4.0         farver_1.1.0        zoo_1.8-6          
    ##  [40] jsonlite_1.6        ica_1.0-2           gtools_3.8.1       
    ##  [43] dplyr_0.8.3         R.oo_1.22.0         magrittr_1.5       
    ##  [46] Matrix_1.2-17       Rcpp_1.0.2          munsell_0.5.0      
    ##  [49] viridis_0.5.1       ape_5.3             reticulate_1.13    
    ##  [52] lifecycle_0.1.0     R.methodsS3_1.7.1   stringi_1.4.3      
    ##  [55] yaml_2.2.0          gbRd_0.4-11         MASS_7.3-51.4      
    ##  [58] gplots_3.0.1.1      Rtsne_0.15          plyr_1.8.4         
    ##  [61] grid_3.6.1          parallel_3.6.1      gdata_2.18.0       
    ##  [64] listenv_0.7.0       ggrepel_0.8.1       crayon_1.3.4       
    ##  [67] lattice_0.20-38     graphlayouts_0.5.0  splines_3.6.1      
    ##  [70] SDMTools_1.1-221.1  zeallot_0.1.0       knitr_1.24         
    ##  [73] pillar_1.4.2        igraph_1.2.4.1      future.apply_1.3.0 
    ##  [76] reshape2_1.4.3      codetools_0.2-16    leiden_0.3.1       
    ##  [79] glue_1.3.1          evaluate_0.14       lsei_1.2-0         
    ##  [82] metap_1.1           RcppParallel_4.4.3  data.table_1.12.2  
    ##  [85] tweenr_1.0.1        vctrs_0.2.0         png_0.1-7          
    ##  [88] Rdpack_0.11-0       polyclip_1.10-0     gtable_0.3.0       
    ##  [91] RANN_2.6.1          purrr_0.3.2         tidyr_1.0.0        
    ##  [94] future_1.14.0       assertthat_0.2.1    ggforce_0.3.1      
    ##  [97] xfun_0.9            rsvd_1.0.2          tidygraph_1.1.2    
    ## [100] RSpectra_0.15-0     survival_2.44-1.1   viridisLite_0.3.0  
    ## [103] tibble_2.1.3        cluster_2.1.0       globals_0.12.4     
    ## [106] fitdistrplus_1.0-14 ROCR_1.0-7

``` r
{                                                                                                                                                                                                           
sink(file=paste0(OUTDIR,"/sessionInfo.txt"))
print(sessionInfo())
sink()
}
```
