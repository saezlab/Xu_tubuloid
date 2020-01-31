CK121 sorted CD24+ cells: initial clustering
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
SeuratObject <- getSeuratObject(path = "../data/sc/CK121_CD24/filtered_feature_bc_matrix/",
                     project_name = "CK121_CD24", 
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

![](1_initial_clustering_CK121_CD24_files/figure-gfm/QC-1.png)<!-- -->

``` r
print(CombinePlots(plots = list(FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "percent.mt"),
                                FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "nFeature_RNA")))
)
```

![](1_initial_clustering_CK121_CD24_files/figure-gfm/QC-2.png)<!-- -->

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

Doing so, 3571 out of an initial 3655 total were retrieve for the
analysis, thus discarding a 84 cells.

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

![](1_initial_clustering_CK121_CD24_files/figure-gfm/sel-1.png)<!-- -->

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
    ## Positive:  S100A8, FTL, S100A4, TYROBP, PLAUR, S100A12, CXCL8, IGKC, CD69, PLEK 
    ##     AIF1, MME, IGHM, DUSP2, FCER1G, APOBEC3A, MMP9, LYZ, AC079767.4, AREG 
    ##     IGHD, IL1B, IGHA1, AL928768.3, HSPA6, FCRL5, IGHG1, JCHAIN, GPR18, HCST 
    ## Negative:  IGFBP7, KRT18, TNFRSF12A, S100A10, KRT8, HSPB1, CD151, DSTN, ATP1B1, GSTP1 
    ##     APP, CLDN3, EPCAM, S100A13, SDC4, TPM1, CLU, MGST3, PTMS, CD59 
    ##     SOX4, LMNA, TSPAN1, CLDN4, ANXA2, MYO6, FXYD2, TMEM176B, ITGB8, CD9 
    ## PC_ 2 
    ## Positive:  LGALS2, GATM, MT1G, PDZK1, PDZK1IP1, NAT8, MGST1, UGT2B7, GLYAT, FXYD2 
    ##     CMBL, BHMT, FBP1, MT1F, CXCL14, SHMT1, BBOX1, BHMT2, FTCD, MIOX 
    ##     DPYS, ALDOB, PBLD, TMEM27, ECHS1, HRSP12, SMIM24, USH1C, GLYATL1, ACMSD 
    ## Negative:  PTPRB, ENG, RAMP2, FLT1, A2M, RAMP3, SOX18, SLC9A3R2, MEIS2, EPAS1 
    ##     IGFBP5, PLAT, CLEC14A, EMCN, ADGRL4, LDB2, KDR, GIMAP7, SLCO2A1, EGFL7 
    ##     MMRN2, IFI27, CDH5, PODXL, HEG1, PLPP1, ECSCR.1, CD34, RNASE1, KANK3 
    ## PC_ 3 
    ## Positive:  GATM, LGALS2, RBP5, NAT8, ALDOB, ASS1, PDZK1, GLYAT, BHMT, DPYS 
    ##     PDZK1IP1, MIOX, UGT2B7, PEPD, AK4, MGST1, MT1G, GPX3, FTCD, DDC 
    ##     DAB2, GLYATL1, KHK, GSTA1, FBP1, AKR1C3, IL32, APOE, HRSP12, EHHADH 
    ## Negative:  MUC1, ERBB4, MAL, WFDC2, KCNJ1, TMPRSS4, DEFB1, S100A2, TACSTD2, CLCNKB 
    ##     ITGB6, CDH16, ITM2C, TSPAN8, SLC12A1, MECOM, CLDN7, SCNN1A, ATP6V1B1, CA12 
    ##     PAPPA2, SFRP1, CLDN16, TFCP2L1, TMEM213, IGFBP2, TFAP2B, ITGA2, CASR, BCAM 
    ## PC_ 4 
    ## Positive:  BGN, COL8A1, PTN, CLDN1, ITGB8, CDH6, PDLIM3, SULT1C4, VCAM1, SBSPON 
    ##     SPON2, UGT2A3, VIM, FSTL3, CFH, COL6A2, C1S, CTGF, C1orf186, NSG1 
    ##     CDH2, BAMBI, MXRA8, NBL1, PALLD, COL12A1, PRUNE2, TNC, PAPLN, HRH1 
    ## Negative:  ALDOB, GSTA1, BBOX1, PCK1, MIOX, HPD, GSTA2, FABP1, CLCNKB, CA12 
    ##     OGDHL, PTH1R, TMEM213, GLYATL1, ATP6V1B1, PRODH2, ALDH6A1, SUCLG1, ATP6V0A4, APOE 
    ##     KCNJ1, ALB, GPD1, AK4, PEPD, DPYS, KHK, SLC13A3, MECOM, ECHS1 
    ## PC_ 5 
    ## Positive:  HIST1H4C, CD69, IGKC, IGHM, MT-CO1, ZNF331, MT-ATP6, MT-ND5, TMSB10, MT-ND3 
    ##     TNFSF9, AC079767.4, MT-ND4, MT-CO2, MT-CO3, MT-ND1, MT-CYB, JUN, IGHD, GADD45B 
    ##     NR4A1, IGHA1, AREG, JCHAIN, FCRL5, ARID5B, ID3, AL928768.3, GPR18, DUSP2 
    ## Negative:  S100A8, S100A11, ITM2B, TYROBP, PLAUR, S100A6, FTL, S100A12, CXCL8, S100A4 
    ##     MME, IER3, FGL2, AIF1, SERPINA1, EGR1, CEBPD, FCER1G, APOBEC3A, PLEK 
    ##     MMP9, VMP1, LYZ, CXCL1, RBP7, HSPA6, IL1B, KCNJ15, TIMP2, ADM

``` r
print(ElbowPlot(SeuratObject,ndims = 50) + geom_vline(xintercept = 25, col="red"))
```

![](1_initial_clustering_CK121_CD24_files/figure-gfm/pca-1.png)<!-- -->

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
    ## Number of nodes: 3571
    ## Number of edges: 155609
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9670
    ## Number of communities: 6
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3571
    ## Number of edges: 155609
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9394
    ## Number of communities: 10
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3571
    ## Number of edges: 155609
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9142
    ## Number of communities: 10
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3571
    ## Number of edges: 155609
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8914
    ## Number of communities: 10
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3571
    ## Number of edges: 155609
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8709
    ## Number of communities: 11
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3571
    ## Number of edges: 155609
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8511
    ## Number of communities: 11
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3571
    ## Number of edges: 155609
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8329
    ## Number of communities: 12
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3571
    ## Number of edges: 155609
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8180
    ## Number of communities: 12
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3571
    ## Number of edges: 155609
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8045
    ## Number of communities: 12
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3571
    ## Number of edges: 155609
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7927
    ## Number of communities: 14
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

![](1_initial_clustering_CK121_CD24_files/figure-gfm/clustree-1.png)<!-- -->

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

    ## 22:27:41 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 22:27:41 Read 3571 rows and found 25 numeric columns

    ## 22:27:41 Using Annoy for neighbor search, n_neighbors = 30

    ## 22:27:41 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 22:27:42 Writing NN index file to temp file /tmp/RtmpmUm5lH/file3801171563c4
    ## 22:27:42 Searching Annoy index using 1 thread, search_k = 3000
    ## 22:27:43 Annoy recall = 100%
    ## 22:27:43 Commencing smooth kNN distance calibration using 1 thread
    ## 22:27:44 Initializing from normalized Laplacian + noise
    ## 22:27:44 Commencing optimization for 500 epochs, with 157668 positive edges
    ## 22:27:51 Optimization finished

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
