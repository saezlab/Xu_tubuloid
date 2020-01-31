CK5 early organoid: initial clustering
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
SeuratObject <- getSeuratObject(path = "../data/sc/Organoid/filtered_feature_bc_matrix/",
                     project_name = "CK5_organoid", 
                     mt.pattern = "^MT-", min.cells = 5, 
                     min.features = 200)
```

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

![](1_initial_clustering_CK5_organoid_files/figure-gfm/QC-1.png)<!-- -->

``` r
print(CombinePlots(plots = list(FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "percent.mt"),
                                FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "nFeature_RNA")))
)
```

![](1_initial_clustering_CK5_organoid_files/figure-gfm/QC-2.png)<!-- -->

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

Doing so, 2736 out of an initial 2792 total were retrieve for the
analysis, thus discarding a 56 cells.

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

![](1_initial_clustering_CK5_organoid_files/figure-gfm/sel-1.png)<!-- -->

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
    ## Positive:  NEAT1, KCNQ1OT1, CLU, C3, LTF, CP, MUC1, ADAMTS1, SPP1, SLC34A2 
    ##     PIGR, SLC3A1, C1R, SERPING1, CDH16, COL4A3, NCOA7, C1S, SAA1, SAA2 
    ##     TNFAIP2, SLPI, IFI44L, IL4I1, SOD2, SCTR, WFDC2, BIRC3, GPNMB, LCN2 
    ## Negative:  TUBA1B, CFL1, RAN, H2AFZ, TUBB, S100A16, JPT1, HMGA1, TXN, ACTB 
    ##     TMSB4X, PFN1, OAZ1, SNRPB, CALM1, KRT18, ARPC2, RANBP1, HMGB1, C12orf75 
    ##     LSM4, RPS26, EIF5A, TUBB4B, MYL12A, PHLDA2, GSTP1, STMN1, CCND1, SNRPD1 
    ## PC_ 2 
    ## Positive:  FN1, IL32, SPARC, CDH6, MAP1B, TPM1, RAB3B, FLG, TMSB10, KRT23 
    ##     TUBA1A, SLIT3, THY1, CDH2, MFGE8, DCBLD2, S100A6, PLA2G16, CALD1, AKAP12 
    ##     EMP3, RPS27L, SH3BGRL3, OCIAD2, RHOC, CRYAB, SPOCK2, ITIH5, CYTL1, GLIPR2 
    ## Negative:  MKI67, HMGB2, CENPF, TOP2A, NUSAP1, CDC20, UBE2C, PBK, CDCA3, DLGAP5 
    ##     GTSE1, PLK1, HMMR, DEPDC1, TPX2, ASPM, BIRC5, AURKB, CCNB1, CCNB2 
    ##     CDKN3, CCNA2, SPC25, RRM2, ANLN, NUF2, CEP55, NEK2, CENPE, TK1 
    ## PC_ 3 
    ## Positive:  SPARC, FLG, THY1, CDH2, SLIT3, CENPF, PTGIS, ITIH5, CD70, MKI67 
    ##     CYTL1, MAP1B, ASPM, SEPT4, CRB2, IGFBP4, VCAM1, TOP2A, THBS2, NUSAP1 
    ##     GTSE1, FN1, DEPDC1, CDC20, DPYSL3, CD14, TM4SF4, PLK1, BCAT1, HMMR 
    ## Negative:  TSPAN1, RAB25, CLDN7, S100A14, EPCAM, AGR2, TACSTD2, KRT7, WFDC2, TMPRSS4 
    ##     SAT1, GPRC5A, CLDN4, MUC1, KITLG, PRSS8, UCP2, SFTA2, KRT19, GATA3 
    ##     CST6, SMIM22, ITGB6, CDH1, MPZL2, ITGA2, SFN, KRTCAP3, MAL2, KLK6 
    ## PC_ 4 
    ## Positive:  FLG, THY1, SLIT3, KRT23, PTGIS, ITIH5, SPARC, CYTL1, CD14, C6orf99 
    ##     CRB2, IGFBP4, SEPT4, AP1M2, UPK3B, FGF18, PRRX2, KRT19, TMOD1, BCAT1 
    ##     UPK1B, SNCA, CFTR, GATA3, RGS5, SOX17, LHX1, NGFR, AC005482.1, TM4SF4 
    ## Negative:  SERPINE2, MT2A, GYPC, PLAU, EMP3, CXCL14, VIM, GMDS, MT1E, PHLDA1 
    ##     GLRX, CLDN16, XKR4, TRIM55, ADIRF, TM4SF1, LGALS1, TPM1, PKIB, IGFBP6 
    ##     SAA1, HTRA1, ANPEP, GCHFR, SLC34A2, MSMP, TSPAN8, CRYAB, NNMT, S100A6 
    ## PC_ 5 
    ## Positive:  DBI, LDHB, ANXA4, SPP1, CLU, NNMT, SERPINA1, WFDC2, FTL, LDHA 
    ##     SLPI, PDZK1IP1, RPS2, RPL21, SOD2, GSTM3, TMEM37, DEFB1, TNFSF10, LCN2 
    ##     HLA-DMB, MYL12B, ENO1, PRDX1, HIST1H4C, AGR2, S100A11, PGK1, ID2, CLDN10 
    ## Negative:  PLAU, PLEC, AREG, LAMC2, DCBLD2, KLF6, F3, MFGE8, WNT7A, AKAP12 
    ##     UCA1, AL161431.1, EREG, KRTAP2-3, HMGA2, HSPG2, SCEL, RAB3B, LAMA3, PXDN 
    ##     PLAT, CCND1, ITGB8, ITGB4, CST6, PMEPA1, ARL4C, TRBC2, GAS6, ITGA2

``` r
print(ElbowPlot(SeuratObject,ndims = 50) + geom_vline(xintercept = 25, col="red"))
```

![](1_initial_clustering_CK5_organoid_files/figure-gfm/pca-1.png)<!-- -->

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
    ## Number of nodes: 2736
    ## Number of edges: 95837
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9422
    ## Number of communities: 4
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2736
    ## Number of edges: 95837
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9133
    ## Number of communities: 6
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2736
    ## Number of edges: 95837
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8872
    ## Number of communities: 8
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2736
    ## Number of edges: 95837
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8673
    ## Number of communities: 9
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2736
    ## Number of edges: 95837
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8519
    ## Number of communities: 10
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2736
    ## Number of edges: 95837
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8389
    ## Number of communities: 10
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2736
    ## Number of edges: 95837
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8269
    ## Number of communities: 11
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2736
    ## Number of edges: 95837
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8152
    ## Number of communities: 12
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2736
    ## Number of edges: 95837
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8053
    ## Number of communities: 12
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2736
    ## Number of edges: 95837
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7950
    ## Number of communities: 13
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

![](1_initial_clustering_CK5_organoid_files/figure-gfm/clustree-1.png)<!-- -->

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

    ## 15:11:15 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 15:11:15 Read 2736 rows and found 25 numeric columns

    ## 15:11:15 Using Annoy for neighbor search, n_neighbors = 30

    ## 15:11:15 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 15:11:16 Writing NN index file to temp file /tmp/RtmpKkEMZP/file61bb2e8fde1f
    ## 15:11:16 Searching Annoy index using 1 thread, search_k = 3000
    ## 15:11:16 Annoy recall = 100%
    ## 15:11:17 Commencing smooth kNN distance calibration using 1 thread
    ## 15:11:17 Initializing from normalized Laplacian + noise
    ## 15:11:17 Commencing optimization for 500 epochs, with 111076 positive edges
    ## 15:11:23 Optimization finished

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
