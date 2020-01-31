CK119 late organoid: initial clustering
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
SeuratObject <- getSeuratObject(path = "../data/sc/CK119_organoid/filtered_feature_bc_matrix/",
                     project_name = "CK119_organoid", 
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

![](1_initial_clustering_CK119_organoid_files/figure-gfm/QC-1.png)<!-- -->

``` r
print(CombinePlots(plots = list(FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "percent.mt"),
                                FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "nFeature_RNA")))
)
```

![](1_initial_clustering_CK119_organoid_files/figure-gfm/QC-2.png)<!-- -->

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

Doing so, 5328 out of an initial 5681 total were retrieve for the
analysis, thus discarding a 353 cells.

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

![](1_initial_clustering_CK119_organoid_files/figure-gfm/sel-1.png)<!-- -->

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
    ## Positive:  H2AFZ, TUBA1B, TUBA1C, TUBB, RAN, RANBP1, UBE2S, STMN1, HMGB1, PTTG1 
    ##     KIAA0101, SNRPD1, CKS1B, DTYMK, HSPD1, KPNA2, NME1, TYMS, CKS2, TUBB4B 
    ##     C1QBP, HNRNPAB, DEK, S100A16, TNFRSF12A, NCL, CACYBP, PHF19, HMGN2, CDKN3 
    ## Negative:  NEAT1, C3, CP, MUC1, LTF, COL4A3, AC019181.2, SCN2A, NDRG1, KCNQ1OT1 
    ##     IGFBP3, TNFAIP2, SAT1, FAM20A, PIGR, CA12, CMYA5, FOS, SORL1, LINC01320 
    ##     ELF3, CLU, CPAMD8, CFTR, MUC16, TMEM91, RNF213, OLFM4, C8orf4, ANGPTL4 
    ## PC_ 2 
    ## Positive:  CLDN7, KLK6, TACSTD2, LCN2, KRT7, PERP, MAL2, RAB25, CST6, C15orf48 
    ##     RAB11FIP1, ISG20, KRT19, S100A14, SAT1, HEBP2, PTGES, SERPINB1, SDCBP2, CLDN4 
    ##     SH3BGRL3, SMIM5, C9orf16, S100A9, GCNT3, LGALS3, CSTB, TPM4, KRT8, KRT18 
    ## Negative:  VCAN, SPP1, VIM, PLEKHA1, ADAMTS1, CLU, MYLK, HTRA1, TFAP2B, SFRP1 
    ##     AC019181.2, CENPF, CRNDE, LINC01320, KCNQ1OT1, LTF, IGFBP2, MT1E, NNMT, MKI67 
    ##     TPM2, GYPC, SCN2A, ASPM, CP, FGFR1, NTRK2, PEG10, TINAG, COL4A3 
    ## PC_ 3 
    ## Positive:  TUBA1A, PLAU, EMP3, VIM, NEFL, AKAP12, OCIAD2, SFRP1, IGFBP7, UCHL1 
    ##     ADIRF, LGALS1, AKR1B1, CRYM, GLRX, IGFBP6, FILIP1L, HTRA1, GYPC, DCBLD2 
    ##     PPP1R1A, MAP1B, TPM1, S100A6, MAL, CAV2, CCND1, SERPINE2, FTL, CAV1 
    ## Negative:  MKI67, TOP2A, CENPF, ASPM, NUSAP1, HMGB2, DLGAP5, CCNA2, CENPE, GTSE1 
    ##     HMMR, TPX2, UBE2C, PRC1, PLK1, CCNB1, ANLN, CDC20, SMC4, CCNB2 
    ##     CASC5, NUF2, NDC80, NCAPG, KIF23, CDCA3, CDK1, RRM2, DEPDC1, CENPA 
    ## PC_ 4 
    ## Positive:  FGFBP1, LAMC2, KRT8, AREG, F3, KRT18, EMP1, PLEC, MALT1, LAMA3 
    ##     EREG, PLAT, GPRC5A, SCEL, ARHGAP29, ITGA2, WNT7A, MALL, ITGB4, S100A14 
    ##     EZR, FLNB, PADI1, PHLDA2, IL1RL1, C6orf132, THBS1, CYR61, DCBLD2, GAL 
    ## Negative:  BTG1, TXNIP, GDF15, FTL, TMEM140, SQSTM1, PNRC1, HLA-B, SELM, IFIT3 
    ##     PSMB9, RSAD2, IFIT1, PIK3IP1, HERC5, MX1, H1FX, IFIT2, GBP2, FTH1 
    ##     CDKN1A, AKR1C3, BBC3, IFI6, G0S2, YPEL3, UBE2L6, CEBPB, B2M, DDIT3 
    ## PC_ 5 
    ## Positive:  SPP1, WFDC2, AGR2, CLU, RARRES2, CPVL, FAM162A, ATP1B1, LDHB, PDZK1IP1 
    ##     CYB5A, SERPINA1, GSTM3, OLFM4, TNFRSF11B, BBOX1, DEFB1, SLPI, KITLG, KRT19 
    ##     RARRES1, TSPAN1, NNMT, TNFSF10, CP, TCN1, RNASET2, SPAG4, IGFBP7, ID2 
    ## Negative:  RSAD2, ATF3, PMAIP1, ISG15, IL32, MX1, KLF6, IFIT3, CDKN2B, IFIT2 
    ##     ISG20, CDKN1A, IFIT1, TMEM140, PIK3IP1, HERC5, PLAU, THSD7A, C15orf48, SQSTM1 
    ##     CMPK2, IFI27, OAS1, OPTN, OAS2, XAF1, FN1, MAP1B, IL1RN, OAS3

``` r
print(ElbowPlot(SeuratObject,ndims = 50) + geom_vline(xintercept = 25, col="red"))
```

![](1_initial_clustering_CK119_organoid_files/figure-gfm/pca-1.png)<!-- -->

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
    ## Number of nodes: 5328
    ## Number of edges: 187682
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9501
    ## Number of communities: 6
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5328
    ## Number of edges: 187682
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9288
    ## Number of communities: 7
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5328
    ## Number of edges: 187682
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9122
    ## Number of communities: 9
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5328
    ## Number of edges: 187682
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8977
    ## Number of communities: 10
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5328
    ## Number of edges: 187682
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8849
    ## Number of communities: 10
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5328
    ## Number of edges: 187682
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8726
    ## Number of communities: 10
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5328
    ## Number of edges: 187682
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8610
    ## Number of communities: 11
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5328
    ## Number of edges: 187682
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8500
    ## Number of communities: 13
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5328
    ## Number of edges: 187682
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8417
    ## Number of communities: 15
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5328
    ## Number of edges: 187682
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8329
    ## Number of communities: 16
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

![](1_initial_clustering_CK119_organoid_files/figure-gfm/clustree-1.png)<!-- -->
\#\# Non-linear dim reduction
    (umap)

``` r
SeuratObject <- RunUMAP(SeuratObject, dims = 1:25)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 15:28:19 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 15:28:19 Read 5328 rows and found 25 numeric columns

    ## 15:28:19 Using Annoy for neighbor search, n_neighbors = 30

    ## 15:28:19 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 15:28:19 Writing NN index file to temp file /tmp/RtmpXY95Sb/file7fd6557cf64f
    ## 15:28:19 Searching Annoy index using 1 thread, search_k = 3000
    ## 15:28:21 Annoy recall = 100%
    ## 15:28:21 Commencing smooth kNN distance calibration using 1 thread
    ## 15:28:22 Initializing from normalized Laplacian + noise
    ## 15:28:22 Commencing optimization for 500 epochs, with 224316 positive edges
    ## 15:28:33 Optimization finished

``` r
SeuratObject <- RunTSNE(SeuratObject, dims = 1:25)
```

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
    ##  [82] metap_1.1           RcppParallel_4.4.3  data.table_1.12.8  
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
