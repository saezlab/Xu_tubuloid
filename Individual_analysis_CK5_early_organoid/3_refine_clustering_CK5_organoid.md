CK5 early organoid : Refine cell clustering
================
Javier Perales-Paton - <javier.perales@bioquant.uni-heidelberg.de>

## Load libraries and auxiliar functions

``` r
set.seed(1234)
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(GSEABase))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(genesorteR))
suppressPackageStartupMessages(require(ComplexHeatmap))
suppressPackageStartupMessages(require(clustree))
suppressPackageStartupMessages(require(cowplot))
source("../src/seurat_fx.R")
```

## Load SeuratObject with initial clustering outcome

``` r
SeuratObject <- readRDS("./output/2_cell_assignment/data/SeuratObject.rds")
```

## Define output directory

``` r
# Define output directory
OUTDIR <- paste0("./output/3_refine_clustering/")
if(! dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
```

## Re-clustering and high dim. reduction after filtering contaminant populations

To trace back previous cell boundaries, we store the outcome

``` r
if(!"init_seurat_clusters" %in% colnames(SeuratObject@meta.data))
  SeuratObject$init_seurat_clusters <- SeuratObject$seurat_clusters

if(!"init_assign" %in% colnames(SeuratObject@meta.data))
  SeuratObject$init_assign <- Idents(SeuratObject)
```

We remove two cell populations that present with over-expression of
genes related to bad quality cells. We observed that cluster 5
over-express high mitochondrial genes, whereas cluster 9 over-express
ribosomal genes.

``` r
# Remove original reference to clustering
SeuratObject <- SeuratObject[, !grepl(":(High|Low)", 
                                      Idents(SeuratObject))]
```

We proceed with the standard pipeline for cell clustering, over-written
previous outcome

``` r
## Feature selection
SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 2000)

## PCA 
SeuratObject <- RunPCA(SeuratObject, features = VariableFeatures(SeuratObject), npcs = 50)
```

    ## PC_ 1 
    ## Positive:  TUBA1B, H2AFZ, CENPW, BIRC5, TUBB, TK1, HMGA1, RAN, CDKN3, DTYMK 
    ##     CKS1B, KRT18, JPT1, TUBB4B, RANBP1, ACTB, LSM4, CENPM, CFL1, PCLAF 
    ##     S100A16, TUBA1C, CCND1, CDC20, KRT8, MKI67, SRSF3, ARPC2, CENPF, UBE2S 
    ## Negative:  CLU, SPP1, NEAT1, SLC34A2, LTF, C3, SLPI, MUC1, IFITM3, WFDC2 
    ##     SOD2, SERPING1, ADAMTS1, BCAM, CP, KCNQ1OT1, SAA1, SAA2, SLC3A1, FTL 
    ##     ANXA4, C1R, NCOA7, SERPINA1, TMEM37, CDH16, TFPI, LCN2, PIGR, GPNMB 
    ## PC_ 2 
    ## Positive:  FN1, MAP1B, SPARC, RAB3B, CDH6, IL32, TPM1, TMSB10, DCBLD2, FLG 
    ##     CDH2, KRT23, THY1, MFGE8, AKAP12, TUBA1A, SLIT3, ITIH5, SPOCK2, CALD1 
    ##     C6orf99, PTGIS, CYTL1, CD70, THBS2, IGFBP4, GLIPR2, TM4SF4, S100A6, RHOC 
    ## Negative:  WFDC2, HMGB2, AGR2, DEFB1, MKI67, NUSAP1, CENPF, BCAM, TOP2A, MNS1 
    ##     CDCA3, SPP1, PBK, CDC20, LCN2, DLGAP5, AURKB, CPVL, UBE2C, PLK1 
    ##     HMMR, DEPDC1, GTSE1, ANP32E, CCNB1, BIRC5, TPX2, MUC1, ASPM, AGR3 
    ## PC_ 3 
    ## Positive:  CLDN7, TSPAN1, RAB25, S100A14, KRT7, TACSTD2, GPRC5A, TMPRSS4, CLDN4, ITGA2 
    ##     EPCAM, PRSS8, CST6, CDH1, KRT19, UCA1, MAL2, MACC1, SCEL, UCP2 
    ##     SFTA2, F3, PROM2, ITGB6, ST14, FOLR3, VTCN1, CDA, SFN, GATA3 
    ## Negative:  SPARC, THY1, FLG, CENPF, CDH2, ITIH5, MKI67, CD70, NUSAP1, VCAM1 
    ##     SLIT3, ASPM, DEPDC1, HMGB2, GTSE1, TOP2A, PTGIS, CDC20, PTTG1, PLK1 
    ##     CRB2, C1R, CYTL1, SPC25, CDKN3, KIF4A, HMMR, UBE2C, PBK, IGFBP4 
    ## PC_ 4 
    ## Positive:  FLG, THY1, SLIT3, ITIH5, PTGIS, KRT23, SPARC, CD14, CYTL1, CRB2 
    ##     IGFBP4, C6orf99, PRRX2, FGF18, AP1M2, SEPT4, UPK3B, TMOD1, KRT19, BCAT1 
    ##     SNCA, NGFR, CFTR, MAF, BCAM, RGS5, CD74, AC005482.1, TM4SF4, UPK1B 
    ## Negative:  SERPINE2, MT2A, PLAU, GYPC, EMP3, VIM, CXCL14, MT1E, PHLDA1, GMDS 
    ##     XKR4, TRIM55, CLDN16, ANPEP, PKIB, GLRX, TPM1, TM4SF1, LGALS1, DCBLD2 
    ##     ADIRF, HTRA1, S100A6, SPON1, ITGA3, TNC, GCHFR, FGFR1, MSMP, TGFBI 
    ## PC_ 5 
    ## Positive:  LDHB, ATP5F1B, PRDX1, DBI, NNMT, CPVL, ACAT2, MYL12B, SLC25A5, NQO1 
    ##     FDPS, ENO1, PSMA3, ANXA4, FDFT1, PSMA7, NDUFA9, HSPE1, UCHL1, NME1 
    ##     ANXA2, HSPD1, MYL12A, CCT5, CCT2, GSTM3, FHL2, HSPA8, FH, GGH 
    ## Negative:  NEAT1, NDRG1, TOP2A, ASPM, UBE2C, MKI67, PLK1, TMEM132A, HSPG2, GTSE1 
    ##     CDCA8, KIF23, AURKB, CENPF, ANLN, HJURP, KNL1, BTG1, TMSB10, AURKA 
    ##     TPX2, CCNA2, DEPDC1, KIF2C, RAP2B, SOX4, LAMB3, KLF6, MMP14, DLGAP5

``` r
print(ElbowPlot(SeuratObject,ndims = 50) + geom_vline(xintercept = 25, col="red"))
```

![](3_refine_clustering_CK5_organoid_files/figure-gfm/re_clust-1.png)<!-- -->

``` r
## Cell clustering
SeuratObject <- FindNeighbors(SeuratObject, dims = 1:25)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
SeuratObject <- FindClusters(SeuratObject, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2291
    ## Number of edges: 78732
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8293
    ## Number of communities: 8
    ## Elapsed time: 0 seconds

``` r
## Agreement with previous clustering
table("initial"=SeuratObject$init_seurat_clusters,
      "final"=SeuratObject$seurat_clusters)
```

    ##        final
    ## initial   0   1   2   3   4   5   6   7
    ##       0 551   0   0   1   8   4   0   0
    ##       1   6   4   0 348   1   1   0   0
    ##       2   0 349   3   1   0   0   0   0
    ##       3   0   1 348   0   1   0   0   0
    ##       4   2   0   0   0 279   0   0   0
    ##       5   0   8   0   1   0 231   0   0
    ##       6   0   0   0   0   1   0  87   1
    ##       7   0   0   0   0   0   0   0   0
    ##       8   0   0   0   0   0   0   0  54
    ##       9   0   0   0   0   0   0   0   0

``` r
table("Assigned"=SeuratObject$init_assign,
      "final"=SeuratObject$seurat_clusters)
```

    ##                  final
    ## Assigned            0   1   2   3   4   5   6   7
    ##   DCT-like_1        0 349   3   1   0   0   0   0
    ##   DCT-like_2        0   1 348   0   1   0   0   0
    ##   Epith.:HighMT     0   0   0   0   0   0   0   0
    ##   Epith.:HighRibo   0   0   0   0   0   0   0   0
    ##   PEC-like_1        0   0   0   0   1   0  87   1
    ##   PEC-like_2        0   0   0   0   0   0   0  54
    ##   Prolif.TPC        0   8   0   1   0 231   0   0
    ##   TPC_1           551   0   0   1   8   4   0   0
    ##   TPC_2             6   4   0 348   1   1   0   0
    ##   TPC_3             2   0   0   0 279   0   0   0

``` r
## UMAP
SeuratObject <- RunUMAP(SeuratObject, dims = 1:25)
```

    ## 08:04:38 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 08:04:38 Read 2291 rows and found 25 numeric columns

    ## 08:04:38 Using Annoy for neighbor search, n_neighbors = 30

    ## 08:04:38 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 08:04:38 Writing NN index file to temp file /tmp/RtmpRbknEF/file12cf623a5fbf
    ## 08:04:38 Searching Annoy index using 1 thread, search_k = 3000
    ## 08:04:38 Annoy recall = 100%
    ## 08:04:39 Commencing smooth kNN distance calibration using 1 thread
    ## 08:04:40 Initializing from normalized Laplacian + noise
    ## 08:04:40 Commencing optimization for 500 epochs, with 90478 positive edges
    ## 08:04:44 Optimization finished

``` r
d1 <- DimPlot(SeuratObject, group.by = "init_assign") + ggtitle("Initial Cell assignment")
d2 <- DimPlot(SeuratObject) + ggtitle("Re-classification")

print(CombinePlots(list(d1,d2)))
```

![](3_refine_clustering_CK5_organoid_files/figure-gfm/re_clust-2.png)<!-- -->

## Diagnostics of unsupervised clustering

We will investigate how the selection of multiple resolutions affects
the partition into individual cell
clusters.

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.")
```

![](3_refine_clustering_CK5_organoid_files/figure-gfm/clustree-1.png)<!-- -->

We also checked whether mithochondrial expression drives the clustering

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.", 
     node_colour="percent.mt", node_colour_aggr="mean")
```

![](3_refine_clustering_CK5_organoid_files/figure-gfm/clustree_mt-1.png)<!-- -->

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.", 
     node_colour="nFeature_RNA", node_colour_aggr="median")
```

![](3_refine_clustering_CK5_organoid_files/figure-gfm/clustree_nFeat-1.png)<!-- -->

``` r
VlnPlot(SeuratObject, 
    feature=c("percent.mt", "nFeature_RNA", 
        "Dissociation", "G2M.Score"), 
    pt.size=0.4,
    ncol=2)
```

![](3_refine_clustering_CK5_organoid_files/figure-gfm/vln_badQual_mtNfeat-1.png)<!-- -->

## Archive processed data for downstream analysis

``` r
DATA_DIR <- paste0(OUTDIR,"/data")
if(!dir.exists(DATA_DIR)) dir.create(DATA_DIR)
```

``` r
# 1 Clustering outcome
saveClusteringOutcome(SeuratObject , assay="RNA", fl=paste0(OUTDIR,"/init_clustering.tsv"))
```

    ## [WARN] Selected metacols for Clustering outcome:RNA_snn_res.0.1,RNA_snn_res.0.2,RNA_snn_res.0.3,RNA_snn_res.0.4,RNA_snn_res.0.5,RNA_snn_res.0.6,RNA_snn_res.0.7,RNA_snn_res.0.8,RNA_snn_res.0.9,RNA_snn_res.1,seurat_clusters

``` r
# 2 Initial idents (same as seurat_clusters)
saveActiveIdents(SeuratObject, fl=paste0(OUTDIR,"/active_idents.tsv"))
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
    ##  [1] grid      stats4    parallel  stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] cowplot_1.0.0        clustree_0.4.1       ggraph_2.0.0.9000   
    ##  [4] ComplexHeatmap_2.0.0 genesorteR_0.3.1     Matrix_1.2-17       
    ##  [7] dplyr_0.8.3          GSEABase_1.46.0      graph_1.62.0        
    ## [10] annotate_1.62.0      XML_3.98-1.20        AnnotationDbi_1.46.1
    ## [13] IRanges_2.18.2       S4Vectors_0.22.1     Biobase_2.44.0      
    ## [16] BiocGenerics_0.30.0  ggplot2_3.2.1        Seurat_3.1.0        
    ## [19] rmarkdown_1.15       nvimcom_0.9-82      
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] backports_1.1.4     circlize_0.4.7      plyr_1.8.4         
    ##   [4] igraph_1.2.4.1      lazyeval_0.2.2      splines_3.6.1      
    ##   [7] listenv_0.7.0       digest_0.6.21       htmltools_0.3.6    
    ##  [10] viridis_0.5.1       gdata_2.18.0        checkmate_1.9.4    
    ##  [13] magrittr_1.5        memoise_1.1.0       cluster_2.1.0      
    ##  [16] ROCR_1.0-7          globals_0.12.4      graphlayouts_0.5.0 
    ##  [19] RcppParallel_4.4.3  R.utils_2.9.0       colorspace_1.4-1   
    ##  [22] blob_1.2.0          ggrepel_0.8.1       xfun_0.9           
    ##  [25] crayon_1.3.4        RCurl_1.95-4.12     jsonlite_1.6       
    ##  [28] zeallot_0.1.0       survival_2.44-1.1   zoo_1.8-6          
    ##  [31] ape_5.3             glue_1.3.1          polyclip_1.10-0    
    ##  [34] gtable_0.3.0        leiden_0.3.1        GetoptLong_0.1.7   
    ##  [37] future.apply_1.3.0  shape_1.4.4         scales_1.0.0       
    ##  [40] pheatmap_1.0.12     DBI_1.0.0           bibtex_0.4.2       
    ##  [43] Rcpp_1.0.2          metap_1.1           viridisLite_0.3.0  
    ##  [46] xtable_1.8-4        clue_0.3-57         reticulate_1.13    
    ##  [49] bit_1.1-14          rsvd_1.0.2          mclust_5.4.5       
    ##  [52] SDMTools_1.1-221.1  tsne_0.1-3          htmlwidgets_1.3    
    ##  [55] httr_1.4.1          gplots_3.0.1.1      RColorBrewer_1.1-2 
    ##  [58] ica_1.0-2           pkgconfig_2.0.3     R.methodsS3_1.7.1  
    ##  [61] farver_1.1.0        uwot_0.1.4          tidyselect_0.2.5   
    ##  [64] labeling_0.3        rlang_0.4.0         reshape2_1.4.3     
    ##  [67] munsell_0.5.0       tools_3.6.1         RSQLite_2.1.2      
    ##  [70] ggridges_0.5.1      evaluate_0.14       stringr_1.4.0      
    ##  [73] yaml_2.2.0          npsurv_0.4-0        knitr_1.24         
    ##  [76] bit64_0.9-7         fitdistrplus_1.0-14 tidygraph_1.1.2    
    ##  [79] caTools_1.17.1.2    purrr_0.3.2         RANN_2.6.1         
    ##  [82] pbapply_1.4-2       future_1.14.0       nlme_3.1-141       
    ##  [85] R.oo_1.22.0         compiler_3.6.1      plotly_4.9.0       
    ##  [88] png_0.1-7           lsei_1.2-0          tibble_2.1.3       
    ##  [91] tweenr_1.0.1        stringi_1.4.3       RSpectra_0.15-0    
    ##  [94] lattice_0.20-38     vctrs_0.2.0         pillar_1.4.2       
    ##  [97] lifecycle_0.1.0     Rdpack_0.11-0       lmtest_0.9-37      
    ## [100] GlobalOptions_0.1.0 RcppAnnoy_0.0.13    data.table_1.12.8  
    ## [103] bitops_1.0-6        irlba_2.3.3         gbRd_0.4-11        
    ## [106] R6_2.4.0            KernSmooth_2.23-16  gridExtra_2.3      
    ## [109] codetools_0.2-16    MASS_7.3-51.4       gtools_3.8.1       
    ## [112] assertthat_0.2.1    rjson_0.2.20        withr_2.1.2        
    ## [115] sctransform_0.2.0   tidyr_1.0.0         Rtsne_0.15         
    ## [118] ggforce_0.3.1

``` r
{                                                                                                                                                                                                           
sink(file=paste0(OUTDIR,"/sessionInfo.txt"))
print(sessionInfo())
sink()
}
```
