CK119 late healthy organoid : Refine cell clustering
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
                                      Idents(SeuratObject),
                      ignore.case=TRUE)]
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
    ## Positive:  PTTG1, TYMS, CENPF, KIAA0101, UBE2S, TUBA1B, MKI67, KPNA2, ANLN, TPX2 
    ##     PRC1, CENPM, MAD2L1, MT2A, CKS2, RRM2, CDKN3, PHF19, TOP2A, CEP55 
    ##     CKS1B, STMN1, KIF20B, RANBP1, ZWINT, DTYMK, NUSAP1, UBE2T, PBK, HMMR 
    ## Negative:  SAT1, LCN2, SLPI, C15orf48, CLDN7, KLK6, S100A9, SAA2, ERO1A, PERP 
    ##     CD55, RAB11FIP1, RNASET2, MUC1, ISG20, GCNT3, KLK8, SDCBP2, PTGES, RAB25 
    ##     ELF3, PRSS8, FXYD3, MMP7, KLK7, SMIM5, TSPAN1, TACSTD2, SPINT1, HLA-B 
    ## PC_ 2 
    ## Positive:  IGFBP7, CLU, SPP1, VIM, VCAN, ANXA4, CRYM, PPP1R1A, HTRA1, NNMT 
    ##     GYPC, SLC34A2, SPON1, FXYD2, CYS1, HLA-DMB, CPVL, STC1, EMP3, AKR1C3 
    ##     TUBA1A, NTRK2, BNIP3, CTSC, FILIP1L, FAM162A, KCNJ16, ADAMTS1, FAM134B, TNFSF10 
    ## Negative:  MAL2, CLDN4, AREG, KRT7, CLDN7, KRT18, TACSTD2, SFN, KLK6, EZR 
    ##     HN1, S100A14, TPM4, GPRC5A, ITGA2, CST6, KRT8, GIPC1, RAB11FIP1, HES4 
    ##     PTGES, RAB25, GALNT3, EHF, MACC1, RTKN2, GCNT3, SMIM5, EREG, C6orf132 
    ## PC_ 3 
    ## Positive:  PLAU, AKAP12, CCND1, DCBLD2, S100A16, TPM1, OCIAD2, TUBA1A, CAV2, RHOF 
    ##     MAL, CLDN1, LAMC2, GLRX, PHLDA2, TNFRSF12A, HMGA1, EIF6, WNT7A, SCEL 
    ##     ANXA3, ADIRF, ARL4C, CAP1, EZR, SH3BGRL3, S100A6, RAB32, ITGB8, RHOD 
    ## Negative:  WFDC2, OLFM4, RARRES2, CP, LTF, BCAM, CLU, CFI, RARRES1, DEFB1 
    ##     C3, ATP1B1, MNS1, IGFBP3, SPP1, SERPINA1, PIGR, TCN1, AGR2, BBOX1 
    ##     ST8SIA4, HMGB2, CFTR, SLPI, TOP2A, PPP1R1B, TMEM91, MKI67, CA12, NUSAP1 
    ## PC_ 4 
    ## Positive:  SQSTM1, TMEM140, RSAD2, FTL, GDF15, IFIT3, PIK3IP1, FTH1, MX1, BTG1 
    ##     IFIT1, CDKN1A, ISG15, HERC5, IFIT2, HLA-B, ATF3, ISG20, OAS1, BBC3 
    ##     PSMB9, OPTN, DDIT3, HIST1H2AC, H1FX, IFI6, THSD7A, YPEL3, CMPK2, XAF1 
    ## Negative:  KRT8, KRT18, FGFBP1, KITLG, ARHGAP29, MALT1, EMP1, S100A6, S100A14, ACTB 
    ##     TNFRSF11B, SPP1, AREG, BCAM, TSPAN1, S100A10, HIGD1A, GPRC5A, TRNP1, AGR2 
    ##     MALL, MMP7, ANXA2, F3, KRT19, PHLDA2, MPZL2, MUC1, PROM2, LIMA1 
    ## PC_ 5 
    ## Positive:  NEAT1, ARHGAP29, ATP1B1, TXNIP, BTG1, IER3, PNRC1, TNFSF10, UBC, JUND 
    ##     ELF3, CEBPB, SOX4, TSC22D1, DCDC2, ZFP36L2, IER2, ANXA4, TNFRSF11B, UGCG 
    ##     JUN, CEBPD, ARID5B, KRT19, CITED4, IVNS1ABP, DDIT4, WSB1, ACSL4, NNMT 
    ## Negative:  SH3BGRL3, CD59, KLK10, HLA-G, ECM1, S100A10, CD99, S100A6, IL32, CST6 
    ##     MSLN, FXYD5, C12orf75, UPK1B, LAMC2, KLK6, KLK11, C15orf48, CEACAM6, S100P 
    ##     C9orf16, B2M, PERP, ADIRF, ACTB, MFI2, ALDH1A3, HMGA1, SPRR2A, MFGE8

``` r
print(ElbowPlot(SeuratObject,ndims = 50) + geom_vline(xintercept = 25, col="red"))
```

![](3_refine_clustering_files/figure-gfm/re_clust-1.png)<!-- -->

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
    ## Number of nodes: 3693
    ## Number of edges: 128997
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8812
    ## Number of communities: 8
    ## Elapsed time: 0 seconds

``` r
## Agreement with previous clustering
table("initial"=SeuratObject$init_seurat_clusters,
      "final"=SeuratObject$seurat_clusters)
```

    ##        final
    ## initial   0   1   2   3   4   5   6   7
    ##       0 820   0   3   9   0   0   2   0
    ##       1   3 617   0   0   0   1   0   0
    ##       2   2   0 526   2   1   0   0   0
    ##       3   0   0   1   9 480   0   1   0
    ##       4   2   0   1 471   3   2   0   0
    ##       5   0   1   0   0   0 404   0   0
    ##       6   1   0   0   0   1   0 253   0
    ##       7   0   0   0   0   0   0   0   0
    ##       8   0   0   0   0   0   0   0   0
    ##       9   0   0   0   0   0   0   1  76

``` r
table("Assigned"=SeuratObject$init_assign,
      "final"=SeuratObject$seurat_clusters)
```

    ##                       final
    ## Assigned                 0   1   2   3   4   5   6   7
    ##   Epith.Unknown          1   0   0   0   1   0 253   0
    ##   PEC-like               0   0   0   0   0   0   1  76
    ##   Prolif.TPC_1           3 617   0   0   0   1   0   0
    ##   Prolif.TPC_2           0   1   0   0   0 404   0   0
    ##   Prolif.TPC_3:LowQual   0   0   0   0   0   0   0   0
    ##   TPC_1                820   0   3   9   0   0   2   0
    ##   TPC_2                  2   0 526   2   1   0   0   0
    ##   TPC_3                  0   0   1   9 480   0   1   0
    ##   TPC_4                  2   0   1 471   3   2   0   0
    ##   TPC_5:LowQual          0   0   0   0   0   0   0   0

``` r
## UMAP
SeuratObject <- RunUMAP(SeuratObject, dims = 1:25)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 11:56:18 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 11:56:18 Read 3693 rows and found 25 numeric columns

    ## 11:56:18 Using Annoy for neighbor search, n_neighbors = 30

    ## 11:56:18 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 11:56:19 Writing NN index file to temp file /tmp/RtmpG0v0fe/file52b4dd9486f
    ## 11:56:19 Searching Annoy index using 1 thread, search_k = 3000
    ## 11:56:20 Annoy recall = 100%
    ## 11:56:20 Commencing smooth kNN distance calibration using 1 thread
    ## 11:56:21 Initializing from normalized Laplacian + noise
    ## 11:56:21 Commencing optimization for 500 epochs, with 151124 positive edges
    ## 11:56:30 Optimization finished

``` r
d1 <- DimPlot(SeuratObject, group.by = "init_assign") + ggtitle("Initial Cell assignment")
d2 <- DimPlot(SeuratObject) + ggtitle("Re-classification")

print(CombinePlots(list(d1,d2)))
```

![](3_refine_clustering_files/figure-gfm/re_clust-2.png)<!-- -->

## Diagnostics of unsupervised clustering

We will investigate how the selection of multiple resolutions affects
the partition into individual cell clusters.

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.")
```

![](3_refine_clustering_files/figure-gfm/clustree-1.png)<!-- -->

We also checked whether mithochondrial expression drives the clustering

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.", 
     node_colour="percent.mt", node_colour_aggr="mean")
```

![](3_refine_clustering_files/figure-gfm/clustree_mt-1.png)<!-- -->

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.", 
     node_colour="nFeature_RNA", node_colour_aggr="median")
```

![](3_refine_clustering_files/figure-gfm/clustree_nFeat-1.png)<!-- -->

``` r
VlnPlot(SeuratObject, 
    feature=c("percent.mt", "nFeature_RNA", 
        "Dissociation", "G2M.Score"), 
    pt.size=0.4,
    ncol=2)
```

![](3_refine_clustering_files/figure-gfm/vln_badQual_mtNfeat-1.png)<!-- -->

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
