PKD1-KO organoid derived from human CD24+ cells: Refine cell clustering
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
    ## Positive:  FTH1, KRT18, MT2A, GAPDH, SH3BGRL3, KRT19, S100A10, LCN2, C12orf75, IGFBP7 
    ##     S100A9, TUBA1B, KRT8, RARRES3, ANXA2, FTL, KRT7, H2AFZ, LDHA, MMP7 
    ##     SLC25A5, LGALS1, SLPI, HMGA1, RHOC, S100A4, FGFBP1, STMN1, PHLDA2, ANXA3 
    ## Negative:  XIST, NEAT1, POLR2J3.1, FTX, WSB1, ONECUT2, CCNL1, AC019197.1, KCNQ1OT1, GOLGA8A 
    ##     ZBTB20, NFAT5, LINC01320, RNF19A, SCN2A, GABPB1-AS1, ARHGEF38, TSIX, KIAA1551, GOLGA8B 
    ##     MUC16, SLC5A3, CEP95, CP, TNFAIP2, DGKH, HCG18, GLS, USP53, REL 
    ## PC_ 2 
    ## Positive:  AC073957.3, FTX, GOLGA8B, AL021368.2, AC019117.2, AC092683.1, TTN, GOLGA8A, POLR2J3.1, LINC00621 
    ##     HCG18, KCNQ1OT1, AP000542.2, TSIX, MIR29B2CHG, AC114760.2, ZNF891, AC003681.1, SCN2A, AL354733.3 
    ##     AF117829.1, ZMAT1, AL022323.4, LRRTM2, RIMKLB, DNM3OS, GABPB1-AS1, SPIN2A, XIST, MIR34AHG 
    ## Negative:  C4orf3, BNIP3L, PGK1, FAM162A, BTG1, GAPDH, HLA-C, TXNIP, IFITM1, VKORC1 
    ##     FTL, HLA-B, P4HA1, GATA3, EGLN3, BNIP3, JUNB, PYURF, YWHAH, MARCKS 
    ##     METRNL, MARCKSL1, SLPI, HLA-DRB1, NDRG1, ID2, CITED2, ID1, SERPINA1, FXYD3 
    ## PC_ 3 
    ## Positive:  VCAN, HLA-DMB, TMEM101, THSD7A, CXCL2, PKHD1, ABCC4, LINC02532, LINC01320, CYP27A1 
    ##     SOX6, RXFP1, IGFBP7, SERPING1, SCN2A, PROM1, FILIP1, SCTR, ADH1C, TNFRSF11B 
    ##     STIM2, TSPAN8, IFI44L, GABRP, SHKBP1, RNF144B, AIFM1, ADGRL2, COL4A1, DSEL 
    ## Negative:  ANGPTL4, NDRG1, P4HA1, HILPDA, ENO2, KISS1R, IGFBP5, EGLN3, PCSK6, VEGFA 
    ##     ADGRE2, FUT11, TCN1, FAM13A, ERO1A, SLC2A1, TMEM91, FABP5, GPI, TFRC 
    ##     AHNAK2, PTPRR, STRA6, VIM, KANK4, CA9, PAG1, PLOD2, SLC40A1, TXNIP 
    ## PC_ 4 
    ## Positive:  PLAT, PLAU, ECM1, MACC1, UPK1B, F3, PIM1, KRT23, AKAP12, RNF39 
    ##     PMEPA1, CYTOR, PLAUR, LOXL2, CAV2, F2RL1, IL1RL1, CLTRN, CDKN2B, LAMA3 
    ##     ITGB1, C15orf48, IL32, DCBLD2, RND3, SPAG9, CAV1, FAM83A, PTGES, SEMA3B 
    ## Negative:  SPP1, SLC40A1, SERPINB3, C4BPA, NOSTRIN, SLC34A2, LYPD1, SERPINA1, FAM20A, KANK4 
    ##     KCNQ2, PCSK1, TMEM101, PPP1R1B, MEIS1, CEACAM7, ANKS1B, ADAMTS1, CP, CEACAM5 
    ##     AGR3, MAP2K6, CD74, TESC, C2orf88, HABP2, LRG1, RARRES1, LCN2, CFTR 
    ## PC_ 5 
    ## Positive:  SPP1, VCAN, FXYD2, CYR61, S100A1, IGFBP7, COL4A2, CHST2, MYADM, COL4A1 
    ##     SLC3A1, SOD3, FBXO17, PROM1, ERBB4, CITED2, HIPK2, BMP2, FOXJ1, CREB5 
    ##     PTAFR, LINC02532, LINC01320, CYP27A1, FILIP1, MFAP2, SMAD7, HCFC1R1, IRX3, FAM84A 
    ## Negative:  CEACAM5, RARRES1, SPINK1, CEACAM7, MUC4, SLC40A1, KANK4, MEIS1, SERPINB3, STXBP6 
    ##     GPX2, PSCA, CEACAM6, CYP3A5, FAM83A, C15orf48, KLK11, PTGES, SPATA6L, PLEKHS1 
    ##     PCSK1, RAB27B, PARM1, SLC12A2, ASS1, ADH1C, FAM3D, ANKS1B, MAP2K6, FBP1

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
    ## Number of nodes: 496
    ## Number of edges: 20862
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7299
    ## Number of communities: 4
    ## Elapsed time: 0 seconds

``` r
## Agreement with previous clustering
table("initial"=SeuratObject$init_seurat_clusters,
      "final"=SeuratObject$seurat_clusters)
```

    ##        final
    ## initial   0   1   2   3
    ##       0 242   0  56   0
    ##       1   0   0   0   0
    ##       2   6 180   0   1
    ##       3   0   1   0  10

``` r
table("Assigned"=SeuratObject$init_assign,
      "final"=SeuratObject$seurat_clusters)
```

    ##               final
    ## Assigned         0   1   2   3
    ##   DCT-like_1   242   0  56   0
    ##   DCT-like_2     0   1   0  10
    ##   TPC_1:LowCov   0   0   0   0
    ##   TPC_2          6 180   0   1

``` r
## UMAP
SeuratObject <- RunUMAP(SeuratObject, dims = 1:25)
```

    ## 12:40:44 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 12:40:44 Read 496 rows and found 25 numeric columns

    ## 12:40:44 Using Annoy for neighbor search, n_neighbors = 30

    ## 12:40:44 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 12:40:44 Writing NN index file to temp file /tmp/RtmpXPT7MS/file154362744d60
    ## 12:40:44 Searching Annoy index using 1 thread, search_k = 3000
    ## 12:40:44 Annoy recall = 100%
    ## 12:40:44 Commencing smooth kNN distance calibration using 1 thread
    ## 12:40:45 Initializing from normalized Laplacian + noise
    ## 12:40:45 Commencing optimization for 500 epochs, with 18942 positive edges
    ## 12:40:46 Optimization finished

``` r
d1 <- DimPlot(SeuratObject, group.by = "init_assign") + ggtitle("Initial Cell assignment")
d2 <- DimPlot(SeuratObject) + ggtitle("Re-classification")

print(CombinePlots(list(d1,d2)))
```

![](3_refine_clustering_files/figure-gfm/re_clust-2.png)<!-- -->

## Diagnostics of unsupervised clustering

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
