PKD1-KO organoid derived from human CD24+ cells: initial clustering
================
Javier Perales-Patón - <javier.perales@bioquant.uni-heidelberg.de> -
ORCID: 0000-0003-0780-6683

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
SeuratObject <- getSeuratObject(path = "../data/sc/JX2_PKD1KO_organoid/filtered_feature_bc_matrix/",
                     project_name = "JX2_PKD1KO_organoid", 
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

![](1_initial_clustering_JX2_CD24_PKD1KO_files/figure-gfm/QC-1.png)<!-- -->

``` r
print(CombinePlots(plots = list(FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "percent.mt"),
                                FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "nFeature_RNA")))
)
```

![](1_initial_clustering_JX2_CD24_PKD1KO_files/figure-gfm/QC-2.png)<!-- -->

## Data cleasing: gene and cell filtering

To avoid unnecessary sparsity and noise in downstream analysis,  
\* We discard cells with low (\<200) and high (\>6000) number of genes
detected to avoid bad quality and doublets, repectively.  
\* We discard genes not detected across cells.  
For this we use cutoffs in concordance with kidney tissue, where
Proximal tubule cells might have a high content of mitochondrial genes
(\<20%). Later we will apply diagnostics to see if this is related.

``` r
nSamples_before <- ncol(SeuratObject)
SeuratObject <- subset(SeuratObject, nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
nSamples_after <- ncol(SeuratObject)
```

Doing so, 759 out of an initial 1329 total were retrieve for the
analysis, thus discarding a 570 cells.

## Pre-processing the data for cell clustering and cell-type assignment

First, we normalize the data using default parameters

``` r
SeuratObject <- NormalizeData(SeuratObject)
```

Second, we perform a feature selection based on `vst` method to select
high variable genes (top 2000) for cell clustering.

``` r
SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 2000)

print(LabelPoints(VariableFeaturePlot(SeuratObject),
                  points=VariableFeatures(SeuratObject)[1:10],repel = TRUE))
```

    ## Warning: Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
    ## Please use `as_label()` or `as_name()` instead.
    ## This warning is displayed once per session.

    ## When using repel, set xnudge and ynudge to 0 for optimal results

    ## Warning: Transformation introduced infinite values in continuous x-axis

![](1_initial_clustering_JX2_CD24_PKD1KO_files/figure-gfm/sel-1.png)<!-- -->

Third, we center and scale the data prior PCA:

``` r
SeuratObject <- Seurat_scaledata(SeuratObject)
```

    ## Centering and scaling data matrix

Finally we perform PCA, selecting a number of PCs that are relevant
based on the elbow plot. We tested that our results are not sensitive to
this parameter, and the arbitrary selection of number of PCs do not
change final conclusions.

``` r
SeuratObject <- RunPCA(SeuratObject, features = VariableFeatures(SeuratObject), npcs = 50)
```

    ## PC_ 1 
    ## Positive:  FTH1, MT2A, KRT18, KRT19, GAPDH, SH3BGRL3, LCN2, C12orf75, S100A9, FGFBP1 
    ##     IGFBP7, TUBA1B, KRT8, STMN1, RARRES3, LGALS1, ANXA3, S100A4, PCLAF, FTL 
    ##     LDHA, SLC25A5, HMGA1, ANXA2, PHLDA2, CXCL5, SLPI, KRT7, H2AFZ, C19orf33 
    ## Negative:  XIST, NEAT1, ONECUT2, WSB1, CCNL1, RNF19A, ZBTB20, NFAT5, POLR2J3.1, FTX 
    ##     AC019197.1, LINC01320, KIAA1551, DST, MUC16, CP, SCN2A, KCNQ1OT1, GABPB1-AS1, SLC5A3 
    ##     GOLGA8A, ARHGEF38, VEGFA, RAB11FIP1, TNFAIP2, GLS, CFTR, KLF6, COBLL1, SPAG9 
    ## PC_ 2 
    ## Positive:  GAPDH, MMP7, SLPI, LDHA, FTL, RHOC, PDZK1IP1, SH3BGRL3, KRT19, PGK1 
    ##     HLA-C, C4orf3, HLA-B, ANXA2, H2AFZ, MARCKSL1, KRT7, FTH1, IFITM1, MYDGF 
    ##     CSTB, HSPA8, CLDN7, SERPINA1, CD99, FAM162A, VDAC2, LDHB, KLK6, SFN 
    ## Negative:  FTX, XIST, GOLGA8B, GOLGA8A, KCNQ1OT1, TSIX, POLR2J3.1, NEAT1, AC092683.1, AC073957.3 
    ##     AL021368.2, TTN, GABPB1-AS1, AP000542.2, HCG18, SCN2A, WSB1, AC019117.2, AC003681.1, AL354733.3 
    ##     MIR29B2CHG, TMEM67, CBLB, MIR34AHG, AC019197.1, LINC00621, ARHGEF38, RIMKLB, AC114760.2, AF117829.1 
    ## PC_ 3 
    ## Positive:  ANGPTL4, NDRG1, P4HA1, EGLN3, KISS1R, HILPDA, TCN1, ENO2, PCSK6, TXNIP 
    ##     IGFBP5, SLC40A1, GPI, TMEM91, FUT11, AC097534.2, VEGFA, FABP5, ADGRE2, GATA3 
    ##     FAM162A, SLC2A1, PLOD2, KANK4, FAM13A, STRA6, TFRC, PPP1R1B, CA9, ID1 
    ## Negative:  VCAN, CCND1, CXCL2, SCN2A, IGFBP7, PKHD1, FILIP1, CPLANE1, THSD7A, SOX6 
    ##     HLA-DMB, HELLS, LINC01320, RXFP1, C12orf75, AF117829.1, ADI1, RNF19A, ABCC4, CXCL8 
    ##     TNFRSF11B, MAL, AC019197.1, KRT7, STIM2, PLAU, CYP27A1, MAPK10, HMGA1, RPGR 
    ## PC_ 4 
    ## Positive:  PLAT, MACC1, PLAU, KRT23, F3, AKAP12, ERO1A, DCBLD2, ITGB1, UPK1B 
    ##     ECM1, CYTOR, CAV2, PLAUR, LOXL2, PIM1, RND3, MAL2, PMEPA1, IL32 
    ##     SPAG9, PADI1, VEGFA, SLC7A5, AL133453.1, CDKN2B, ENO2, CAV1, IGFBP5, SLC5A3 
    ## Negative:  SPP1, TMEM101, CFTR, HLA-DMB, LCN2, CD74, C4BPA, SLC34A2, SERPINA1, TESC 
    ##     C2CD4A, PSAT1, AGR3, FAM20A, TNFSF10, LYPD1, NCOA7, LRG1, CP, AIFM1 
    ##     CXCL2, SERPING1, FAM181B, JUNB, HLA-DRA, ADAMTS1, SCTR, IFITM1, SOD2, HABP2 
    ## PC_ 5 
    ## Positive:  CEACAM5, RARRES1, SLC40A1, CEACAM7, SPINK1, MEIS1, KANK4, SERPINB3, MUC4, STXBP6 
    ##     GPX2, CYP3A5, PCSK1, PSCA, PLEKHS1, ANKS1B, CA8, MAP2K6, SLC12A2, C15orf48 
    ##     LCN2, KLK11, CEACAM6, SPATA6L, C2orf88, HP, ADH1C, PTGES, PARM1, FAM83A 
    ## Negative:  VCAN, SPP1, FXYD2, CYR61, IGFBP7, COL4A2, COL4A1, MYADM, SOD3, CHST2 
    ##     ERBB4, S100A1, HIPK2, PROM1, SLC3A1, BMP2, CITED2, FBXO17, CREB5, DUSP6 
    ##     TRIB1, MYLK, CYP27A1, FILIP1, VIM, PTAFR, TRIM55, SMAD7, IRX3, TINAG

``` r
print(ElbowPlot(SeuratObject,ndims = 50) + geom_vline(xintercept = 25, col="red"))
```

![](1_initial_clustering_JX2_CD24_PKD1KO_files/figure-gfm/pca-1.png)<!-- -->

We decided to use 25 PCs for the clustering.

We also calculate a dissociation score

``` r
dissociation <- read.csv("../data/Prior/scrnaseq-digestion-paper_coregene_df-FALSE-v3.csv",
             stringsAsFactors=FALSE)
dissTop40 <- intersect(head(dissociation$gene_symbol, 40), rownames(SeuratObject))
SeuratObject <- AddModuleScore(SeuratObject, features = list("dissociation"=dissTop40))
colnames(SeuratObject@meta.data)[which(colnames(SeuratObject@meta.data)=="Cluster1")] <- "Dissociation"
```

And cellcycle phase

``` r
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
SeuratObject <- CellCycleScoring(SeuratObject,
                 s.features = s.genes,
                 g2m.features = g2m.genes)
```

## Cell clustering

We are going to use Shared-Nearest Neighbour with Graph partitioning for
cell clustering. We run the algorithm with multiple resolutions. Later
this is read-out to give an idea how these are consistent from a lower
to high number of expected cell populations in the sample.

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
    ## Number of nodes: 759
    ## Number of edges: 35999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9276
    ## Number of communities: 2
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 759
    ## Number of edges: 35999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8740
    ## Number of communities: 2
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 759
    ## Number of edges: 35999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8203
    ## Number of communities: 2
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 759
    ## Number of edges: 35999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7819
    ## Number of communities: 3
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 759
    ## Number of edges: 35999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7465
    ## Number of communities: 4
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 759
    ## Number of edges: 35999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7130
    ## Number of communities: 5
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 759
    ## Number of edges: 35999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.6823
    ## Number of communities: 5
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 759
    ## Number of edges: 35999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.6528
    ## Number of communities: 6
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 759
    ## Number of edges: 35999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.6242
    ## Number of communities: 7
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 759
    ## Number of edges: 35999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.5978
    ## Number of communities: 7
    ## Elapsed time: 0 seconds

``` r
# We redefine the final partitioning with resolution 0.5
SeuratObject$seurat_clusters <- SeuratObject$RNA_snn_res.0.5 
Idents(SeuratObject) <- SeuratObject$RNA_snn_res.0.5 
```

## Diagnostics of unsupervised clustering

We will investigate how the selection of multiple resolutions affects
the partition into individual cell clusters.

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.")
```

![](1_initial_clustering_JX2_CD24_PKD1KO_files/figure-gfm/clustree-1.png)<!-- -->

We also checked whether mithochondrial expression drives the clustering

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.", 
     node_colour="percent.mt", node_colour_aggr="mean")
```

![](1_initial_clustering_JX2_CD24_PKD1KO_files/figure-gfm/clustree_mt-1.png)<!-- -->

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.", 
     node_colour="nFeature_RNA", node_colour_aggr="median")
```

![](1_initial_clustering_JX2_CD24_PKD1KO_files/figure-gfm/clustree_nFeat-1.png)<!-- -->

``` r
VlnPlot(SeuratObject, 
    feature=c("percent.mt", "nFeature_RNA", 
        "Dissociation", "G2M.Score"), 
    pt.size=0.4,
    ncol=2)
```

![](1_initial_clustering_JX2_CD24_PKD1KO_files/figure-gfm/vln_badQual_mtNfeat-1.png)<!-- -->

## Clustering decision

We decided to chose resolution 0.5 for the initial clustering. We will
save this output to archive the outcome. Herein, first column is the
cell barcode, and rest of columns are clusters where each cell belong to
across multiple resolutions. The ones that are active for the study are
those related to resolution 0.5 (`RNA_snn_res.0.5` and
`seurat_clusters`).

``` r
# 1 Clustering outcome
saveClusteringOutcome(SeuratObject , assay="RNA", fl=paste0(OUTDIR,"/init_clustering.tsv"))
```

    ## [WARN] Selected metacols for Clustering outcome:RNA_snn_res.0.1,RNA_snn_res.0.2,RNA_snn_res.0.3,RNA_snn_res.0.4,RNA_snn_res.0.5,RNA_snn_res.0.6,RNA_snn_res.0.7,RNA_snn_res.0.8,RNA_snn_res.0.9,RNA_snn_res.1,seurat_clusters

``` r
# 2 Initial idents (same as seurat_clusters)
saveActiveIdents(SeuratObject, fl=paste0(OUTDIR,"/active_idents.tsv"))
```

## Non-linear dim reduction (umap)

``` r
SeuratObject <- RunUMAP(SeuratObject, dims = 1:25)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 06:53:28 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 06:53:28 Read 759 rows and found 25 numeric columns

    ## 06:53:28 Using Annoy for neighbor search, n_neighbors = 30

    ## 06:53:28 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 06:53:28 Writing NN index file to temp file /tmp/RtmpocdQIE/file74b0761f8055
    ## 06:53:28 Searching Annoy index using 1 thread, search_k = 3000
    ## 06:53:28 Annoy recall = 100%
    ## 06:53:29 Commencing smooth kNN distance calibration using 1 thread
    ## 06:53:29 Initializing from normalized Laplacian + noise
    ## 06:53:29 Commencing optimization for 500 epochs, with 31948 positive edges
    ## 06:53:32 Optimization finished

``` r
DimPlot(SeuratObject, reduction="umap") + ggtitle(Project(SeuratObject))
```

![](1_initial_clustering_JX2_CD24_PKD1KO_files/figure-gfm/umap-1.png)<!-- -->

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
    ##  [1] grid      stats4    parallel  stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] rmarkdown_1.15       openxlsx_4.2.3       cowplot_1.0.0       
    ##  [4] clustree_0.4.1       ggraph_2.0.0.9000    ComplexHeatmap_2.0.0
    ##  [7] genesorteR_0.3.1     Matrix_1.2-17        dplyr_0.8.3         
    ## [10] GSEABase_1.46.0      graph_1.62.0         annotate_1.62.0     
    ## [13] XML_3.98-1.20        AnnotationDbi_1.46.1 IRanges_2.18.2      
    ## [16] S4Vectors_0.22.1     Biobase_2.44.0       BiocGenerics_0.30.0 
    ## [19] ggplot2_3.2.1        Seurat_3.1.0         nvimcom_0.9-82      
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
    ##  [61] farver_1.1.0        uwot_0.1.4          labeling_0.3       
    ##  [64] tidyselect_0.2.5    rlang_0.4.0         reshape2_1.4.3     
    ##  [67] munsell_0.5.0       tools_3.6.1         RSQLite_2.1.2      
    ##  [70] ggridges_0.5.1      evaluate_0.14       stringr_1.4.0      
    ##  [73] yaml_2.2.0          npsurv_0.4-0        knitr_1.24         
    ##  [76] bit64_0.9-7         fitdistrplus_1.0-14 tidygraph_1.1.2    
    ##  [79] zip_2.1.1           caTools_1.17.1.2    purrr_0.3.2        
    ##  [82] RANN_2.6.1          pbapply_1.4-2       future_1.14.0      
    ##  [85] nlme_3.1-141        R.oo_1.22.0         compiler_3.6.1     
    ##  [88] plotly_4.9.0        png_0.1-7           lsei_1.2-0         
    ##  [91] tibble_2.1.3        tweenr_1.0.1        stringi_1.4.3      
    ##  [94] RSpectra_0.15-0     lattice_0.20-38     vctrs_0.2.0        
    ##  [97] pillar_1.4.2        lifecycle_0.1.0     Rdpack_0.11-0      
    ## [100] lmtest_0.9-37       GlobalOptions_0.1.0 RcppAnnoy_0.0.13   
    ## [103] data.table_1.12.8   bitops_1.0-6        irlba_2.3.3        
    ## [106] gbRd_0.4-11         R6_2.4.0            KernSmooth_2.23-16 
    ## [109] gridExtra_2.3       codetools_0.2-16    MASS_7.3-51.4      
    ## [112] gtools_3.8.1        assertthat_0.2.1    rjson_0.2.20       
    ## [115] withr_2.1.2         sctransform_0.2.0   tidyr_1.0.0        
    ## [118] Rtsne_0.15          ggforce_0.3.1

``` r
{                                                                                                                                                                                                           
sink(file=paste0(OUTDIR,"/sessionInfo.txt"))
print(sessionInfo())
sink()
}
```
