Control2 whole tissue biopsy, healthy: initial clustering
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
SeuratObject <- getSeuratObject(path = "../data/sc/Control2/filtered_feature_bc_matrix/",
                     project_name = "Control2_kidney", 
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

![](1_initial_clustering_files/figure-gfm/QC-1.png)<!-- -->

``` r
print(CombinePlots(plots = list(FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "percent.mt"),
                                FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "nFeature_RNA")))
)
```

![](1_initial_clustering_files/figure-gfm/QC-2.png)<!-- -->

## Data cleasing: gene and cell filtering

To avoid unnecessary sparsity and noise in downstream analysis,  
\* We discard cells with low (\<200) and high (\>4000) number of genes
detected to avoid bad quality and doublets, repectively.  
\* We discard genes not detected across cells.  
In addition, we consider a thresholds of 10% for MT gene expression
since this is single-nucleus RNAseq but ambient RNA from cytoplasm has
contaminated the library. Later we will ignore remaining MT gene
expression because this cannot be derived from the nucleus.

``` r
nSamples_before <- ncol(SeuratObject)
SeuratObject <- subset(SeuratObject, nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
nSamples_after <- ncol(SeuratObject)
# Ignore MT genes
SeuratObject <- SeuratObject[grep("^MT-", rownames(SeuratObject), invert=TRUE),]
```

Doing so, 4735 out of an initial 5039 total were retrieve for the
analysis, thus discarding a 304 cells.

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

    ## Warning: Transformation introduced infinite values in continuous x-axis

![](1_initial_clustering_files/figure-gfm/sel-1.png)<!-- -->

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
    ## Positive:  ERBB4, MECOM, PRKG1, KCNIP4, EFNA5, NEBL, AC019197.1, GPC5, KAZN, NAV2 
    ##     ST6GALNAC3, RBFOX1, ATP6V0A4, PODXL, CA10, LINC01606, SLC12A1, CACNB2, NOS1AP, AC138305.1 
    ##     CASR, NTNG1, PLCB1, PDE1A, SGIP1, COBLL1, ATP10A, CA12, PTPRQ, CACNA2D3 
    ## Negative:  SLC13A3, ACSM2B, ACSM2A, LRP2, AC087762.1, SLC28A1, CUBN, SLC17A1, CNTNAP3B, SLC4A4 
    ##     BNC2, CDH6, NOX4, SLC5A12, SLC36A2, RASSF4, LINC01060, TINAG, SLC47A2, LINC00871 
    ##     GDA, AGXT2, CLDN10, SLC5A2, SLC47A1, FAM20A, DPYS, TRPM3, SMIM2-AS1, AC096577.1 
    ## PC_ 2 
    ## Positive:  PTPRQ, ST6GALNAC3, NPHS2, NPHS1, NTNG1, PTPRO, CLIC5, PODXL, ATP10A, PLA2R1 
    ##     FGF1, PLCE1, NFASC, ADAMTS19, FMN2, SPOCK2, PCOLCE2, TARID, TYRO3, ALS2CL 
    ##     DACH2, CTGF, WT1, HTRA1, AC109466.1, FYN, NPAS3, NDNF, SRGAP2, AC092813.2 
    ## Negative:  MECOM, CCSER1, KCNIP4, AC019197.1, COBLL1, EFNA5, CACNA2D3, CA12, MAGI1, KAZN 
    ##     LINC01606, ERBB4, PLCL1, CASR, AC138305.1, SLC12A1, TFCP2L1, GPC5, ATP6V0A4, ITPR1 
    ##     KCTD1, UMOD, PTGER3, ESRRB, EGF, AC005208.1, LINC01762, TMEM52B, PDE1A, ENOX1 
    ## PC_ 3 
    ## Positive:  PTPRQ, NTNG1, PTPRO, NPHS1, NPHS2, CLIC5, FGF1, ATP10A, PLA2R1, ADAMTS19 
    ##     PLCE1, COL4A3, NEBL, FMN2, MAGI2, PARD3B, DPP6, NFASC, TARID, DACH2 
    ##     NPAS3, SPOCK2, DACH1, PCOLCE2, AC109466.1, CCBE1, COL4A4, RBFOX1, WT1, AC092813.2 
    ## Negative:  EMCN, LDB2, MEIS2, PTPRB, ENG, SLCO2A1, FLT1, TEK, RAPGEF4, ZEB1 
    ##     PECAM1, HEG1, SHANK3, IGFBP5, EPAS1, PLPP1, NAV1, EGFL7, HECW2, SH3RF3 
    ##     CACNA1C, TGFBR2, GNA14, PLAT, ELMO1, ERG, CEACAM1, EBF1, PRKCH, ID1 
    ## PC_ 4 
    ## Positive:  SLC12A1, CCSER1, GPC5, CACNA2D3, CASR, PLCB1, UMOD, KCNIP4, ENOX1, HIP1 
    ##     ESRRB, LINC01606, ADAMTS9-AS2, PCDH9, SIM2, LINC01762, HS6ST2, AC092078.2, AC005208.1, EGF 
    ##     ERBB4, KNG1, AC058822.1, SGIP1, LINC00970, PHACTR1, GP2, ACPP, CLDN16, CPM 
    ## Negative:  CLNK, PDE1C, SLC26A7, IL18, LINC01187, ATP6V1C2, CA8, CELF2, NXPH2, HS6ST3 
    ##     LEF1, NRXN3, AC092422.1, ADGRF5, GALNT17, CNTNAP5, LAMA2, STAP1, KIT, VWA5B1 
    ##     PSD3, AC026333.3, ADAMTSL1, RCAN2, ADGRF1, PACRG, SNTB1, SYT17, AQP6, ATP6V0D2 
    ## PC_ 5 
    ## Positive:  SLC12A1, PLCB1, UMOD, ENOX1, SIM2, CASR, GPC5, CACNB2, CLDN10-AS1, PDE1A 
    ##     SLC16A12, LINC01606, GP2, AC058822.1, ACPP, PPM1E, CLDN16, RBFOX1, HIP1, CACNA2D3 
    ##     SCHLAP1, CABP1, CLDN14, RP1, SGIP1, THSD4, SCD5, PHACTR1, INPP4B, ESRRB 
    ## Negative:  PIK3C2G, SLC8A1, LINC01099, TEX41, LINC01098, LHX1, BMPR1B, SLC8A1-AS1, GATA3, GRIP1 
    ##     SCNN1G, TRPM6, SCNN1B, SOX5, TMTC2, MYO1B, SLC12A3, AC019197.1, TBC1D9, CADPS2 
    ##     SCN2A, ATF3, TOX3, AF233439.1, KRT19, PWRN1, SNTG1, SCN3A, FOSB, MID1

``` r
print(ElbowPlot(SeuratObject,ndims = 50) + geom_vline(xintercept = 25, col="red"))
```

![](1_initial_clustering_files/figure-gfm/pca-1.png)<!-- -->

We decided to use 25 PCs for the clustering.

We also calculate a dissociation
score

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
    ## Number of nodes: 4735
    ## Number of edges: 182364
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9737
    ## Number of communities: 10
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 4735
    ## Number of edges: 182364
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9560
    ## Number of communities: 11
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 4735
    ## Number of edges: 182364
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9415
    ## Number of communities: 13
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 4735
    ## Number of edges: 182364
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9292
    ## Number of communities: 14
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 4735
    ## Number of edges: 182364
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9182
    ## Number of communities: 15
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 4735
    ## Number of edges: 182364
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9092
    ## Number of communities: 16
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 4735
    ## Number of edges: 182364
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9004
    ## Number of communities: 16
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 4735
    ## Number of edges: 182364
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8921
    ## Number of communities: 16
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 4735
    ## Number of edges: 182364
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8839
    ## Number of communities: 18
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 4735
    ## Number of edges: 182364
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8757
    ## Number of communities: 19
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

![](1_initial_clustering_files/figure-gfm/clustree-1.png)<!-- -->

We also checked whether mithochondrial expression drives the clustering

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.", 
     node_colour="percent.mt", node_colour_aggr="mean")
```

![](1_initial_clustering_files/figure-gfm/clustree_mt-1.png)<!-- -->

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.", 
     node_colour="nFeature_RNA", node_colour_aggr="median")
```

![](1_initial_clustering_files/figure-gfm/clustree_nFeat-1.png)<!-- -->

``` r
VlnPlot(SeuratObject, 
    feature=c("percent.mt", "nFeature_RNA", 
        "Dissociation", "G2M.Score"), 
    pt.size=0.4,
    ncol=2)
```

![](1_initial_clustering_files/figure-gfm/vln_badQual_mtNfeat-1.png)<!-- -->

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

    ## 07:38:00 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 07:38:00 Read 4735 rows and found 25 numeric columns

    ## 07:38:00 Using Annoy for neighbor search, n_neighbors = 30

    ## 07:38:00 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 07:38:01 Writing NN index file to temp file /tmp/RtmpyF69L0/file3e9d2e5946c9
    ## 07:38:01 Searching Annoy index using 1 thread, search_k = 3000
    ## 07:38:03 Annoy recall = 100%
    ## 07:38:03 Commencing smooth kNN distance calibration using 1 thread
    ## 07:38:03 Initializing from normalized Laplacian + noise
    ## 07:38:04 Commencing optimization for 500 epochs, with 203726 positive edges
    ## 07:38:15 Optimization finished

``` r
DimPlot(SeuratObject, reduction="umap") + ggtitle(Project(SeuratObject))
```

![](1_initial_clustering_files/figure-gfm/umap-1.png)<!-- -->

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
    ## [5] Seurat_3.1.0      rmarkdown_1.15    nvimcom_0.9-82   
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
    ##  [31] R.utils_2.9.0       pkgconfig_2.0.3     htmltools_0.3.6    
    ##  [34] bibtex_0.4.2        htmlwidgets_1.3     rlang_0.4.0        
    ##  [37] farver_1.1.0        zoo_1.8-6           jsonlite_1.6       
    ##  [40] ica_1.0-2           gtools_3.8.1        dplyr_0.8.3        
    ##  [43] R.oo_1.22.0         magrittr_1.5        Matrix_1.2-17      
    ##  [46] Rcpp_1.0.2          munsell_0.5.0       viridis_0.5.1      
    ##  [49] ape_5.3             reticulate_1.13     lifecycle_0.1.0    
    ##  [52] R.methodsS3_1.7.1   stringi_1.4.3       yaml_2.2.0         
    ##  [55] gbRd_0.4-11         MASS_7.3-51.4       gplots_3.0.1.1     
    ##  [58] Rtsne_0.15          plyr_1.8.4          grid_3.6.1         
    ##  [61] parallel_3.6.1      gdata_2.18.0        listenv_0.7.0      
    ##  [64] ggrepel_0.8.1       crayon_1.3.4        lattice_0.20-38    
    ##  [67] graphlayouts_0.5.0  splines_3.6.1       SDMTools_1.1-221.1 
    ##  [70] zeallot_0.1.0       knitr_1.24          pillar_1.4.2       
    ##  [73] igraph_1.2.4.1      future.apply_1.3.0  reshape2_1.4.3     
    ##  [76] codetools_0.2-16    leiden_0.3.1        glue_1.3.1         
    ##  [79] evaluate_0.14       lsei_1.2-0          metap_1.1          
    ##  [82] RcppParallel_4.4.3  data.table_1.12.8   tweenr_1.0.1       
    ##  [85] vctrs_0.2.0         png_0.1-7           Rdpack_0.11-0      
    ##  [88] polyclip_1.10-0     gtable_0.3.0        RANN_2.6.1         
    ##  [91] purrr_0.3.2         tidyr_1.0.0         future_1.14.0      
    ##  [94] assertthat_0.2.1    ggforce_0.3.1       xfun_0.9           
    ##  [97] rsvd_1.0.2          tidygraph_1.1.2     RSpectra_0.15-0    
    ## [100] survival_2.44-1.1   viridisLite_0.3.0   tibble_2.1.3       
    ## [103] cluster_2.1.0       globals_0.12.4      fitdistrplus_1.0-14
    ## [106] ROCR_1.0-7

``` r
{                                                                                                                                                                                                           
sink(file=paste0(OUTDIR,"/sessionInfo.txt"))
print(sessionInfo())
sink()
}
```
