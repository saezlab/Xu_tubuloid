Organoid derived from CD24+ cells, using H.Clevers protocol: initial
clustering
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
SeuratObject <- getSeuratObject(path = "../data/sc/JX1_HC_organoid/filtered_feature_bc_matrix/",
                     project_name = "JX1_HC_organoid", 
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

![](1_initial_clustering_JX1_HC_organoid_files/figure-gfm/QC-1.png)<!-- -->

``` r
print(CombinePlots(plots = list(FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "percent.mt"),
                                FeatureScatter(SeuratObject, 
                           feature1 = "nCount_RNA", 
                           feature2 = "nFeature_RNA")))
)
```

![](1_initial_clustering_JX1_HC_organoid_files/figure-gfm/QC-2.png)<!-- -->

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

Doing so, 7802 out of an initial 8779 total were retrieve for the
analysis, thus discarding a 977 cells.

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

![](1_initial_clustering_JX1_HC_organoid_files/figure-gfm/sel-1.png)<!-- -->

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
    ## Positive:  MT2A, ACTB, KRT8, TUBA1B, LGALS1, CCND1, HMGA1, COTL1, ARL4C, DCBLD2 
    ##     EIF5A, PHLDA2, CTNNAL1, PCLAF, KRT18, TNFRSF12A, C12orf75, PTTG1, HSP90AA1, SFRP1 
    ##     ACAT2, S100A6, TUBB, S100A16, PHLDA1, JPT1, FDPS, CDKN3, SQLE, H2AFZ 
    ## Negative:  WFDC2, CLU, SPP1, IGFBP3, MALAT1, TXNIP, NEAT1, CP, SLPI, AGR3 
    ##     PIGR, CFTR, NTRK2, ADAMTS1, CD74, LTF, GBP2, NDRG1, GPNMB, C3 
    ##     SOD2, ERBB4, HPN, SAT1, SLC3A1, TMEM101, BTG2, SERPING1, TCIM, PIK3R1 
    ## PC_ 2 
    ## Positive:  S100A6, TMSB10, PLAU, HPGD, SH3BGRL3, LGALS1, SFRP1, FHL2, EMP3, DPP4 
    ##     SYTL2, ANXA2, ANPEP, VIM, S100A10, ADIRF, TUBA1A, FXYD2, FGFBP1, DCBLD2 
    ##     CST6, CDA, S100A16, ANXA1, TPM1, KRT19, PKIB, CCND1, ITGA3, PDZK1IP1 
    ## Negative:  CENPF, MKI67, ASPM, TOP2A, HMGB2, HMMR, DLGAP5, NUSAP1, UBE2C, CCNA2 
    ##     PBK, TPX2, NUF2, GTSE1, KIF23, SGO2, CENPA, CDC20, NCAPG, CDCA8 
    ##     CENPE, DEPDC1, TTK, NDC80, PLK1, SPC25, CEP55, KIF4A, ANLN, KIF2C 
    ## PC_ 3 
    ## Positive:  DBI, ATP1B1, TXNIP, TNFSF10, TNFRSF11B, DEFB1, MARCKSL1, KITLG, ZNF503, CFTR 
    ##     CEBPD, CITED4, PFN1, CAMK2N1, HOXA10, HIPK2, SOD2, CPVL, JUNB, ZFP36L2 
    ##     VCAN, ERBB4, CD74, UTRN, ACSL4, TRAM1, LINC02532, DUSP6, DEPTOR, AGR3 
    ## Negative:  MT2A, TUBA1B, KRT18, LGALS1, TUBB4B, SFN, ANXA2, HMGA1, FGFBP1, VIM 
    ##     ITGA3, CCNB1, ANGPTL4, S100A16, TMSB10, AREG, SERPINE2, KRT7, PLAU, CENPF 
    ##     KLK5, KLK6, EMP3, CDC20, CKS2, TOP2A, TUBA1C, MKI67, CAVIN3, SFRP1 
    ## PC_ 4 
    ## Positive:  KRT18, LCN2, KRT7, KRT19, KLK6, KRT8, AGR2, CFL1, S100A9, MMP7 
    ##     CPVL, PPIA, TNFRSF11B, ACTB, C12orf75, FDPS, SPP1, ANXA2, TUBA1B, TESC 
    ##     CLU, SOD3, DBI, TUBB4B, PHGDH, HNRNPA2B1, S100A16, WFDC2, CD74, HLA-DMB 
    ## Negative:  NDRG1, NEAT1, P4HA1, BTG1, HILPDA, MALAT1, VEGFA, ERO1A, IGFBP5, ENO2 
    ##     SLC2A1, RAP2B, SLC5A3, MARCKS, RND3, LOX, KLF6, CAV1, SNHG12, ARRDC3 
    ##     HK2, AHNAK2, PLAUR, CDKN1A, SOX4, ZFAS1, PMAIP1, TENT5A, ANKRD37, BHLHE40 
    ## PC_ 5 
    ## Positive:  ONECUT2, MALAT1, FAM111B, MEIS1, EFNB2, CP, C2orf88, KANK4, MCM10, CDC6 
    ##     SLC7A11, DTL, CLSPN, MCM3, HNRNPAB, MUC16, GPRC5A, CDCA7, MCM4, HELLS 
    ##     ATAD2, IGFBP5, UHRF1, HES4, GINS2, CADM1, KPNB1, MACC1, SLC1A5, KCNQ1OT1 
    ## Negative:  SPP1, IGFBP7, FXYD2, TIMP1, PDZK1IP1, VIM, CLU, DBI, LGALS1, PPIA 
    ##     HTRA1, ADIRF, CKB, PGK1, BNIP3, PPP1R1A, C12orf75, AKR1C3, TMSB10, DLGAP5 
    ##     NNMT, HMMR, TNFSF10, RARRES2, FTL, CENPF, SERPINA1, PLK1, NEK2, IGFBP6

``` r
print(ElbowPlot(SeuratObject,ndims = 50) + geom_vline(xintercept = 25, col="red"))
```

![](1_initial_clustering_JX1_HC_organoid_files/figure-gfm/pca-1.png)<!-- -->

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
    ## Number of nodes: 7802
    ## Number of edges: 272887
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9479
    ## Number of communities: 5
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7802
    ## Number of edges: 272887
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9259
    ## Number of communities: 7
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7802
    ## Number of edges: 272887
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9091
    ## Number of communities: 9
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7802
    ## Number of edges: 272887
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8952
    ## Number of communities: 12
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7802
    ## Number of edges: 272887
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8853
    ## Number of communities: 14
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7802
    ## Number of edges: 272887
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8768
    ## Number of communities: 16
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7802
    ## Number of edges: 272887
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8683
    ## Number of communities: 17
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7802
    ## Number of edges: 272887
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8598
    ## Number of communities: 18
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7802
    ## Number of edges: 272887
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8515
    ## Number of communities: 18
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 7802
    ## Number of edges: 272887
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8438
    ## Number of communities: 19
    ## Elapsed time: 0 seconds

``` r
# We redefine the final partitioning with resolution 0.5
SeuratObject$seurat_clusters <- SeuratObject$RNA_snn_res.0.5 
Idents(SeuratObject) <- SeuratObject$RNA_snn_res.0.5 
```

## Diagnostics of unsupervised clustering

We will investigate how the selection of multiple resolutions affects
the partition into individual cell
clusters.

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.")
```

![](1_initial_clustering_JX1_HC_organoid_files/figure-gfm/clustree-1.png)<!-- -->

We also checked whether mithochondrial expression drives the clustering

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.", 
     node_colour="percent.mt", node_colour_aggr="mean")
```

![](1_initial_clustering_JX1_HC_organoid_files/figure-gfm/clustree_mt-1.png)<!-- -->

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.", 
     node_colour="nFeature_RNA", node_colour_aggr="median")
```

![](1_initial_clustering_JX1_HC_organoid_files/figure-gfm/clustree_nFeat-1.png)<!-- -->

``` r
VlnPlot(SeuratObject, 
    feature=c("percent.mt", "nFeature_RNA", 
        "Dissociation", "G2M.Score"), 
    pt.size=0.4,
    ncol=2)
```

![](1_initial_clustering_JX1_HC_organoid_files/figure-gfm/vln_badQual_mtNfeat-1.png)<!-- -->

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

    ## 07:11:20 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 07:11:20 Read 7802 rows and found 25 numeric columns

    ## 07:11:20 Using Annoy for neighbor search, n_neighbors = 30

    ## 07:11:20 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 07:11:21 Writing NN index file to temp file /tmp/RtmpkYccrC/file7d3120756622
    ## 07:11:21 Searching Annoy index using 1 thread, search_k = 3000
    ## 07:11:23 Annoy recall = 100%
    ## 07:11:24 Commencing smooth kNN distance calibration using 1 thread
    ## 07:11:24 Initializing from normalized Laplacian + noise
    ## 07:11:24 Commencing optimization for 500 epochs, with 331636 positive edges
    ## 07:11:41 Optimization finished

``` r
DimPlot(SeuratObject, reduction="umap") + ggtitle(Project(SeuratObject))
```

![](1_initial_clustering_JX1_HC_organoid_files/figure-gfm/umap-1.png)<!-- -->

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
