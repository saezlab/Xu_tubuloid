ADPKD3 whole tissue biopsy, disease: initial clustering
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
SeuratObject <- getSeuratObject(path = "../data/sc/ADPKD3/filtered_feature_bc_matrix/",
                     project_name = "ADPKD3_kidney", 
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

Doing so, 5733 out of an initial 6422 total were retrieve for the
analysis, thus discarding a 689 cells.

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
    ## Positive:  LDB2, MEIS2, SH3RF3, DOCK4, PLPP1, TMTC1, EBF1, GNA14, CALD1, ERG 
    ##     PTPRB, CACNA1C, PLAT, SORBS1, TCF4, TIMP3, INSR, PREX2, SPARCL1, DLC1 
    ##     TACC1, NID1, FLT1, LHFPL6, TNS1, ELMO1, CRISPLD2, CDH13, ST6GALNAC3, FBXL7 
    ## Negative:  MACC1, SGMS2, HRH1, LAMA3, RP1, NCEH1, GRM8, AP002989.1, NRG3, GALNT17 
    ##     LAMC2, SH3RF1, LINC01320, MAL2, ERBB4, LMO7, MYO5B, PTH2R, AC012593.1, CAMK1D 
    ##     FAM135A, NEDD4L, CABLES1, KCTD8, SLC12A1, PAQR5, FGF13, TMEM178B, EDIL3, HDAC9 
    ## PC_ 2 
    ## Positive:  ACSM2A, SLC17A1, ACSM2B, KCNJ15, LRP2, CDH6, CUBN, SLC13A1, SLC2A9, SLC17A3 
    ##     SLC47A2, TINAG, SLC28A1, AC018467.1, CDH9, SLC47A1, GLYAT, PDZK1, ZNF804B, RHOBTB1 
    ##     AL589669.1, ADAMTS9-AS1, KCTD16, SLC16A9, SLC3A1, RASSF4, SLC16A12, CLDN10, DPYS, TRPM3 
    ## Negative:  LDB2, MEIS2, NFATC2, EBF1, CACNA1C, TIMP3, TMTC1, ERG, PLAT, EMP1 
    ##     TACC1, PLPP1, PREX2, PTPRB, ELMO1, LAMA2, SPARCL1, MT2A, CDH13, FKBP5 
    ##     TCF4, FILIP1, PDE3A, GNA14, FLT1, LHFPL6, HIF3A, SOX5, PDZRN4, DPYSL3 
    ## PC_ 3 
    ## Positive:  MECOM, AC019197.1, SCN2A, SCN3A, TFCP2L1, KCNIP4, ERBB4, CA12, CCSER1, LINC01606 
    ##     ATP6V0A4, PIK3C2G, AC138305.1, TOX, COBLL1, RHOBTB3, RALYL, AF233439.1, PAPPA2, PRDM16 
    ##     LHX1, DACH1, COL4A3, ACSL4, EFNA5, GATA3, DEPTOR, HSD11B2, SLC8A1, LINC01762 
    ## Negative:  NRG3, GALNT17, HRH1, SEMA5A, ARHGAP42, SGMS2, RP1, NEDD9, EDIL3, AP002989.1 
    ##     ITPKC, EMP1, RIPOR2, CAMK1D, MAMLD1, MYO5B, LAMA3, TACC1, SLC44A5, AKAP12 
    ##     TMEM178B, AC019117.1, CYP3A5, GULP1, CPNE4, NCEH1, EREG, TSHZ2, PRUNE2, MACC1 
    ## PC_ 4 
    ## Positive:  CRISPLD2, NEGR1, SLIT3, COL12A1, C7, LAMA2, GFPT2, PRKG1, DMD, CALD1 
    ##     DCLK1, SYNPO2, SLC19A2, TNC, PID1, ATRNL1, NID1, MMP19, CLMP, MICAL2 
    ##     NFASC, ADGRD1, SVEP1, PRR16, PALLD, CDH11, ABCA9, PDZRN4, LTBP1, CARMN 
    ## Negative:  PTPRB, FLT1, PLAT, EGFL7, VWF, PECAM1, ST6GALNAC3, ADGRL4, SLCO2A1, RAPGEF4 
    ##     ERG, PLEKHG1, SH2D3C, EPAS1, NOTCH4, EMCN, ENPP2, ANO2, TCF4, SEC14L1 
    ##     PRKCH, RALGAPA2, CD93, BTNL9, MYRIP, PLXNA2, ENG, DOCK9, PREX2, LRRC32 
    ## PC_ 5 
    ## Positive:  AC019197.1, MECOM, SCN2A, KCNIP4, ERBB4, PRKG1, KAZN, PLCB4, SAMD4A, DCDC2 
    ##     DACH1, MYO1E, CNKSR3, TMTC2, COL4A3, COBLL1, CCSER1, PLEKHG1, NHS, COL4A4 
    ##     PCSK5, KALRN, WWTR1, EFNA5, OSMR, SCN3A, EPAS1, SYTL2, DOCK9, SASH1 
    ## Negative:  IKZF1, PTPRC, CHST11, PRKCB, PIK3R5, ARHGAP15, AOAH, CD163, STAT4, PDE3B 
    ##     LAPTM5, SRGN, LCP2, SLC1A3, IQGAP2, SAMSN1, HCLS1, ITK, DOCK8, FYB1 
    ##     THEMIS, CELF2, DOCK10, TBXAS1, RHOH, IL7R, MRC1, INPP5D, IFNG-AS1, CD96

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
    ## Number of nodes: 5733
    ## Number of edges: 209799
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9678
    ## Number of communities: 11
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5733
    ## Number of edges: 209799
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9497
    ## Number of communities: 12
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5733
    ## Number of edges: 209799
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9321
    ## Number of communities: 14
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5733
    ## Number of edges: 209799
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9149
    ## Number of communities: 15
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5733
    ## Number of edges: 209799
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9014
    ## Number of communities: 18
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5733
    ## Number of edges: 209799
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8904
    ## Number of communities: 18
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5733
    ## Number of edges: 209799
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8796
    ## Number of communities: 18
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5733
    ## Number of edges: 209799
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8699
    ## Number of communities: 20
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5733
    ## Number of edges: 209799
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8607
    ## Number of communities: 21
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5733
    ## Number of edges: 209799
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8521
    ## Number of communities: 21
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

    ## 07:57:42 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 07:57:42 Read 5733 rows and found 25 numeric columns

    ## 07:57:42 Using Annoy for neighbor search, n_neighbors = 30

    ## 07:57:42 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 07:57:43 Writing NN index file to temp file /tmp/RtmponWKEX/file45d013ebe308
    ## 07:57:43 Searching Annoy index using 1 thread, search_k = 3000
    ## 07:57:44 Annoy recall = 100%
    ## 07:57:45 Commencing smooth kNN distance calibration using 1 thread
    ## 07:57:45 Initializing from normalized Laplacian + noise
    ## 07:57:45 Commencing optimization for 500 epochs, with 242004 positive edges
    ## 07:57:58 Optimization finished

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
