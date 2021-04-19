CK5 early healthy organoid: initial clustering
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
(\<20%). Later we will apply diagnostics to see if this is related.

``` r
nSamples_before <- ncol(SeuratObject)
SeuratObject <- subset(SeuratObject, nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
nSamples_after <- ncol(SeuratObject)
```

Doing so, 2403 out of an initial 2792 total were retrieve for the
analysis, thus discarding a 389 cells.

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
    ## Positive:  TUBA1B, H2AFZ, RAN, TUBB, CFL1, HMGA1, JPT1, KRT18, ACTB, CENPW 
    ##     S100A16, BIRC5, TK1, RANBP1, ARPC2, LSM4, CKS1B, DTYMK, TUBB4B, SRSF3 
    ##     CCND1, KRT8, CDKN3, OAZ1, TUBA1C, TMSB4X, SNRPB, PHLDA2, TXN, PFN1 
    ## Negative:  CLU, NEAT1, SPP1, LTF, C3, SLC34A2, SLPI, MUC1, BCAM, KCNQ1OT1 
    ##     CP, SERPING1, SOD2, ADAMTS1, WFDC2, IFITM3, SAA1, SLC3A1, SAA2, C1R 
    ##     NCOA7, PIGR, CDH16, SERPINA1, GPNMB, ANXA4, LCN2, TMEM37, C1S, FTL 
    ## PC_ 2 
    ## Positive:  FN1, IL32, MAP1B, TPM1, RAB3B, SPARC, CDH6, TMSB10, DCBLD2, FLG 
    ##     CDH2, MFGE8, TUBA1A, KRT23, AKAP12, THY1, SLIT3, CALD1, S100A6, SPOCK2 
    ##     ITIH5, RHOC, C6orf99, EMP3, GLIPR2, THBS2, SH3BGRL3, CD70, CYTL1, PTGIS 
    ## Negative:  WFDC2, HMGB2, MKI67, CENPF, NUSAP1, TOP2A, CDCA3, CDC20, PBK, UBE2C 
    ##     AGR2, MNS1, DLGAP5, AURKB, PLK1, HMMR, BCAM, DEFB1, TPX2, CCNB1 
    ##     DEPDC1, CCNB2, GTSE1, BIRC5, ASPM, CDKN3, CCNA2, SPC25, ANP32E, TK1 
    ## PC_ 3 
    ## Positive:  SPARC, FLG, THY1, CENPF, CDH2, ITIH5, MKI67, CD70, SLIT3, VCAM1 
    ##     ASPM, NUSAP1, PTGIS, DEPDC1, HMGB2, TOP2A, GTSE1, CDC20, CYTL1, CRB2 
    ##     PTTG1, SPC25, PLK1, IGFBP4, CDKN3, UBE2C, HMMR, SEPT4, CCNA2, PBK 
    ## Negative:  TSPAN1, CLDN7, RAB25, S100A14, KRT7, TACSTD2, GPRC5A, TMPRSS4, EPCAM, CLDN4 
    ##     CST6, ITGA2, PRSS8, KRT19, CDH1, MAL2, UCA1, UCP2, SFTA2, MACC1 
    ##     PROM2, SAT1, SCEL, ITGB6, CDA, AGR2, VTCN1, ST14, SFN, GATA3 
    ## PC_ 4 
    ## Positive:  FLG, THY1, SLIT3, ITIH5, KRT23, PTGIS, SPARC, CD14, CYTL1, CRB2 
    ##     IGFBP4, C6orf99, PRRX2, FGF18, SEPT4, AP1M2, UPK3B, KRT19, TMOD1, BCAT1 
    ##     SNCA, NGFR, BCAM, CFTR, MAF, RGS5, CD74, AC005482.1, UPK1B, DEFB1 
    ## Negative:  SERPINE2, MT2A, PLAU, GYPC, EMP3, VIM, CXCL14, PHLDA1, MT1E, TRIM55 
    ##     XKR4, GMDS, PKIB, ANPEP, CLDN16, GLRX, TPM1, TM4SF1, LGALS1, DCBLD2 
    ##     HTRA1, ADIRF, S100A6, SPON1, TNC, GCHFR, ITGA3, MSMP, TGFBI, IGFBP6 
    ## PC_ 5 
    ## Positive:  LDHB, DBI, NNMT, ATP5F1B, PRDX1, MYL12B, ENO1, ANXA4, CPVL, SLC25A5 
    ##     PSMA7, SPP1, HSPE1, RPS2, ACAT2, RPL21, LDHA, FDPS, TNFSF10, GSTM3 
    ##     CD9, NQO1, UCHL1, HLA-DMB, FHL2, GSTP1, NDUFA9, ANXA2, FDFT1, CCT2 
    ## Negative:  NEAT1, SOX4, TOP2A, MKI67, UBE2C, ASPM, NDRG1, HSPG2, AREG, PLK1 
    ##     MFGE8, GTSE1, KIF23, TPX2, UCA1, KIF2C, KLF6, KNL1, CDCA8, AURKA 
    ##     HJURP, ANLN, AURKB, CENPF, CCNA2, CENPA, RAB3B, TMEM132A, MMP14, DEPDC1

``` r
print(ElbowPlot(SeuratObject,ndims = 50) + geom_vline(xintercept = 25, col="red"))
```

![](1_initial_clustering_CK5_organoid_files/figure-gfm/pca-1.png)<!-- -->

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
    ## Number of nodes: 2403
    ## Number of edges: 82017
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9399
    ## Number of communities: 4
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2403
    ## Number of edges: 82017
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9057
    ## Number of communities: 5
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2403
    ## Number of edges: 82017
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8759
    ## Number of communities: 7
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2403
    ## Number of edges: 82017
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8554
    ## Number of communities: 9
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2403
    ## Number of edges: 82017
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8402
    ## Number of communities: 10
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2403
    ## Number of edges: 82017
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8272
    ## Number of communities: 11
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2403
    ## Number of edges: 82017
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8152
    ## Number of communities: 12
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2403
    ## Number of edges: 82017
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8032
    ## Number of communities: 12
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2403
    ## Number of edges: 82017
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7923
    ## Number of communities: 13
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2403
    ## Number of edges: 82017
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7824
    ## Number of communities: 15
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

![](1_initial_clustering_CK5_organoid_files/figure-gfm/clustree-1.png)<!-- -->

We also checked whether mithochondrial expression drives the clustering

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.", 
     node_colour="percent.mt", node_colour_aggr="mean")
```

![](1_initial_clustering_CK5_organoid_files/figure-gfm/clustree_mt-1.png)<!-- -->

``` r
clustree(SeuratObject, prefix = "RNA_snn_res.", 
     node_colour="nFeature_RNA", node_colour_aggr="median")
```

![](1_initial_clustering_CK5_organoid_files/figure-gfm/clustree_nFeat-1.png)<!-- -->

``` r
VlnPlot(SeuratObject, 
    feature=c("percent.mt", "nFeature_RNA", 
        "Dissociation", "G2M.Score"), 
    pt.size=0.4,
    ncol=2)
```

![](1_initial_clustering_CK5_organoid_files/figure-gfm/vln_badQual_mtNfeat-1.png)<!-- -->

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

    ## 14:56:30 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 14:56:30 Read 2403 rows and found 25 numeric columns

    ## 14:56:30 Using Annoy for neighbor search, n_neighbors = 30

    ## 14:56:30 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 14:56:31 Writing NN index file to temp file /tmp/RtmpXkpO7H/file1c9d17a32fd3
    ## 14:56:31 Searching Annoy index using 1 thread, search_k = 3000
    ## 14:56:31 Annoy recall = 100%
    ## 14:56:32 Commencing smooth kNN distance calibration using 1 thread
    ## 14:56:32 Initializing from normalized Laplacian + noise
    ## 14:56:32 Commencing optimization for 500 epochs, with 95292 positive edges
    ## 14:56:38 Optimization finished

``` r
DimPlot(SeuratObject, reduction="umap") + ggtitle(Project(SeuratObject))
```

![](1_initial_clustering_CK5_organoid_files/figure-gfm/umap-1.png)<!-- -->

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
