CK121 sorted CD24+ cells : final cell assignment
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
SeuratObject <- readRDS("./output/3_refine_clustering//data/SeuratObject.rds")
```

## Define output directory

``` r
# Define output directory
OUTDIR <- paste0("./output/4_final_assignment/")
if(! dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
```

## Identify the proliferating PT-like cell population

``` r
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

SeuratObject <- CellCycleScoring(SeuratObject, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
d3 <- DimPlot(SeuratObject, group.by = "Phase") + ggtitle("Cell Cycle Phase")
d4 <- DimPlot(SeuratObject) + ggtitle("Re-classification")

print(CombinePlots(list(d3,d4)))
```

![](4_final_assignment_CK121_CD24_files/figure-gfm/cell_phase-1.png)<!-- -->

## Cell marker extraction

We are going to use two methods for the extraction of cell markers.  
1\. Differential gene expression using wilcox test from Seurat. Only
positive and consistent markers will be tested (with default
parameters).  
2\. GenesorteR to get posterior probabilities for each gene of observing
a certain cell population.

### 1 wilcox test

``` r
up <- setNames(vector("list",length=length(levels(SeuratObject))), 
               levels(SeuratObject))
for(idx in names(up)) {
  up.idx <- FindMarkers(SeuratObject,ident.1 = idx, 
                        ident.2 = setdiff(levels(SeuratObject), idx), only.pos=T)
  cols_names <- colnames(up.idx)
  
  # Add two extra cols
  up.idx$cluster <- idx
  up.idx$gene <- rownames(up.idx)
  
  up[[idx]] <- up.idx
}
```

### 2 Gene sorter

``` r
sg <- sortGenes(SeuratObject@assays$RNA@data, Idents(SeuratObject))
```

    ## Warning in sortGenes(SeuratObject@assays$RNA@data, Idents(SeuratObject)):
    ## A Friendly Warning: Some genes were removed because they were zeros in all
    ## cells after binarization. You probably don't need to do anything but you
    ## might want to look into this. Maybe you forgot to pre-filter the genes? You
    ## can also use a different binarization method. Excluded genes are available
    ## in the output under '$removed'.

``` r
#define a small set of markers
#mm = getMarkers(sg, quant = 0.975)

#cluster genes and make a heatmap
#pp = plotMarkerHeat(sg$inputMat, sg$inputClass, mm$markers, clusterGenes=TRUE, outs = TRUE)

#pp$gene_class_info #gene clusters

#the top 25 genes for each cluster by specificity scores
top_markers = apply(sg$specScore, 2, function(x) names(head(sort(x, decreasing = TRUE), n = 25)))
```

Finally we save the markers for manual exploration:

``` r
MARKERS_OUTDIR <- paste0(OUTDIR,"/Markers")
if(! dir.exists(MARKERS_OUTDIR)) dir.create(MARKERS_OUTDIR, recursive = TRUE)

# Wilcox
## TSV
write.table(do.call("rbind",up)[c("cluster", "gene", cols_names)],
        file = paste0(MARKERS_OUTDIR,"/",Project(SeuratObject),"_wilcox_DEGs_up.csv"),
        sep=",", col.names=TRUE, row.names = FALSE, quote=FALSE
        )
## bin
saveRDS(up, file=paste0(MARKERS_OUTDIR,"/wilcox_up.rds"))

# Genesorter
saveRDS(sg, file=paste0(MARKERS_OUTDIR,"/genesorter_out.rds"))
write.table(as.matrix(sg$specScore), 
        paste0(MARKERS_OUTDIR,"/",Project(SeuratObject),"_specScore.csv"),
            sep=",",row.names = TRUE,col.names = NA,quote=FALSE)

write.table(as.matrix(sg$condGeneProb), 
        paste0(MARKERS_OUTDIR,"/",Project(SeuratObject),"_condGeneProb.csv"),
            sep=",",row.names = TRUE,col.names = NA,quote=FALSE)
```

## Final cell assignment

The final table of markers on the second round of cell clustering give
the following results:

| Cluster   | rational                             |
| :-------- | :----------------------------------- |
| Cluster 0 | TPC : 43% AKAP12                     |
| Cluster 1 | TPC : 77% AKAP12, 67% VCAM           |
| Cluster 2 | LOH, TAL cells : PAX8, SLC12A1, UMOD |
| Cluster 3 | LOH, TAL cells : PAX8, 29% SLC12A1   |
| Cluster 4 | DCT                                  |
| Cluster 5 | IC-A                                 |
| Cluster 6 | PEC cells                            |

``` r
ren_id <- c("0"="TPC_2",
            "1"="TPC_1",
            "2"="LOH, TAL cells_1",
            "3"="LOH, TAL cells_2",
            "4"="DCT",
            "5"="IC-A",
            "6"="PEC")
SeuratObject <- RenameIdents(SeuratObject, ren_id)
```

## TSNE

``` r
DimPlot(SeuratObject, reduction="tsne")
```

![](4_final_assignment_CK121_CD24_files/figure-gfm/tsne_final-1.png)<!-- -->

## UMAP

``` r
DimPlot(SeuratObject, reduction="umap")
```

![](4_final_assignment_CK121_CD24_files/figure-gfm/umap_final-1.png)<!-- -->

## Archive processed data for downstream analysis

``` r
write.table(data.frame("assign"=ren_id),
        file=paste0(OUTDIR,"/assignment.csv"),
        sep=",", col.names = NA, row.names=TRUE, quote=FALSE)
# final idents
write.table(data.frame("Ident"=SeuratObject@active.ident,
               "seurat_clusters"=SeuratObject$seurat_clusters),
            file=paste0(OUTDIR,"/active_idents.csv"),
            sep=",", col.names = NA, row.names = TRUE, quote=TRUE)
```

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
    ##  [1] cowplot_1.0.0        clustree_0.4.1       ggraph_2.0.0.9000   
    ##  [4] ComplexHeatmap_2.0.0 genesorteR_0.3.1     Matrix_1.2-17       
    ##  [7] dplyr_0.8.3          GSEABase_1.46.0      graph_1.62.0        
    ## [10] annotate_1.62.0      XML_3.98-1.20        AnnotationDbi_1.46.1
    ## [13] IRanges_2.18.2       S4Vectors_0.22.1     Biobase_2.44.0      
    ## [16] BiocGenerics_0.30.0  ggplot2_3.2.1        Seurat_3.1.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.15          colorspace_1.4-1    rjson_0.2.20       
    ##   [4] ggridges_0.5.1      mclust_5.4.5        circlize_0.4.7     
    ##   [7] GlobalOptions_0.1.0 clue_0.3-57         farver_1.1.0       
    ##  [10] leiden_0.3.1        listenv_0.7.0       npsurv_0.4-0       
    ##  [13] graphlayouts_0.5.0  ggrepel_0.8.1       bit64_0.9-7        
    ##  [16] codetools_0.2-16    splines_3.6.1       R.methodsS3_1.7.1  
    ##  [19] lsei_1.2-0          knitr_1.24          polyclip_1.10-0    
    ##  [22] zeallot_0.1.0       jsonlite_1.6        ica_1.0-2          
    ##  [25] cluster_2.1.0       png_0.1-7           R.oo_1.22.0        
    ##  [28] pheatmap_1.0.12     uwot_0.1.4          ggforce_0.3.1      
    ##  [31] sctransform_0.2.0   compiler_3.6.1      httr_1.4.1         
    ##  [34] backports_1.1.4     assertthat_0.2.1    lazyeval_0.2.2     
    ##  [37] tweenr_1.0.1        htmltools_0.3.6     tools_3.6.1        
    ##  [40] rsvd_1.0.2          igraph_1.2.4.1      gtable_0.3.0       
    ##  [43] glue_1.3.1          RANN_2.6.1          reshape2_1.4.3     
    ##  [46] Rcpp_1.0.2          vctrs_0.2.0         gdata_2.18.0       
    ##  [49] ape_5.3             nlme_3.1-141        gbRd_0.4-11        
    ##  [52] lmtest_0.9-37       xfun_0.9            stringr_1.4.0      
    ##  [55] globals_0.12.4      lifecycle_0.1.0     irlba_2.3.3        
    ##  [58] gtools_3.8.1        future_1.14.0       MASS_7.3-51.4      
    ##  [61] zoo_1.8-6           scales_1.0.0        tidygraph_1.1.2    
    ##  [64] RColorBrewer_1.1-2  yaml_2.2.0          memoise_1.1.0      
    ##  [67] reticulate_1.13     pbapply_1.4-2       gridExtra_2.3      
    ##  [70] stringi_1.4.3       RSQLite_2.1.2       caTools_1.17.1.2   
    ##  [73] bibtex_0.4.2        shape_1.4.4         Rdpack_0.11-0      
    ##  [76] SDMTools_1.1-221.1  rlang_0.4.0         pkgconfig_2.0.3    
    ##  [79] bitops_1.0-6        evaluate_0.14       lattice_0.20-38    
    ##  [82] ROCR_1.0-7          purrr_0.3.2         labeling_0.3       
    ##  [85] htmlwidgets_1.3     bit_1.1-14          tidyselect_0.2.5   
    ##  [88] RcppAnnoy_0.0.13    plyr_1.8.4          magrittr_1.5       
    ##  [91] R6_2.4.0            gplots_3.0.1.1      DBI_1.0.0          
    ##  [94] pillar_1.4.2        withr_2.1.2         fitdistrplus_1.0-14
    ##  [97] survival_2.44-1.1   RCurl_1.95-4.12     tibble_2.1.3       
    ## [100] future.apply_1.3.0  tsne_0.1-3          crayon_1.3.4       
    ## [103] KernSmooth_2.23-16  plotly_4.9.0        rmarkdown_1.15     
    ## [106] viridis_0.5.1       GetoptLong_0.1.7    data.table_1.12.8  
    ## [109] blob_1.2.0          metap_1.1           digest_0.6.21      
    ## [112] xtable_1.8-4        tidyr_1.0.0         R.utils_2.9.0      
    ## [115] RcppParallel_4.4.3  munsell_0.5.0       viridisLite_0.3.0

``` r
{                                                                                                                                                                                                           
sink(file=paste0(OUTDIR,"/sessionInfo.txt"))
print(sessionInfo())
sink()
}
```
