Number of cells during filtering in sorted
================
Javier Perales-Paton

``` r
set.seed(1234)
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(reshape2))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(ggrepel))
suppressPackageStartupMessages(require(cowplot))
suppressPackageStartupMessages(require(ComplexHeatmap))
suppressPackageStartupMessages(require(cowplot))
source("../src/seurat_fx.R")
```

## Prepare environment

Define the output directory for figures and tables

``` r
OUTDIR <- "./NumberOfCells/"
if(!dir.exists(OUTDIR)) dir.create(OUTDIR)
```

## Read raw data is in original analysis

``` r
sorted <- c("CK120_CD13", "CK121_CD24")

cnts <- c("Initial", "afterQC", "Final")

res <- matrix(NA, nrow=length(cnts), ncol=length(sorted),
          dimnames=list(cnts, sorted))
```

``` r
rawDat_dirs <- c("CK120_CD13"="../data/sc/CK120_CD13/filtered_feature_bc_matrix/",
         "CK121_CD24"="../data/sc/CK121_CD24/filtered_feature_bc_matrix/")
stopifnot(all(names(rawDat_dirs)==colnames(res)))

res[1:2,] <- sapply(names(rawDat_dirs), function(sid) {
           cat(paste0("Reading:",sid,"\n"), file=stdout())
        S <- getSeuratObject(path=rawDat_dirs[sid],
                project_name = sid,
                mt.pattern = "^MT-", min.cells=5,
                min.features=200)
        nCells <- ncol(S)
        QC <- ncol(subset(S, nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 80))
        c(nCells,QC)
         })
```

    ## Reading:CK120_CD13

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Reading:CK121_CD24

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

``` r
final_list <- c("CK120_CD13"="../Individual_analysis_CK120_CD13/output/4_final_assignment/active_idents.csv",
          "CK121_CD24"="../Individual_analysis_CK121_CD24/output/4_final_assignment/active_idents.csv")
stopifnot(all(names(final_list)==colnames(res)))
res[3,] <- sapply(names(final_list), function(sid) {
        cat(paste0("Reading:", sid, "\n"), file=stdout())
        tab <- read.table(final_list[sid], sep=",", header=TRUE, stringsAsFactors=FALSE)
        return(nrow(tab))
          })
```

    ## Reading:CK120_CD13
    ## Reading:CK121_CD24

## Save table

``` r
write.table(res, file=paste0(OUTDIR,"/NumberOfCells_sorted.csv"), 
        sep=",", row.names=TRUE, col.names=NA, quote=FALSE)
```
