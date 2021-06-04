Number of cells during filtering in organoids
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
organoids <- c("CK5_organoid",
           "CK119_organoid",
           "JX1_HC_organoid",
           "JX2_PKD1KO_organoid",
           "JX3_PKD2KO_organoid")

cnts <- c("Initial", "afterQC", "Final")

res <- matrix(NA, nrow=length(cnts), ncol=length(organoids),
          dimnames=list(cnts, organoids))
```

``` r
rawDat_dirs <- c("CK5_organoid"="../data/sc/Organoid/filtered_feature_bc_matrix/",
         "CK119_organoid"="../data/sc/CK119_organoid/filtered_feature_bc_matrix/",
         "JX1_HC_organoid"="../data/sc/JX1_HC_organoid/filtered_feature_bc_matrix/",
         "JX2_PKD1KO_organoid"="../data/sc/JX2_PKD1KO_organoid/filtered_feature_bc_matrix/",
         "JX3_PKD2KO_organoid"="../data/sc/JX3_PKD2KO_organoid/filtered_feature_bc_matrix/")
stopifnot(all(names(rawDat_dirs)==colnames(res)))

res[1:2,] <- sapply(names(rawDat_dirs), function(sid) {
           cat(paste0("Reading:",sid,"\n"), file=stdout())
        S <- getSeuratObject(path=rawDat_dirs[sid],
                project_name = sid,
                mt.pattern = "^MT-", min.cells=5,
                min.features=200)
        nCells <- ncol(S)
        QC <- ncol(subset(S, nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20))
        c(nCells,QC)
         })
```

    ## Reading:CK5_organoid
    ## Reading:CK119_organoid

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Reading:JX1_HC_organoid
    ## Reading:JX2_PKD1KO_organoid
    ## Reading:JX3_PKD2KO_organoid

``` r
final_list <- c("CK5_organoid"="../Individual_analysis_CK5_early_organoid/output/4_final_assignment/active_idents.csv",
          "CK119_organoid"="../Individual_analysis_CK119_late_organoid/output/4_final_assignment/active_idents.csv",
          "JX1_HC_organoid"="../Individual_analysis_JX1_HC_organoid/output/4_final_assignment/active_idents.csv",
          "JX2_PKD1KO_organoid"="../Individual_analysis_JX2_PKD1KO_organoid/output/4_final_assignment/active_idents.csv",
          "JX3_PKD2KO_organoid"="../Individual_analysis_JX3_PKD2KO_organoid/output/4_final_assignment/active_idents.csv")
stopifnot(all(names(final_list)==colnames(res)))
res[3,] <- sapply(names(final_list), function(sid) {
        cat(paste0("Reading:", sid, "\n"), file=stdout())
        tab <- read.table(final_list[sid], sep=",", header=TRUE, stringsAsFactors=FALSE)
        return(nrow(tab))
          })
```

    ## Reading:CK5_organoid
    ## Reading:CK119_organoid
    ## Reading:JX1_HC_organoid
    ## Reading:JX2_PKD1KO_organoid
    ## Reading:JX3_PKD2KO_organoid

## Save table

``` r
write.table(res, file=paste0(OUTDIR,"/NumberOfCells_organoids.csv"), 
        sep=",", row.names=TRUE, col.names=NA, quote=FALSE)
```
