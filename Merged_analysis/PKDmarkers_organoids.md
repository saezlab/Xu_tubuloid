PKD marker expression from PKD organoids
================
Javier Perales-Paton - <javier.perales@bioquant.uni-heidelberg.de>

## Load libraries and auxiliar functions

``` r
set.seed(1234)
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(cowplot))
suppressPackageStartupMessages(require(ComplexHeatmap))
suppressPackageStartupMessages(require(GSEABase))
source("../src/seurat_fx.R")
```

## Load SeuratObject with initial clustering outcome

``` r
CK120_CD13 <- readRDS("../Individual_analysis_CK120_CD13/output/4_final_assignment/data/SeuratObject.rds")

CK121_CD24 <- readRDS("../Individual_analysis_CK121_CD24/output/4_final_assignment/data/SeuratObject.rds")

early <- readRDS(paste0("../Individual_analysis_CK5_early_organoid/",
            "output/4_final_assignment/data/SeuratObject.rds")
)
early$final_Ident <- Idents(early)

late <- readRDS(paste0("../Individual_analysis_CK119_late_organoid/",
            "output/4_final_assignment/data/SeuratObject.rds")
)
late$final_Ident <- Idents(late)

HC <- readRDS(paste0("../Individual_analysis_JX1_HC_organoid/",
             "output/4_final_assignment/data/SeuratObject.rds")
)
HC$final_Ident <- Idents(HC)

ADPKD_PKD1KO <- readRDS(paste0("../Individual_analysis_JX2_PKD1KO_organoid/",
                "output/4_final_assignment/data/SeuratObject.rds"))
ADPKD_PKD1KO$final_Ident <- Idents(ADPKD_PKD1KO)

ADPKD_PKD2KO <- readRDS(paste0("..//Individual_analysis_JX3_PKD2KO_organoid/",
                "output/4_final_assignment/data/SeuratObject.rds"))
ADPKD_PKD2KO$final_Ident <- Idents(ADPKD_PKD2KO)
```

## Merge

``` r
SL <- list(CK120_CD13, CK121_CD24, early, late, HC, ADPKD_PKD1KO, ADPKD_PKD2KO)
rm(CK120_CD13, CK121_CD24, early, late, HC, ADPKD_PKD1KO, ADPKD_PKD2KO)
S <- merge(SL[[1]], SL[-1])
```

    ## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names
    ## are duplicated across objects provided. Renaming to enforce unique cell
    ## names.

## Heatmap

``` r
cols <- readRDS(file="./output/color_scheme.rds")
```

``` r
# Curated list of PKD markers from literature
genes <- c("STAT3", "TGFB1", "MET", "LGALS3", "NDRG1")
```

``` r
Sx<- S[intersect(rownames(S), genes), ]
Sx<- ScaleData(Sx, verbose = FALSE)

GSC <- GeneSetCollection(GeneSet(genes, setName=" "))

CK5_hp <- DoHeatmap2(SeuratObject = subset(Sx, orig.ident=="CK5_organoid"),
        row_names_fontisze=12,
           res=NULL, cols=cols,
           assay="RNA", width = unit(41,"mm"),name="",
           show_hr = FALSE,
           GSC=GSC)

CK119_hp <- DoHeatmap2(SeuratObject = subset(Sx, orig.ident=="CK119_organoid"),
        row_names_fontisze=12,
           res=NULL, cols=cols,
           assay="RNA", width = unit(41,"mm"),name="",
           show_hr = FALSE,
           GSC=GSC)
JX2_hp <- DoHeatmap2(SeuratObject = subset(Sx, orig.ident=="JX2_PKD1KO_organoid"),
        row_names_fontisze=12,
           res=NULL, cols=cols,
           assay="RNA", width = unit(41,"mm"),name="",
           show_hr = FALSE,
           GSC=GSC)
JX3_hp <- DoHeatmap2(SeuratObject = subset(Sx, orig.ident=="JX3_PKD2KO_organoid"),
        row_names_fontisze=12,
           res=NULL, cols=cols,
           assay="RNA", width = unit(41,"mm"),name="",
           show_hr = FALSE,
           GSC=GSC)



ht_list <-  JX2_hp + JX3_hp + CK5_hp + CK119_hp 
```

    ## Warning: Heatmap/annotation names are duplicated:

    ## Warning: Heatmap/annotation names are duplicated: ,

    ## Warning: Heatmap/annotation names are duplicated: , ,

``` r
# draw(ht_list, ht_gap = unit(10, "mm"))
draw(ht_list, ht_gap = unit(10, "mm"), 
    column_title = "PKD1-KO and PKD2-KO (PKD), Early and late control organoids", 
    column_title_gp = gpar(fontsize = 16))
```

![](./PKDmarkers_organoids//figures/heatmap_PKDmarkers_organoids-1.png)<!-- -->

## Pseudobulking

We collapse single-cell information to express it as a pseudobulk
profile for simplicity.

``` r
mat <- scater::sumCountsAcrossCells(as.matrix(S@assays$RNA@counts),
                         ids = S$orig.ident
                         )
CPM <- edgeR::cpm(mat)[genes, ]

# NOTE: this is same outcome as 
#   genesorteR::sortGenes(S@assays$RNA@counts, S$orig.ident)$condGeneProb[genes,]
PCT <- lapply(SplitObject(S, split.by="orig.ident"), function(Sj) { 
        rowSums(Sj@assays$RNA@counts[genes, ] > 0)/ncol(Sj)
                         })
PCT <- t(do.call("rbind", PCT)) * 100

# Same order
ord <- c("JX2_PKD1KO_organoid", "JX3_PKD2KO_organoid",
     "CK5_organoid", "CK119_organoid") # drop JX1, PKD-ctrl order
CPM <- CPM[, ord]
PCT <- PCT[, ord]
```

``` r
hp1 <-Heatmap(t(scale(t(CPM))),
          name="Gene Expr.\n(CPM, row scaled)",
    column_title = "Pseudobulk",
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_side = "left", column_names_side = "top",
    row_names_gp = gpar(fontsize=12))

col_fun <- circlize::colorRamp2(c(0, 50, 100), c("blue", "#EEEEEE", "red"))
hp2 <- Heatmap(PCT, col=col_fun,
    column_title = "Positive expressing cells",
     cluster_rows = FALSE, cluster_columns = FALSE,
     row_names_side = "left", column_names_side="top", 
     row_names_gp = gpar(fontsize = 10),
     cell_fun = function(j, i, x, y, width, height, fill) {
         grid.text(paste0(sprintf("%.2f", PCT[i, j])),
               x, y, 
              gp = gpar(fontsize = 12))
             } ,
             name="Perc.(%)")
draw(hp1 + hp2, ht_gap = unit(10, "mm"))
```

![](./PKDmarkers_organoids//figures/heatmap_PKDmarkers_Organoids_bulk-1.png)<!-- -->
