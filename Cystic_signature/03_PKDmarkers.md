ADPKD marker expression in cystic populations
================
Javier Perales-Paton - <javier.perales@bioquant.uni-heidelberg.de>

## Load libraries and auxiliar functions

``` r
set.seed(1234)
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(clustree))
suppressPackageStartupMessages(require(cowplot))
suppressPackageStartupMessages(require(GSEABase))
source("../src/seurat_fx.R")
```

## Load SeuratObject with initial clustering outcome

``` r
CK224 <- readRDS(paste0("../Individual_analysis_CK224_PT25_PKD2-/",
            "output/2_cell_assignment/data/SeuratObject.rds"))
CK224$final_Ident <- Idents(CK224)

CK225 <- readRDS(paste0("../Individual_analysis_CK225_PT8_PKD1-/",
            "output/2_cell_assignment/data/SeuratObject.rds"))
CK225$final_Ident <- Idents(CK225)

Control1 <- readRDS(paste0("../Individual_analysis_Control1/",
            "output/2_cell_assignment/data/SeuratObject.rds"))
Control1$final_Ident <- Idents(Control1)

Control2 <- readRDS(paste0("../Individual_analysis_Control2/",
            "output/2_cell_assignment/data/SeuratObject.rds"))
Control2$final_Ident <- Idents(Control2)

ADPKD3 <- readRDS(paste0("../Individual_analysis_ADPKD3/",
            "output/2_cell_assignment/data/SeuratObject.rds"))
ADPKD3$final_Ident <- Idents(ADPKD3)

ADPKD3 <- ADPKD3[ , Idents(ADPKD3)!=":LowQual"]
```

## Merging

``` r
S <- merge(ADPKD3, list(CK224, CK225, Control1, Control2))
```

    ## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names
    ## are duplicated across objects provided. Renaming to enforce unique cell
    ## names.

``` r
rm(CK224, CK225, Control1, Control2, ADPKD3)
```

## Selection of PKD markers

``` r
# Curated list of PKD markers from literature
genes <- c("STAT3", "TGFB1", "MET", "LGALS3", "NDRG1")

plots <- DotPlot_panel(S, assay = "RNA",
                       intersect(genes, rownames(S)), 
               dot.scale = 12,
              scale.min = 0, scale.max = 100, col.min = -2.5, col.max = 2.5)

# Common scale
plots <- lapply(plots, function(gg) {
  gg + coord_flip() + scale_y_discrete(position = "right") +
    
  # This is extremely important to use same scaling color for all samples
  scale_color_gradient(low="lightgrey", high = "blue",limits=c(-1.5, 2.5)) +
      
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        plot.title = element_text(size=18, hjust = 0.5))
  })
```

    ## Scale for 'colour' is already present. Adding another scale for
    ## 'colour', which will replace the existing scale.
    ## Scale for 'colour' is already present. Adding another scale for
    ## 'colour', which will replace the existing scale.
    ## Scale for 'colour' is already present. Adding another scale for
    ## 'colour', which will replace the existing scale.
    ## Scale for 'colour' is already present. Adding another scale for
    ## 'colour', which will replace the existing scale.
    ## Scale for 'colour' is already present. Adding another scale for
    ## 'colour', which will replace the existing scale.

``` r
# Remove xlab
plots <- lapply(plots, function(gg) gg + theme(axis.title.x = element_blank()))

# Remove legend from the first three ones
plots[1:4] <- lapply(plots[1:4], function(gg) gg + NoLegend())
# Remove y axis title from the last three ones
plots[-1] <- lapply(plots[-1], function(gg) gg + theme(axis.title.y = element_blank()))

CombinePlots(plots,
            rel_widths=c(7, 6, 7, 6, 9),
             ncol = 5)
```

    ## Warning: Graphs cannot be vertically aligned unless the axis parameter is
    ## set. Placing graphs unaligned.

![](./03_PKDmarkers//figures/dotplot_PKDmarkers_tissue-1.png)<!-- -->

``` r
S<- S[intersect(rownames(S), genes), ]
S<- ScaleData(S, verbose = FALSE)

GSC <- GeneSetCollection(GeneSet(genes, setName=" "))

Control1_hp <- DoHeatmap2(SeuratObject = subset(S, orig.ident=="Control1_kidney"),
        row_names_fontisze=12,
           res=NULL,# cols=cols,
           assay="RNA", name="",
           show_hr = FALSE,
           GSC=GSC)

Control2_hp <- DoHeatmap2(SeuratObject = subset(S, orig.ident=="Control1_kidney"),
        row_names_fontisze=12,
           res=NULL,# cols=cols,
           assay="RNA", name="",
           show_hr = FALSE,
           GSC=GSC)

CK224_hp <- DoHeatmap2(SeuratObject = subset(S, orig.ident=="CK224_kidney_PKD2-"),
        row_names_fontisze=12,
           res=NULL,# cols=cols,
           assay="RNA", name="",
           show_hr = FALSE,
           GSC=GSC)

CK225_hp <- DoHeatmap2(SeuratObject = subset(S, orig.ident=="CK225_kidney_PKD1-"),
        row_names_fontisze=12,
           res=NULL,# cols=cols,
           assay="RNA", name="",
           show_hr = FALSE,
           GSC=GSC)

ADPKD3_hp <- DoHeatmap2(SeuratObject = subset(S, orig.ident=="ADPKD3_kidney"),
        row_names_fontisze=12,
           res=NULL,# cols=cols,
           assay="RNA", name="",
           show_hr = FALSE,
           GSC=GSC)



ht_list <- ADPKD3_hp + CK224_hp + CK225_hp + Control1_hp + Control2_hp
```

    ## Warning: Heatmap/annotation names are duplicated:

    ## Warning: Heatmap/annotation names are duplicated: ,

    ## Warning: Heatmap/annotation names are duplicated: , ,

    ## Warning: Heatmap/annotation names are duplicated: , , ,

``` r
# draw(ht_list, ht_gap = unit(10, "mm"))
draw(ht_list, ht_gap = unit(10, "mm"), 
    column_title = "ADPKD3, CK224, CK225 (ADPKD) and Control1, Control2 (healthy) tissue", 
    column_title_gp = gpar(fontsize = 16))
```

![](./03_PKDmarkers//figures/heatmap_PKDmarkers_tissue-1.png)<!-- -->
