#' Visualize Symphony mapping results
#' Javier Perales-Patón (c)

## Environment reproducibility
set.seed(1234)
OUTDIR <- "./02_symphony_vis/"
if(!dir.exists(OUTDIR)) dir.create(OUTDIR);

## Load libraries
library(Seurat)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(ComplexHeatmap)
library(tibble)
source("utils_seurat.R")

## Visualization settings
# Color scheme for REFERENCE atlas (tissue)
gg_color_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}
# Color scheme for merged Organoids
col_org <- readRDS("../Merged_analysis/output/color_scheme.rds")

# Color scheme for human clusters
# - List of human clusters annotated
ren_id <- c("0"="PC-CD/CNT",
	    "1"="PT_1",
	    "2"="TAL_1",
	    "3"="TAL_2",
	    "4"="DCT_1",
	    "5"="EC",
	    "6"="PT_3/PEC",
	    "7"="PT_4",
	    "8"="PT_2",
	    "9"="vSMC/Fib",
	    "10"="T-Cell",
	    "11"="IC-A",
	    "12"="Pod", # Podocytes
	    "13"="B-Cell",
	    "14"="Mac",
	    "15"="IC-B",
	    "16"="DCT_2",
	    "17"="Mast-Cell",
	    "18"="Unk_1",
	    "19"="Unk_2")

col_map <- setNames(gg_color_hue(length(ren_id)),
		    sort(ren_id))

## Load data
## Organoid data
list_rds <- c("CK5_organoid"="../Individual_analysis_CK5_early_organoid/output/4_final_assignment/data/SeuratObject.rds",
	      "CK119_organoid"="../Individual_analysis_CK119_late_organoid/output/4_final_assignment/data/SeuratObject.rds",
	      "JX1_HC_organoid"="../Individual_analysis_JX1_HC_organoid/output/4_final_assignment/data/SeuratObject.rds",
	      "JX2_PKD1KO_organoid"="../Individual_analysis_JX2_PKD1KO_organoid/output/4_final_assignment/data/SeuratObject.rds",
	      "JX3_PKD2KO_organoid"="../Individual_analysis_JX3_PKD2KO_organoid/output/4_final_assignment/data/SeuratObject.rds")

## Arrange plots
list_dummy <- setNames(vector("list", length=length(list_rds)),
		       names(list_rds))

list_ORIumap <- list_QUEumap <- list_nxn <- list_cnt <- list_hp <- list_dummy

for(idx in names(list_rds)) {
	ORG <- readRDS(list_rds[idx])
	ann <- read.table(file=paste0("./01_symphony_mapping/",idx,".csv"), 
			  sep=",", header=TRUE, stringsAsFactors=FALSE, row.names=1) 
	ORG <- AddMetaData(ORG, ann)
	# compare annotations
	nxn <- table(ind=Idents(ORG), pred=ORG$ann)

	# Record annotation counts from prediction
	cnt <- table(ORG$ann)

	# Arrange visualizations
	# - Heatmap
	hp <- Heatmap(nxn, column_title=idx,
		      column_names_side="top", 
		      row_names_side="left",
		      name="Ncells")
	# - Original umap
	oriUMAP <- DimPlot(ORG) + ggtitle(Project(ORG)) + 
		scale_color_manual(values=col_org[intersect(names(col_org), levels(ORG))])

	# - ... Mapped into human clusters
	queUMAP <- DimPlot(ORG, group.by="ann") + ggtitle(Project(ORG)) + 
		scale_color_manual(values=col_map[intersect(names(col_map), unique(ORG$ann))])

	# Store results
	list_nxn[[idx]] <- nxn 
	list_cnt[[idx]] <- cnt 
	list_QUEumap[[idx]] <- queUMAP 
	list_ORIumap[[idx]] <- oriUMAP
	list_hp[[idx]] <- hp
}


# Build up the grid of plots
sbs <- plot_grid(plot_grid(plotlist=list_ORIumap, nrow=1, align="v"),
	  plot_grid(plotlist=list_QUEumap, nrow=1, align="v"),
	plot_grid(plotlist=lapply(list_hp, function(z) grid.grabExpr(draw(z))), nrow=1, align="v"),
	  nrow=3)

pdf(file=paste0(OUTDIR,"Symphony_sidebyside.pdf"), width=26, height=12)
print(sbs)
dev.off()

# Record annotation counts
uni_ann <- sort(unique(unlist(lapply(list_cnt, function(z) names(z)))))
mat_cnt <- lapply(list_cnt, function(z) {
			  res <- setNames(z[uni_ann], uni_ann)
			  res[is.na(res)] <- 0
			return(res)
	      })
mat_cnt <- do.call("cbind", mat_cnt)

write.table(mat_cnt, file=paste0(OUTDIR,"Symphony_OrganoidMapping_counts.tsv"),
	    sep=",", row.names=TRUE, col.names=NA, quote=FALSE)
