### TITLE : Seurat-attached R functions
### AUTHOR : Perales-Paton, Javier - javier.perales@bioquant.uni-heidelberg.de
### DESCRIPTION : These functions are pretty similar to those used in the tutorial
###       but just a wrapper to facilitate readibility of the scripts.

# Get a SeuratObject from 10x data
# Also calculates the percentage of mitocondrial genes
getSeuratObject <- function(path, project_name, min.cells=3, min.features=200, mt.pattern="^MT-") {
  mat10x <-  Read10X(data.dir = path)
  Seurat.obj <- CreateSeuratObject(counts = mat10x, project = project_name, min.cells = min.cells, min.features = min.features)
  
  # Add mitocondrial genes
  Seurat.obj[["percent.mt"]] <- PercentageFeatureSet(Seurat.obj, pattern = mt.pattern)
  
  # Return object
  return(Seurat.obj)
}

# Save the typical violin plot for the lime QC from Seurat.
QCVln <- function(SeuratObject, outdir="./figs/QC") {
  png(paste0(outdir,"/QCvln_",Project(SeuratObject),".png"), width = 1600*3, height = 800*3, res=280)
  print(VlnPlot(SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
}

# Save the Scatter plot to decide thresholds to remove doublets and death cells (mitochondrial genes)
# before subsetting
QCscatter <- function(SeuratObject, outdir="./figs/QC") {
  png(paste0(outdir,"/QCscatter_",Project(SeuratObject),".png"), width = 1600*3, height = 800*3, res=280)
  print(CombinePlots(plots = list(FeatureScatter(SeuratObject, feature1 = "nCount_RNA", feature2 = "percent.mt"),
                                  FeatureScatter(SeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")))
  )
  dev.off()
}

# Scale data for all genes
Seurat_scaledata <- function(SeuratObject) {
  all_genes <- rownames(SeuratObject)
  SeuratObject <- ScaleData(SeuratObject, features = all_genes)
  
  return(SeuratObject)
}

# Show clustertree given a set of set of resolutions.
saveClusterTree <- function(SeuratObject, outdir="./figs/Cluster", prefix="RNA_snn_res.") {
  png(paste0(outdir,"/clustree_",Project(SeuratObject),".png"), width = 2000*3, height = 1000*3, res=280)
  print(clustree(SeuratObject, prefix= prefix))
  dev.off()
}

VarFeatPlot <- function(SeuratObject, Ngenes, outdir="./figs/Cluster") {
  png(filename = paste0(outdir,"/VariableGenes_",Project(SeuratObject),".png"), width = 800*3, height = 800*3, res=280)
  print(LabelPoints(VariableFeaturePlot(SeuratObject), points=VariableFeatures(S)[1:Ngenes],repel = TRUE))
  dev.off()
}

# Plot ModuleScore given a set of genes
PlotModuleScore <- function(SeuratObject, genes, reduction="tsne", tag=NULL) {
  cat("[INFO] : Intersect space with set of genes:\n", file = stdout())
  cat(paste0("\t #",sum(genes %in% rownames(SeuratObject))," out of ",length(genes),"\n",file=stdout()))
  isec <- intersect(rownames(SeuratObject), genes)
  if(length(isec)==0) {
    cat("[ERROR] : No genes are overlapping.\n",file = stdout())
  } else {
  
    SeuratObject <- AddModuleScore(SeuratObject, features=list(c(isec)), name = tag)
    # if(!is.null(tag)) {
    #   colnames(SeuratObject@meta.data)[which(colnames(tmp@meta.data)=="Cluster1")] <- tag
    # } else {
    #   tag <- "Cluster1"
    # }
    print(FeaturePlot(SeuratObject, features=grep(tag,colnames(SeuratObject@meta.data),value=TRUE), reduction=reduction))
  }
  
}

# Enhanced DoHeatmap function
DoHeatmap2 <- function(SeuratObject, GSC, assay="RNA", res=0.5, 
                       show_hr=TRUE, cols=NULL, width =NULL,name="Expr.",
                       row_names_fontisze=12,
                       legend.dir="horizontal") {
  library(ComplexHeatmap)
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  gg_color_hue2 <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 35, c = 100)[1:n]
  }
  
  
  genes <- unlist(geneIds(GSC))
  genes.cols <- unlist(sapply(names(GSC), function(abbn) rep(abbn, length(geneIds(GSC)[[abbn]]))))
  genes.cols <- gsub("__","\n",genes.cols)
  
  if(all(dim(SeuratObject@assays$RNA@scale.data) == c(0,0))) {
    SeuratObject <- ScaleData(SeuratObject, verbose = FALSE)
  }
  mat <- SeuratObject@assays[[assay]]@scale.data
  
  if(is.null(res)) {
    cl <- as.character(Idents(SeuratObject))
    cl_num <- FALSE
  } else {
    cl <- as.character(SeuratObject@meta.data[,paste0(assay,"_snn_res.",res)]) 
    cl_num <- TRUE
  }
  
  # Reorder
  if(cl_num) {
	ord <- order(as.numeric(cl), decreasing = FALSE)
  } else {
	ord <- order(cl, decreasing = FALSE)
  }
  mat <- mat[,ord]
  cl <- cl[ord]
  
  if(is.null(cols)) {
    cl.cols <- gg_color_hue(length(unique(cl)))
  } else {
    cl.cols <- cols[unique(cl)]
  }
  
  if(is.null(names(cl.cols))) {
    if(cl_num) {
      names(cl.cols) <- unique(as.character(sort(as.numeric(cl))))
    } else {
      names(cl.cols) <- unique(as.character(sort(cl)))
    } 
  }
  
  hr_classes <- genes.cols
  
  
  common_genes <- intersect(genes,rownames(mat))
  diff_genes <- setdiff(genes, rownames(mat))
  
  mat2 <- rbind(mat[common_genes,],
                matrix(NA, nrow=length(diff_genes), ncol=ncol(mat),
                       dimnames=list(diff_genes, colnames(mat)))
  )
  mat2 <- mat2[genes,]
  
  hc <- HeatmapAnnotation(df=data.frame("cluster"=cl),col = list("cluster"=cl.cols), show_annotation_name = FALSE,
                          show_legend = FALSE,
			  height=unit(2, "mm"),
                          annotation_legend_param= list(legend_height = unit(8, "cm"),
                                                        grid_width = unit(5, "mm"),
                                                        title_gp=gpar(fontsize=16),
                                                        # at=c(-2.5,-2,-1,0,1,2,2.5),
                                                        # labels = c("","-2","-1","0","1","2",""),
                                                        labels_gp = gpar(fontsize = 14)))
  
  hr <- rowAnnotation(df=data.frame("type"=hr_classes), show_annotation_name = FALSE,
                      col=list("type"=setNames(gg_color_hue2(length(unique(hr_classes))),
                                               unique(hr_classes))),
                      annotation_legend_param= list(legend_height = unit(4, "cm"),
                                                    grid_height = unit(10, "mm"),
                                                    title_gp=gpar(fontsize=16),
                                                    direction = legend.dir, ncol=4,
                                                    # at=c(-2.5,-2,-1,0,1,2,2.5),
                                                    # labels = c("","-2","-1","0","1","2",""),
                                                    labels_gp = gpar(fontsize = 22)))
  
  f1 <-  circlize::colorRamp2(c(-2,0,+2), c("purple", "black", "yellow"))
  hp <- Heatmap(mat2, cluster_rows = FALSE, cluster_columns = FALSE,col = f1,
                name=name,
		use_raster=FALSE,
                top_annotation = hc, 
#  		bottom_annotation = hc,
                # split=factor(genes.cols, levels=unique(genes.cols)),row_title_rot = 0,row_gap = unit(0, "mm"), 
                column_split = factor(cl, levels=unique(cl)), column_title_rot=ifelse(cl_num,0,90), column_gap = unit(0, "mm"),
                # left_annotation = rowAnnotation(foo=anno_block(gpar(fill=table(genes.cols)[unique(genes.cols)]),
                #                                                labels=unique(genes.cols),
                #                                                labels_gp=gpar(col="white",fontsize=10))),
                width=width,
                heatmap_legend_param= list(legend_height = unit(4, "cm"),
                                           title_gp=gpar(fontsize=16),
                                           direction=legend.dir,
                                           # at=c(-2.5,-2,-1,0,1,2,2.5),
                                           # labels = c("","-2","-1","0","1","2",""),
                                           labels_gp = gpar(fontsize = 18)),
                show_column_names = FALSE, row_names_side = "left",
                row_names_gp = gpar(fontsize=row_names_fontisze))
  if(show_hr) {
    hh <- hp + hr 
  } else {
    hh <- hp
  }
  return(hh)
}


## DotPlot function to visualize across samples the average expression and percentage
## of cells per cluster from a set of samples. In contrast to default Seurat function,
## this function provides a list of Seurat-like DotPlots for CombinePlots, so those
## are shown as individual panels.
DotPlot_panel <- function (object=SeuratObject, assay = NULL, features,
                           cols = c("lightgrey", "blue"),
                           col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                           group.by=NULL, split.by = "orig.ident", color.by="avg.exp.scaled",
                           scale.by = "radius", scale.min = NA, scale.max = NA) {
  
  if(!is.null(assay)) {
    if(DefaultAssay(object)!=assay) DefaultAssay(object = object) <- assay
  } else {
    assay <- DefaultAssay(object = object)
  }
  
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  data.features <- FetchData(object = object, vars = features)
  
  # Cell subpopulation groups
  data.features$subpop <- if (is.null(x = group.by)) {
    Idents(object = object)
  } else {
    object[[group.by, drop = TRUE]]
  }
  if (!is.factor(x = data.features$subpop)) {
    data.features$subpop <- factor(x = data.features$subpop)
  }
  data.features$subpop <- as.vector(x = data.features$subpop)
  
  # Grouping for panels
  data.features$split <- object[[split.by, drop = TRUE]]
  
  # Unique ids for subpop per sample
  data.features$id <- paste0(data.features$split,"::",data.features$subpop)
  id.levels <- levels(x = data.features$id)
  
  # Create tidy data.frame for ggplot2
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
                              # 3 are the number of metadata info
                              1:(ncol(x = data.features) - 3), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot ==
                                                     x, "avg.exp"]
                             data.use <- scale(x = data.use)
                             data.use <- MinMax(data = data.use, min = col.min,
                                                max = col.max)
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  
  data.plot$avg.exp.scaled <- avg.exp.scaled
  
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = rev(x = features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  # color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
  #                    no = "colors")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  
  # Split in panels
  data.plot$sname <- sapply(data.plot$id, function(z) strsplit(z, split="::")[[1]][1]) 
  data.plot$id <- sapply(data.plot$id, function(z) strsplit(z, split="::")[[1]][2])
  
  dataset.plot <- split(data.plot, data.plot$sname)
  
  # Create plots
  snames <- unique(object[[split.by, drop=TRUE]])
  plots <- setNames(vector("list",length=length(snames)),
                    snames)
  for(sname in snames) {
    data.plot <- dataset.plot[[sname]]
    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                          y = "id")) + 
      geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
      scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
      labs(x = "Features", y = "Identity") + theme_cowplot()
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2]) + 
      guides(color = guide_colorbar(title = "Average Expression\n(Scaled)")) +
      ggtitle(sname)
    
    plots[[sname]] <- plot
  }
  
  return(plots)
}

#' Save markers as individual CSV files 
#' 
#' @param x a list from loop-wilcox or genesorteR::sortGenes()
#' @param OUTDIR the output directory path
#' @examples
#' saveMarkers.CSV(x, "/output")
saveMarkers.CSV <- function(x, OUTDIR) {
	## Folder to store data
	MARKERS_OUTDIR <- paste0(OUTDIR,"/Markers")
	if(! dir.exists(MARKERS_OUTDIR)) dir.create(MARKERS_OUTDIR, recursive = TRUE)

	if(is.list(x)) {
		if(all(c("condGeneProb", "postClustProb", "specScore") %in% names(x))) {
			# Genesorter
			write.table(as.matrix(x$specScore), paste0(MARKERS_OUTDIR,"/specScore.tsv"),
				    sep="\t",row.names = TRUE,col.names = NA,quote=FALSE)

			write.table(as.matrix(x$condGeneProb), paste0(MARKERS_OUTDIR,"/condGeneProb.tsv"),
				    sep="\t",row.names = TRUE,col.names = NA,quote=FALSE)
		} else {
			# List of DEGs via wilcox
		    for(idx in names(x)) {
			    firstgofirst <- c("cluster", "gene")
			    remaining <- setdiff(colnames(x[[idx]]), firstgofirst)
			    col_names <- c(firstgofirst, remaining)
			write.table(x[[idx]][,col_names],
				    file = paste0(MARKERS_OUTDIR,"/cluster",idx,".tsv"),
		    		    sep="\t",col.names = TRUE, row.names = FALSE, quote=FALSE)
		    }
		}
	} else {
		stop("ERROR: Only list is allowed\n")
	}
}


#' Save markers as in Excel file
#'
#' @param up list with each FindMarkers() each
#' @param sg list object from genesorteR
#' @param SeuratObject with IDents() and $seurat_clusters
#' @param OUTDIR output directory path
#' @examples
#' saveMarkers.Excel(up, sg, SeuratObject, OUTDIR)

saveMarkers.Excel <- function(up, sg, SeuratObject, OUTDIR) {

	## Folder to store data
	MARKERS_OUTDIR <- paste0(OUTDIR,"/Markers")
	if(! dir.exists(MARKERS_OUTDIR)) dir.create(MARKERS_OUTDIR, recursive = TRUE)

	## Excel file
	xlsx_file <- paste0(OUTDIR, "/", Project(SeuratObject), "_markers.xlsx")
	wb <- createWorkbook()

	NclusterNcells <- table(SeuratObject$seurat_clusters)

	IdentClust <- levels(SeuratObject$seurat_clusters)
	IdentSet <- levels(Idents(SeuratObject))

	# 1st sheet
	sheet1 <- data.frame(Cluster=paste0("cluster", IdentClust),
			     Ncells=as.numeric(NclusterNcells),
			     Perc=as.numeric(NclusterNcells)/sum(NclusterNcells)*100,
			     Assignment=rep("?", length(NclusterNcells)),
			     Markers=rep("?", length(NclusterNcells)))
	if(!all(IdentClust==IdentSet)) sheet1$Assignment <- IdentSet;

	addWorksheet(wb, sheetName = "Assignment")
	writeData(wb, sheet ="Assignment", x=sheet1, rowNames=FALSE)

	# GeneSorteR
	specScore <- as.matrix(sg$specScore)
	colnames(specScore) <- paste0("cluster", colnames(specScore))
	addWorksheet(wb, sheetName = "genesorteR_specScore")
	writeData(wb, sheet ="genesorteR_specScore", x=specScore, rowNames=TRUE)

	condGeneProb <- as.matrix(sg$condGeneProb)
	colnames(condGeneProb) <- paste0("cluster", colnames(condGeneProb))
	addWorksheet(wb, sheetName = "genesorteR_condGeneProb")
	writeData(wb, sheet ="genesorteR_condGeneProb", x=condGeneProb, rowNames=TRUE)

	# Wilcox
	for(idx in names(up)) {
		addWorksheet(wb, sheetName = paste0("wilcox_cluster",idx))

			    firstgofirst <- c("cluster", "gene")
			    remaining <- setdiff(colnames(x[[idx]]), firstgofirst)
			    col_names <- c(firstgofirst, remaining)
		writeData(wb, sheet = paste0("wilcox_cluster",idx), 
			  x=up[[idx]][,c("cluster", "gene", col_names)], 
			  rowNames=FALSE)
	}

	saveWorkbook(wb, file =  xlsx_file, overwrite = TRUE)

}
