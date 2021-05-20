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

#' Save a active Idents from SeuratObject
#'
#' @param SeuratObject a seurat object
#' @param fl a path to a file
saveActiveIdents <- function(SeuratObject, fl) {
	if(grepl("\\.csv$", fl)) {
		sepx <- ","
	} else if (grepl("\\.tsv$", fl)) {
		sepx <- "\t"
	} else {
		stop("[ERROR] : file extension not supported\n")
	}

	write.table(data.frame("Ident"=SeuratObject@active.ident),
		    file=fl,
		    sep=sepx, col.names = NA, 
		    row.names = TRUE, quote=TRUE)
}

#' Save a active Idents from SeuratObject
#'
#' @param SeuratObject a seurat object
#' @param assay the name of assay
#' @param fl a path to a file
saveClusteringOutcome <- function(SeuratObject, assay="RNA", fl) {

	if(grepl("\\.csv$", fl)) {
		sepx <- ","
	} else if (grepl("\\.tsv$", fl)) {
		sepx <- "\t"
	} else {
		stop("[ERROR] : file extension not supported\n")
	}

	meta_idx <- grep(paste0("^(",assay,"_snn_res\\.|",
				"seurat_clusters)"),
                        colnames(SeuratObject@meta.data),
                        value=TRUE)

	stopifnot(length(meta_idx)>0)
	cat(paste0("[WARN] Selected metacols for Clustering outcome:",
		   paste(meta_idx, collapse=","), "\n"), file=stdout())

	tab <- SeuratObject@meta.data[, meta_idx, drop=FALSE]
	write.table(tab,
		    file=fl,
		    sep=sepx, col.names = NA, 
		    row.names = TRUE, quote=TRUE)
}


# Enhanced DoHeatmap function
DoHeatmap2 <- function(SeuratObject, GSC, assay="RNA", res=0.5, 
                       show_hr=TRUE, cols=NULL, width =NULL,name="Expr.",
                       row_names_fontisze=12,
		       ttl=character(0),
                       legend.dir="vertical") {
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
  
  if(is(SeuratObject)[1]=="Seurat") {
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
  } else {
	  stop("ERROR: Input is not a SeuratObject\n")

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
                                                        labels_gp = gpar(fontsize = 14)))
  
  hr <- rowAnnotation(df=data.frame("type"=hr_classes), show_annotation_name = FALSE,
                      col=list("type"=setNames(gg_color_hue2(length(unique(hr_classes))),
                                               unique(hr_classes))),
                      annotation_legend_param= list(legend_height = unit(4, "cm"),
                                                    grid_height = unit(10, "mm"),
                                                    title_gp=gpar(fontsize=16),
                                                     direction = legend.dir, # ncol=4,
                                                    labels_gp = gpar(fontsize = 22)))
  
  f1 <-  circlize::colorRamp2(c(-2,0,+2), c("purple", "black", "yellow"))
  hp <- Heatmap(mat2, cluster_rows = FALSE, cluster_columns = FALSE,col = f1,
                name=name,
		use_raster=FALSE,
                top_annotation = hc, 
  		bottom_annotation = hc,
		column_title=ttl,
                split=factor(genes.cols, levels=unique(genes.cols)),row_title_rot = 0,row_gap = unit(2, "mm"), 
                column_split = factor(cl, levels=unique(cl)), column_title_rot=ifelse(cl_num,0,90), column_gap = unit(0, "mm"),
                # left_annotation = rowAnnotation(foo=anno_block(gpar(fill=table(genes.cols)[unique(genes.cols)]),
                #                                                labels=unique(genes.cols),
                #                                                labels_gp=gpar(col="white",fontsize=10))),
                width=width,
		border=TRUE,
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

GsHeatmap <- function(sg, 
		      feature="specScore", maxcol="green4", 
		      GSC, cl, cl_num, legend.dir="vertical",
		      row_names_fontisze=12
		      ) {


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
  
  	hr_classes <- genes.cols

	gene_list <- geneIds(GSC)

	# Reorder
	if(cl_num) {
		ord <- order(as.numeric(cl), decreasing = FALSE)
	} else {
		ord <- order(as.character(cl), decreasing = FALSE)
	}
	mat <- as.matrix(sg[[feature]])[unlist(gene_list), ord]
	cl <- cl[ord]

    	cl.cols <- gg_color_hue(length(unique(cl)))
	if(cl_num) {
		names(cl.cols) <- unique(as.character(sort(as.numeric(cl))))
	} else {
		names(cl.cols) <- unique(as.character(sort(cl)))
	} 


	hc <- HeatmapAnnotation(df=data.frame("cluster"=cl),
				col = list("cluster"=cl.cols), show_annotation_name = FALSE,
                          show_legend = FALSE,
			  height=unit(2, "mm"),
                          annotation_legend_param= list(legend_height = unit(8, "cm"),
                                                        grid_width = unit(5, "mm"),
                                                        title_gp=gpar(fontsize=16),
                                                        labels_gp = gpar(fontsize = 14)))
 
	hr <- rowAnnotation(df=data.frame("type"=hr_classes), show_annotation_name = FALSE,
                      col=list("type"=setNames(gg_color_hue2(length(unique(hr_classes))),
                                               unique(hr_classes))),
                      annotation_legend_param= list(legend_height = unit(4, "cm"),
                                                    grid_height = unit(10, "mm"),
                                                    title_gp=gpar(fontsize=16),
                                                     direction = legend.dir, # ncol=4,
                                                    labels_gp = gpar(fontsize = 22)))

	hp <- Heatmap(mat, name=feature,
		col=c("white",maxcol),row_names_gp = gpar(fontsize=10),
		column_names_side = "top",
                top_annotation = hc, 
  		bottom_annotation = hc,
		cluster_rows = FALSE, cluster_columns = FALSE,
		width=unit(6, "cm"),
		split = factor(unlist(sapply(names(gene_list), function(z) rep(z,length(gene_list[[z]])))),
                       levels=names(gene_list)),
                column_split = factor(cl, levels=unique(cl)), column_title_rot=ifelse(cl_num,0,90), column_gap = unit(0, "mm"),
		show_column_names = FALSE,
        row_title_rot = 0,row_gap = unit(2, "mm"),border=TRUE
        )

	return(hp)
}


#' DoMultiHeatmap multi-heatmap for marker confirmation
#' 
#' @param SeuratObject A seurat object with set Idents, scaled data slot with gene expression
#' @param sg a list output from genesorteR with sparse matrices condGeneProb and specScore
#' @param GSC A GeneSetCollection from GSEABase library
#' @param assay assay slot to take from gene expression
#' @param show_plot whether to show the heatmaps or not 
DoMultiHeatmap <- function(SeuratObject, sg, GSC, assay="RNA", show_plot=TRUE) {
	# Row-Scaled Heatmap expression
	hp1 <- DoHeatmap2(SeuratObject, 
			  GS=GSC, show_hr=FALSE,
			  width=unit(16, "cm"),
			  assay = assay, res = NULL, 
			  )
	# Absolute Heatmap of condGeneProb from genesorteR
	hp2 <- GsHeatmap(sg, feature = "condGeneProb", maxcol = "red",
			 GSC = GSC, cl = levels(SeuratObject), cl_num = FALSE)

	# Absolute Heatmap of condGeneProb from genesorteR
	hp3 <- GsHeatmap(sg, feature = "specScore", maxcol = "green4",
			 GSC = GSC, cl = levels(SeuratObject), cl_num = FALSE)

	# Buildup the list of heatmaps
	ht_list <- hp1 + hp2 + hp3

	if(show_plot) {	
		draw(ht_list, column_title=Project(SeuratObject), 
		     column_title_gp = gpar(fontsize = 30))
		return(NULL)
	} else {
		return(ht_list)
	}
}

#' Dotplot2 adaptation of Seurat::DotPlot for gene sets
#'
#' @param SeuratObject a seurat object
#' @param GSC a GeneSetCollection as defined by GSEABase library
#' @param assay the assay to be defined in the dotplot
#' @param res resolution parameter encoded to set clusters
DotPlot2 <- function(SeuratObject, GSC, assay="RNA", res=NULL) {
 	
	# Extract genes
  	genes <- unlist(geneIds(GSC))
	gene_list <- geneIds(GSC)
	if(is(SeuratObject)[1]=="Seurat") {
	  if(all(dim(SeuratObject@assays$RNA@scale.data) == c(0,0))) {
		  SeuratObject <- ScaleData(SeuratObject, verbose = FALSE)
	  }
	  if(is.null(res)) {
		  cl <- as.character(levels(SeuratObject))
		  cl_num <- FALSE
	  } else {
		  cl <- as.character(levels(SeuratObject@meta.data[,paste0(assay,"_snn_res.",res)]))
		  cl_num <- TRUE
	  }
	 } else {
	  stop("ERROR: Input is not a SeuratObject\n")
	}
  

	# Reorder of clusters
	if(cl_num) {
		ord <- order(as.numeric(cl), decreasing = FALSE)
	} else {
		ord <- order(as.character(cl), decreasing = FALSE)
	}

# Cluster order
	Idents(SeuratObject) <- factor(Idents(SeuratObject), levels=levels(SeuratObject)[ord])

	yfeatures <- data.frame("gene"=unlist(gene_list),
				"geneset"=unlist(sapply(names(gene_list), function(z) rep(z, length(gene_list[[z]]))))
				)
	yfeatures$geneset <- factor(yfeatures$geneset, levels=rev(names(GSC)))
	yfeatures <- yfeatures[order(yfeatures$geneset),]
	col_pairs <- setNames(as.integer(factor(yfeatures$geneset)),
			      yfeatures$gene)


	gg <- DotPlot(SeuratObject, features=genes) + 
 		scale_y_discrete(position="right") + 
 		scale_x_discrete(position="bottom") + 
		coord_flip(clip="off") + 
 		scale_size(limits=c(0,100), breaks=seq(0,100,25)) +
		guides(colour=guide_colourbar(title = "Scaled Average Expression",
					      title.position="top", title.hjust=0.5, 
					      barwidth=10),
		       size = guide_legend(title = "Percent Expressed",
					   title.position="top", title.hjust=0.5, 
					   order=2)) +
		theme(axis.title = element_blank(),
# 		      axis.ticks.y = element_line(size=7, colour=col_pairs),
		      axis.text.x = element_text(angle = 45, hjust=0),
		      legend.position="bottom",
# 		      legend.box="vertical",
		      plot.margin = unit(c(2,3,1,10), "lines"))

	yinit <- 0
	for(ann in levels(yfeatures$geneset)) {
		cnt <- sum(yfeatures$geneset==ann)
		ycoord <- mean(c(yinit, yinit + cnt))+0.5
		yinit <- yinit + cnt
		gg <- gg + annotation_custom(grob=textGrob(label=ann, 
							   hjust=1, 
							   gp=gpar(cex=1.5)), 
			       xmin = ycoord, xmax=ycoord,
			       ymin=-1.5, ymax=-1.5)
		if(ann != levels(yfeatures$geneset)[length(levels(yfeatures$geneset))]) {
			gg <- gg + geom_vline(xintercept=yinit+0.5, 
					      alpha = 0.7, colour= "black") # color=col_pairs[yinit])
		}
	}
	return(gg)

}
 

#' GsDotPlot dotplot for genesoteR ouput grouped by gene sets
#'
#' @param sg a list output from genesoteR with sparse matrices condGeneProb and specScore
#' @param GSC A GeneSetCollection from GSEABase library
#' @param cl clusters identities
#' @param cl_num whether cl is numeric or character
GsDotPlot <- function(sg, GSC, cl, cl_num) {

	# Unfortunately we have to revert both gene sets and order of collection to keep
	# them in the right order as in the gene set since ggplot2 reverse the order of
	# y axis when factors
	GSC <- GeneSetCollection(lapply(GSC, function(GS) {
			      GS@geneIds <- rev(GS@geneIds)
			      return(GS)
	    }))
	GSC <- GSC[rev(names(GSC))]
	
  	genes <- unlist(geneIds(GSC))
	gene_list <- geneIds(GSC)
	# Reorder
	if(cl_num) {
		ord <- order(as.numeric(cl), decreasing = FALSE)
	} else {
		ord <- order(as.character(cl), decreasing = FALSE)
	}
	# Cluster order
	cl <- cl[ord]

	# condGeneProb
	cGP <- reshape2::melt(as.matrix(sg[["condGeneProb"]])[unlist(gene_list), ord])
	colnames(cGP) <- c("gene", "cluster", "condGeneProb")
	cGP$cluster <- factor(as.character(cGP$cluster), levels=cl)

	# specScore
	sS <- reshape2::melt(as.matrix(sg[["specScore"]])[unlist(gene_list), ord])
	colnames(sS) <- c("gene", "cluster", "specScore")
	sS$cluster <- factor(as.character(sS$cluster), levels=cl)

 	stopifnot(all(cGP$cluster==sS$cluster))
	stopifnot(all(cGP$gene==sS$gene))
	dat <- cGP
	dat$specScore <- sS$specScore
	rm(cGP, sS)
	dat$geneset <- sapply(as.character(dat$gene), function(gene) {
			idx <- unlist(lapply(geneIds(GSC), function(gs) any(gene%in%gs)))
			res <- paste0(names(GSC)[idx],collapse=",")
			return(res)
		      })
	dat$geneset <- factor(dat$geneset, levels=names(GSC))

	yfeatures <- unique(dat[, c("gene", "geneset")])
	col_pairs <- setNames(as.integer(factor(yfeatures$geneset)),
			      yfeatures$gene)

	gg <- ggplot(dat, aes(x=cluster, y=gene, size=condGeneProb, color=specScore)) +
		geom_point() + theme_cowplot() +
		scale_color_gradientn(colours= c("lightgrey","lightgreen","red"),
				      values= c(0, 0.05, 0.25, 0.35, 1)) + 
#  		scale_color_gradientn(colours= c("lightgrey","blue"),
#  				      values= c(0, 0.05, 0.35, 1)) + 
 		scale_size(limits=c(0,1), breaks=seq(0,1,0.25)) +
		scale_x_discrete(position="top") + 
		coord_cartesian(clip="off") + 
		guides(colour=guide_colourbar(title="Specificity Score",
					      title.position="top", title.hjust=0.5, barwidth=10),
		       size = guide_legend(title="Conditional Gene Probability",
					   title.position="top", title.hjust=0.5, order=2)) +
		theme(axis.title = element_blank(),
# 		      axis.ticks.y = element_line(size=7,colour=col_pairs),
		      axis.text.x = element_text(angle = 45, hjust=0),
		      legend.position="bottom",
# 		      legend.box="vertical",
		      plot.margin = unit(c(2,3,1,10), "lines"))

	yinit <- 0
	for(ann in levels(yfeatures$geneset)) {
		cnt <- sum(yfeatures$geneset==ann)
		ycoord <- mean(c(yinit, yinit + cnt))+0.5
		yinit <- yinit + cnt
		gg <- gg + annotation_custom(grob=textGrob(label=ann, 
							   hjust=1, 
							   gp=gpar(cex=1.5)), 
			       ymin = ycoord, ymax=ycoord,
			       xmin=-1.5, xmax=-1.5)
		if(ann != levels(yfeatures$geneset)[length(levels(yfeatures$geneset))]) {
			gg <- gg + geom_hline(yintercept=yinit+0.5, 
					      alpha = 0.7, colour= "black") # color=col_pairs[yinit])
		}
	}

	return(gg)
}

#' DoMultiDotPlot for expression and genesorter
#'
#' @param SeuratObject a seurat object with idents active
#' @param sg a list from genesorter output
#' @param GSC a gene set collection with marker genes per set
#' @param assay name of assay for expression in the seurat object
#' @param show_NCells boolean whether to show a table with Ncells per cluster
DoMultiDotPlot <- function(SeuratObject, sg, GSC, assay="RNA", show_NCells=FALSE) {
	require(ggplot2)

	#NOTE: implement different resolution than null
	res <- NULL
	plot1 <- DotPlot2(SeuratObject, GSC=GSC, assay="RNA", res=NULL) 
	plot2 <- GsDotPlot(sg, GSC = GSC, cl = levels(SeuratObject), cl_num = FALSE)

	title <- ggdraw() +
		draw_label(Project(SeuratObject), 
			   fontface="bold", size=22, x=0.5, vjust=0.5, hjust=0) +
		theme(plot.margin=unit(c(0, 0, 0, 4), "mm"))
	if(show_NCells) {
		require(gridExtra)
		cnt <- table(Idents(SeuratObject))
		cnt <- unclass(cnt)
		cnt <- cnt[sort(names(cnt))]
		cnt <- as.matrix(cnt)
		title <- title + annotation_custom(tableGrob(t(cnt)), ymin=-0.8, ymax=1)
	}

	plot_row <- plot_grid(plot1, plot2)
	plot_res <- plot_grid(title, plot_row, ncol=1, rel_heights=c(0.1,1))

	return(plot_res)
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

#' DoClustReport make a report of clusters
#' 
#' @param SeuratObject
#' @param sg Output list from GeneSorteR
#' @param GSC a GeneSetCollection class from GSEABase
#' @param show_NCells whether show Ncells per cluster
DoClustReport <- function(SeuratObject, sg, GSC, show_NCells = TRUE) {
	QC_feats <- c("percent.mt", "nFeature_RNA", 
			"S.Score", "G2M.Score", 
			"Dissociation")
	QC_feats <- intersect(QC_feats, colnames(SeuratObject@meta.data))

	A <- DoMultiDotPlot(SeuratObject, sg, GSC = GSC, show_NCells = show_NCells)
	B <- plot_grid(plotlist = sapply(QC_feats, function(feat) {
			       plot1 <- VlnPlot(SeuratObject, group.by="init_assign", 
						features=feat, pt.size=0.2)
			       plot1 <- plot1 + theme(axis.title.x=element_blank(),
						      axis.text.x=element_text(size=10),
						      plot.margin=unit(c(1,1,1,5),"mm")
						      ) + NoLegend()
			       return(plot1)
			}, simplify = FALSE), ncol=length(QC_feats))

	plot_res <- plot_grid(A,B, ncol=1, rel_heights = c(0.8,0.2))
	return(plot_res)
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
			    remaining <- setdiff(colnames(up[[idx]]), firstgofirst)
			    col_names <- c(firstgofirst, remaining)
		writeData(wb, sheet = paste0("wilcox_cluster",idx), 
			  x=up[[idx]][,c(col_names)], 
			  rowNames=FALSE)
	}

	saveWorkbook(wb, file =  xlsx_file, overwrite = TRUE)

}
