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
DoHeatmap2 <- function(SeuratObject, GSC, assay="RNA", res=0.5, show_hr=TRUE) {
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
  mat <- SeuratObject@assays[[assay]]@scale.data
  
  if(is.null(res)) {
    cl <- as.character(SeuratObject@meta.data[,"seurat_clusters"])
  } else {
    cl <- as.character(SeuratObject@meta.data[,paste0(assay,"_snn_res.",res)]) 
  }
  
  # Reorder
  ord <- order(as.numeric(cl), decreasing = FALSE)
  mat <- mat[,ord]
  cl <- cl[ord]
  
  cl.cols <- setNames(gg_color_hue(length(unique(cl))),unique(as.character(sort(as.numeric(cl)))))
  
  common_genes <- intersect(genes,rownames(mat))
  diff_genes <- setdiff(genes, rownames(mat))
  
  mat2 <- rbind(mat[common_genes,],
                matrix(NA, nrow=length(diff_genes), ncol=ncol(mat),
                       dimnames=list(diff_genes, colnames(mat)))
  )
  mat2 <- mat2[genes,]
  
  hc <- HeatmapAnnotation(df=data.frame("cluster"=cl),col = list("cluster"=cl.cols), show_annotation_name = FALSE,
                          show_legend = FALSE,
                          annotation_legend_param= list(legend_height = unit(8, "cm"),
                                                        grid_width = unit(5, "mm"),
                                                        title_gp=gpar(fontsize=16),
                                                        # at=c(-2.5,-2,-1,0,1,2,2.5),
                                                        # labels = c("","-2","-1","0","1","2",""),
                                                        labels_gp = gpar(fontsize = 14)))
  
  # hr <- rowAnnotation(df=data.frame(markers=genes.cols), show_annotation_name = FALSE)
  hr_classes <- gsub("\\-?[0-9]+$","",genes.cols)
  hr_classes <- gsub("^Modulated_","",hr_classes)
  hr_classes <- gsub("^Healthy_","",hr_classes)
  hr_classes <- gsub("^(CD8|Mixed|CXCR6)_","",hr_classes)
  hr_classes <- gsub("^(M2_?|ResLike|Inflammatory_?|TREM_high)-","",hr_classes)
  
  
  hr <- rowAnnotation(df=data.frame("type"=hr_classes), show_annotation_name = FALSE,
                      col=list("type"=setNames(gg_color_hue2(length(unique(hr_classes))),
                                               unique(hr_classes))),
                      annotation_legend_param= list(legend_height = unit(4, "cm"),
                                                    grid_width = unit(5, "mm"),
                                                    title_gp=gpar(fontsize=16),
                                                    # at=c(-2.5,-2,-1,0,1,2,2.5),
                                                    # labels = c("","-2","-1","0","1","2",""),
                                                    labels_gp = gpar(fontsize = 14)))
  
  f1 <-  circlize::colorRamp2(c(-2,0,+2), c("purple", "black", "yellow"))
  hp <- Heatmap(mat2, cluster_rows = FALSE, cluster_columns = FALSE,col = f1,
                name="Expression",
                top_annotation = hc, bottom_annotation = hc,
                split=factor(genes.cols, levels=unique(genes.cols)),row_title_rot = 0,row_gap = unit(1.5, "mm"), 
                column_split = factor(cl, levels=unique(cl)), column_title_rot=00, column_gap = unit(0.7, "mm"),
                # left_annotation = rowAnnotation(foo=anno_block(gpar(fill=table(genes.cols)[unique(genes.cols)]),
                #                                                labels=unique(genes.cols),
                #                                                labels_gp=gpar(col="white",fontsize=10))),
                heatmap_legend_param= list(legend_height = unit(4, "cm"),
                                           title_gp=gpar(fontsize=16),
                                           # at=c(-2.5,-2,-1,0,1,2,2.5),
                                           # labels = c("","-2","-1","0","1","2",""),
                                           labels_gp = gpar(fontsize = 15)),
                show_column_names = FALSE, row_names_side = "left",
                row_names_gp = gpar(fontsize=7))
  if(show_hr) {
    hh <- hp + hr 
  } else {
    hh <- hp
  }
  return(hh)
}


# Enhanced DoHeatmap function
DoHeatmap3 <- function(SeuratObject, GSC, assay="RNA", res=0.5, show_hr=TRUE) {
  library(ComplexHeatmap)
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  gg_color_hue2 <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 35, c = 100)[1:n]
  }
  
  gg_color_hue3 <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 95, c = 100)[1:n]
  }
  
  
  genes <- unlist(geneIds(GSC))
  genes.cols <- unlist(sapply(names(GSC), function(abbn) rep(abbn, length(geneIds(GSC)[[abbn]]))))
  mat <- SeuratObject@assays[[assay]]@scale.data
  
  if(is.null(res)) {
    cl <- as.character(SeuratObject@meta.data[,"seurat_clusters"])
  } else {
    cl <- as.character(SeuratObject@meta.data[,paste0(assay,"_snn_res.",res)]) 
  }
  
  origin <- SeuratObject@meta.data$orig.ident
  
  # Reorder
  ord <- order(as.numeric(cl), decreasing = FALSE)
  mat <- mat[,ord]
  
  cl <- cl[ord]
  cl.cols <- setNames(gg_color_hue(length(unique(cl))),unique(as.character(sort(as.numeric(cl)))))
  
  origin <- origin[ord]
  origin.cols <- setNames(gg_color_hue3(length(unique(origin))),
                          unique(as.character(origin)))
  
  common_genes <- intersect(genes,rownames(mat))
  diff_genes <- setdiff(genes, rownames(mat))
  
  mat2 <- rbind(mat[common_genes,],
                matrix(NA, nrow=length(diff_genes), ncol=ncol(mat),
                       dimnames=list(diff_genes, colnames(mat)))
  )
  mat2 <- mat2[genes,]
  
  hc <- HeatmapAnnotation(df=data.frame("cluster"=cl, "origin"=origin),
                          col = list("cluster"=cl.cols, "origin"=origin.cols),
                          show_annotation_name = TRUE,annotation_name_side = "left",
                          show_legend = c(FALSE,TRUE),
                          annotation_legend_param= list(legend_height = unit(8, "cm"),
                                                        grid_width = unit(5, "mm"),
                                                        title_gp=gpar(fontsize=16),
                                                        # at=c(-2.5,-2,-1,0,1,2,2.5),
                                                        # labels = c("","-2","-1","0","1","2",""),
                                                        labels_gp = gpar(fontsize = 14)))
  
  # hr <- rowAnnotation(df=data.frame(markers=genes.cols), show_annotation_name = FALSE)
  hr_classes <- gsub("\\-?[0-9]+$","",genes.cols)
  hr_classes <- gsub("^Modulated_","",hr_classes)
  hr_classes <- gsub("^Healthy_","",hr_classes)
  hr_classes <- gsub("^(CD8|Mixed|CXCR6)_","",hr_classes)
  hr_classes <- gsub("^(M2_?|ResLike|Inflammatory_?|TREM_high)-","",hr_classes)
  
  
  hr <- rowAnnotation(df=data.frame("type"=hr_classes), show_annotation_name = FALSE,
                      col=list("type"=setNames(gg_color_hue2(length(unique(hr_classes))),
                                               unique(hr_classes))),
                      annotation_legend_param= list(legend_height = unit(4, "cm"),
                                                    grid_width = unit(5, "mm"),
                                                    title_gp=gpar(fontsize=16),
                                                    # at=c(-2.5,-2,-1,0,1,2,2.5),
                                                    # labels = c("","-2","-1","0","1","2",""),
                                                    labels_gp = gpar(fontsize = 14)))
  
  f1 <-  circlize::colorRamp2(c(-2,0,+2), c("purple", "black", "yellow"))
  hp <- Heatmap(mat2, cluster_rows = FALSE, cluster_columns = FALSE,col = f1,
                name="Expression",
                top_annotation = hc, bottom_annotation = hc,
                split=factor(genes.cols, levels=unique(genes.cols)),row_title_rot = 0,row_gap = unit(1.5, "mm"), 
                column_split = factor(cl, levels=unique(cl)), column_title_rot=00, column_gap = unit(0.7, "mm"),
                # left_annotation = rowAnnotation(foo=anno_block(gpar(fill=table(genes.cols)[unique(genes.cols)]),
                #                                                labels=unique(genes.cols),
                #                                                labels_gp=gpar(col="white",fontsize=10))),
                heatmap_legend_param= list(legend_height = unit(4, "cm"),
                                           title_gp=gpar(fontsize=16),
                                           # at=c(-2.5,-2,-1,0,1,2,2.5),
                                           # labels = c("","-2","-1","0","1","2",""),
                                           labels_gp = gpar(fontsize = 15)),
                show_column_names = FALSE, row_names_side = "left",
                row_names_gp = gpar(fontsize=7))
  if(show_hr) {
    hh <- hp + hr 
  } else {
    hh <- hp
  }
  return(hh)
}

