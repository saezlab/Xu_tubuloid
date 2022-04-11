#' Tubuloid cells mapped on referenced harmony-integrated human tissue
#' Javier Perales-Patón

#' Following : https://github.com/immunogenomics/symphony/blob/main/vignettes/Seurat.ipynb
#' NOTE on reproduciblity : It runs symphony over every organoid with human kidney tissue (ADPKD+control) as reference.
#'  It implies to rerun UMAP and Harmony integration from the reference as it is expected some extra data from the model,
#'   which is is stored from running a newer version of Harmony.

## Environment reproducibility
set.seed(1234)
OUTDIR <- "./01_symphony_mapping/"
if(!dir.exists(OUTDIR)) dir.create(OUTDIR);

## Load libraries
library(symphony)
library(harmony)
library(Seurat)
library(tibble)
source("utils_seurat.R")

## Load human data
S <- readRDS("../Cystic_signature/01_harmony_integration/data/S.rds")

## Re-Annotation on integrated clusters
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
S$idx <- Idents(S)
S$ann <- ren_id[Idents(S)]

subcl <- read.table("./adpkd_subtypes.csv", sep=",", header=TRUE, stringsAsFactors=FALSE, row.names=1)
stopifnot(subcl$index %in% colnames(S))

S$subcluster <- S$ann 
S$subcluster[subcl$index] <- paste0(gsub("[0-9]+: ", "", subcl$cluster_name), ".", as.character(subcl$subcluster_name))
#table(S$ann, S$subcluster)

# Unfortunately Symphony requires some aspects from the umap run that are not stored in Seurat framework by default.
# Thus UMAP must be rerun, probably yielding a slightly different umap outcome. We keep tack of that.
S[['umap']] <- RunUMAP2(Embeddings(S, 'harmony')[, 1:20], 
                          assay='RNA', verbose=FALSE, umap.method='uwot', return.model=TRUE)
## Make symphony ref object
## Re-run harmony because we miss an expected 'obj@reductions$harmony@misc$R'.
## Probably source: different versions of harmony run in the data.
S2 <- RunHarmony.Seurat(S, c("lab", 'orig.ident'), verbose = FALSE)
## Then we transfer this missing piece of data to previous run
S@reductions$harmony@misc <- S2@reductions$harmony@misc
ref <- buildReferenceFromSeurat(S, verbose=TRUE, 
				save_umap = TRUE, 
				save_uwot_path = 'cache_symphony.uwot')

## Organoid data
list_rds <- c("CK5_organoid"="../Individual_analysis_CK5_early_organoid/output/4_final_assignment/data/SeuratObject.rds",
	      "CK119_organoid"="../Individual_analysis_CK119_late_organoid/output/4_final_assignment/data/SeuratObject.rds",
	      "JX1_HC_organoid"="../Individual_analysis_JX1_HC_organoid/output/4_final_assignment/data/SeuratObject.rds",
	      "JX2_PKD1KO_organoid"="../Individual_analysis_JX2_PKD1KO_organoid/output/4_final_assignment/data/SeuratObject.rds",
	      "JX3_PKD2KO_organoid"="../Individual_analysis_JX3_PKD2KO_organoid/output/4_final_assignment/data/SeuratObject.rds")

list_barcodes <- setNames(vector("list", length=length(list_rds)),
		       names(list_rds))
for(idx in names(list_rds)) {
	ORG <- readRDS(list_rds[idx])
	## Map Query
	set.seed(1234)
	query <- mapQuery(ORG@assays$RNA@counts,
			  ORG@meta.data,
			  ref,
			  vars = "orig.ident",
			  return_type = "Seurat")
	# Predict clusters
	query <- knnPredict.Seurat(query, ref, "ann")
	query <- knnPredict.Seurat(query, ref, "subcluster")
	list_barcodes[[idx]] <- query@meta.data[, c("ann", "subcluster")]
}

# Barcodes and referenced annotation
sapply(names(list_barcodes), function(idx) write.csv(as_tibble(list_barcodes[[idx]], 
							       rownames="barcode"), 
						     file=paste0(OUTDIR,idx,".csv"), 
						     quote=FALSE, row.names=FALSE))
