#' Differential gene expression analysis using Queried annotation to a reference human map
#' Javier Perales-Patón (c)

## Environment reproducibility
set.seed(1234)
OUTDIR <- "./03_diffExpr/"
if(!dir.exists(OUTDIR)) dir.create(OUTDIR);

## Load libraries
library(Seurat)
library(gridExtra)

## Set target populations
#NOTE: here we identified two populations in the human data that are recalled in organoids,
#	- PT_4
#	- TAL_2
# Thus we focus our analysis on these, and compare each PKD vs late control organoid
target <- c("PT_4", "TAL_2")

## Load data
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
Idents(S) <- S$ann

S <- S[, Idents(S) %in% target]

## Organoid data
list_rds <- c("CK119_organoid"="../Individual_analysis_CK119_late_organoid/output/4_final_assignment/data/SeuratObject.rds",
	      "JX2_PKD1KO_organoid"="../Individual_analysis_JX2_PKD1KO_organoid/output/4_final_assignment/data/SeuratObject.rds",
	      "JX3_PKD2KO_organoid"="../Individual_analysis_JX3_PKD2KO_organoid/output/4_final_assignment/data/SeuratObject.rds")

## Create a list of Seurat objects to merge them
SL_org <- setNames(vector("list", length=length(list_rds)), names(list_rds))
for(idx in names(list_rds)) {
	ORG <- readRDS(list_rds[idx])
	# Load queried annotations
	ann <- read.table(file=paste0("./01_symphony_mapping/",idx,".csv"), 
			  sep=",", header=TRUE, stringsAsFactors=FALSE, row.names=1) 
	ORG <- AddMetaData(ORG, ann)
	Idents(ORG) <- ORG$ann
	ORG <- ORG[, Idents(ORG) %in% target]
	SL_org[[idx]] <- ORG
}

S_org <- merge(SL_org[[1]], SL_org[-1], merge.data=TRUE)

## Define common genes
# intersect to common genes in both sides
common_genes <- intersect(rownames(S_org), rownames(S))
S_org <- S_org[common_genes, ]
S <- S[common_genes, ]

## Testing for differential expression
control_gr <- "CK119_organoid"
# Merging both PKD organoids 
res_org <- sapply(target, function(cell) {
		S <- S_org[, Idents(S_org)==cell]
		cat("[INFO]: Checking for a correct subset of the data\n", file=stdout())
		cat(paste0("[INFO]: Selected cell:",cell,"\n"), file=stdout())
		print(table(Idents(S), S$orig.ident))
		Idents(S) <- ifelse(S$orig.ident == "CK119_organoid", "control", "PKD")
		S$orig.ident <- "threeOrganoids"
		set.seed(1234)
		res <- FindMarkers(S, ident.1="PKD", ident.2="control",
			       logfc.threshold=0,
			       min.pct=0,
			       method="MAST")
		return(res)
	      }, simplify=FALSE)
# Save it
sapply(names(res_org), function(cell) {
		fl <- paste0(OUTDIR,"DiffExpr_PKDorganoid_",cell,".csv")
		write.table(res_org[[cell]],
			  file=fl, sep=",", 
			  col.names=NA, row.names=TRUE, quote=FALSE)
	      })


########################################################################3
control_gr <- grep("Control", S$orig.ident, value=TRUE)
res_hs <- sapply(target, function(cell) {
			 Si <- S[, Idents(S)==cell]
			cat("[INFO]: Checking for a correct subset of the data\n", file=stdout())
			cat(paste0("[INFO]: Selected cell:",cell,"\n"), file=stdout())
			print(table(Idents(Si), Si$orig.ident))
			 Idents(Si) <- ifelse(grepl("Control", Si$orig.ident), "control", "ADPKD")
			 Si$orig.ident <- "humans"
			 set.seed(1234)
			 res <- FindMarkers(Si, ident.1="ADPKD", ident.2="control",
					    logfc.threshold=0,
					    min.pct=0,
					    method="MAST")
	       		return(res)
	      }, simplify=FALSE)

# Save it
sapply(names(res_hs), function(cell) {
			      fl <- paste0(OUTDIR,"DiffExpr_ADPKDhuman_",cell,".csv")
	      write.table(res_hs[[cell]],
			  file=fl, sep=",", 
			  col.names=NA, row.names=TRUE, quote=FALSE)
	      })
