#' Differential gene expression between ADPKD and control cells using MAST
#' Javier Perales-Patón

## Environment reproducibility
set.seed(1234)

OUTDIR <- "./diffExpr_MAST/"
if(!dir.exists(OUTDIR)) dir.create(OUTDIR);

## Load libraries
library(Seurat)

########################################################################3
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

## Set target populations
target <- ren_id
names(target) <- NULL

# S <- S[, Idents(S) %in% target]

control_gr <- grep("Control", S$orig.ident, value=TRUE)
res_hs <- sapply(target, function(cell) {
			 Si <- S[, Idents(S)==cell]
			cat("[INFO]: Checking for a correct subset of the data\n", file=stdout())
			cat(paste0("[INFO]: Selected cell:",cell,"\n"), file=stdout())
			print(table(Idents(Si), Si$orig.ident))
			 Idents(Si) <- ifelse(grepl("Control", Si$orig.ident), "control", "ADPKD")
			Si$record <- factor(Idents(Si), levels=c("control", "ADPKD"))
			 Si$orig.ident <- "humans"
			if(all(table(Si$record)>2)) {
			 res <- FindMarkers(Si, ident.1="ADPKD", ident.2="control",
					    logfc.threshold=0,
					    min.pct=0,
					    method="MAST")
			} else {
			 res <- data.frame(X="NoTest")
			}
	       		return(res)
	      }, simplify=FALSE)

names(res_hs) <- gsub("/",".",names(res_hs))

# Save it
sapply(names(res_hs), function(cell) {
			      fl <- paste0("MASTde_3xADPKD_",cell,".csv")
	      write.table(res_hs[[cell]],
			  file=paste0(OUTDIR,fl), sep=",", 
			  col.names=NA, row.names=TRUE, quote=FALSE)
	      })

