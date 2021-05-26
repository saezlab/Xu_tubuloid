## Bioinformatics analysis of single-cell and single-nuclei RNAseq data. 
FastQ files were aligned to the human reference genome (GRCh38 assembly) for single-cell 
libraries (scRNA-seq), or a custom reference of the same GRCh38 genome including the introns 
for single-nuclei libraries (snRNA-seq), and the UMI expression was quantified using CellRanger 
(10x Genomics, version 3.1). Raw sequencing data and processed matrices from scRNA-seq were 
deposited in FigShare (https://figshare.com/s/9d06b3df88fb18de653d).

The dataset was analyzed using the Seurat R package (version 3.1.0) [ref+1]. 
Quality control diagnostics were performed by the examination of the total number of detected 
genes per cell, mitochondrial gene expression and library sizes for each sample. For the 
organoid samples (scRNA-seq), cells with less than 200 or greater than 6,000 expressed genes 
were discarded to avoid bad quality cells and cell doublets, respectively. A lenient threshold 
of greater than 20% mitochondrial expression was set to drop bad quality cells, but assuming 
that the threshold retrieves potential proximal tubule-like cells with a characteristic high 
mitochondrial expression. Similar settings were used for the human kidney tissue samples (snRNA-seq), 
with the exception of using more stringent thresholds for doublets (4,000 expressed genes) and 
mitochondrial expression (10%). Raw gene expression was normalized by the library size of each 
cell, and multiplied by a scaling factor of 10,000 afterwards, as it is proposed in the Seurat 
workflow. First, the analysis was done in every sample to find distinct clusters. For this, a 
feature selection prior high dimensionality reduction was performed selecting the 2,000 most 
highly variable genes (HGV) using the ‘vst’ method implemented in Seurat. Mitochondrial genes 
were neglected during this feature selection in snRNA samples. Then, principal component 
analysis was applied to select the first 25 principal components for graph-based clustering 
using Shared-Nearest Neighbour algorithm for cell clustering, and Louvain to find clusters 
(resolution of 0.5). The cell type assignment among the resulting cell clusters was manually 
curated based on cluster-specific markers. For this, gene specificity scores and gene conditional 
probabilities for each cluster were obtained using genesorteR (version 0.3.1) [ref+2]. Additionally, 
differentially over-expressed genes from each cluster were obtained using transcriptome-wide wilcoxon 
rank sum tests implemented in Seurat, but limiting the individual gene tests to only those genes that 
were expressed in 10% of the cells within the cluster with an average over-expression of 0.25 
log2-fold-change. When it is stated in the results, PTPRC and PECAM1 gene expression was used to 
identify contaminant hematopoietic and endothelial cell populations respectively in sorted cells. 
Similarly, the remaining unrelated cell clusters to the biological context under study with high 
content of ribosomal, low transcriptome coverage (200-500 genes) and mitochondrial genes were also 
discarded. When any cluster was discarded in a given sample, a second interaction of unsupervised 
clustering and cluster annotation was performed. Proliferating cell clusters were identified by 
scoring the average expression of cell cycle genes, as compared to a control set as proposed in 
the CellCycleScoring function from Seurat with default parameters [ref+3]. UMAP was used to reduce 
the 25 principal components into 2-dimensional space, showing the cell cluster labels from the 
cell assignment. Dot plots were used to visualize cluster-specific marker expressions, where the 
size of the dot indicates the percentage of cells within the cluster with detected expression, and 
the color shows the average expression of the cluster scaled across the rest of clusters in the 
figure panel.

Gene expression analysis was performed to characterize the polycystic kidney disease (PKD) phenotype 
as compared to controls. For the human tissue samples, PKD and control samples were integrated 
using Harmony (v.1.0) with default parameters and considering every sample library and laboratory 
as an individual batch [ref+4]. The first 20 Harmony-corrected principal components were used 
to find clusters at resolution of 0.5 and calculate the UMAP embedding. Cell clusters overlapping 
cystic populations with controls were identified. Differentially expressed genes were tested by 
pseudo-bulking expression profiles [ref+4], that is summing UMI counts of a given gene from all 
cells from the cluster as a bulk sample. Each subset of pseudo-bulk profiles were processed 
separately by a filtering step of lowly expressed genes, Trimmed Mean of M-values normalization [ref+5]. 
Then the count data was fitted using the negative binomial generalized linear model with 
quasi-likelihood dispersion estimation, and test for differences using the edgeR quasi-likelihood 
pipeline (v.3.26.7) [ref+6]. Differentially expressed genes were considered at a false discovery 
rate of 5% after multiple testing correction. The enrichment of biological pathways were tested 
using fgsea (v.1.0.1) [ref+7], with KEGG, PID and GO terms Biological Processes from MSigDB [ref+8], 
on the ranking of differential expression of each cystic cell population. Footprint-based pathway 
and transcription factor activity were estimated using PROGENy and DoRothEA, respectively [ref+9], 
applied on the ranking of differential expression. For the organoid samples no biological replicates 
were available, thus a more simple and robust exploratory analysis was done. Cell subclusters for a 
given cell type were pooled and compared against (healthy) controls. Seurat implementation of 
Wilcoxon rank sum tests, for up-regulated genes in PKD (adjusted p-value < 0.05, at least 0.25 
log2-fold-change and 10% expression within cluster) was used to explore differences between PKD 
and healthy organoids. Differentially expressed genes were used for over-representation analysis 
of KEGG pathways and GO terms Biological Processes using clusterProfiler (v3.12.0) [ref+10].

Code for the data analysis is available in GitHub at https://github.com/saezlab/Xu_tubuloid .


---
## REFERENCES:
[ref+1] Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM, Hao Y, Stoeckius M, Smibert P, Satija R. Comprehensive Integration of Single-Cell Data. Cell. 2019 Jun 13;177(7):1888-1902.e21

[ref+2] MM Ibrahim and  Kramann R. genesorteR: Feature Ranking in Clustered Single Cell Data. bioRxiv

[ref+3] Tirosh I, Izar B, Prakadan SM, Wadsworth MH, Treacy D, Trombetta JJ, Rotem A, Rodman C, Lian C, Murphy G, Fallahi-Sichani M, Dutton-Regester K, Lin JR, Cohen O, Shah P, Lu D, Genshaft AS, Hughes TK, Ziegler CG, Kazer SW, Gaillard A, Kolb KE, Villani AC, Johannessen CM, Andreev AY, Van Allen EM, Bertagnolli M, Sorger PK, Sullivan RJ, Flaherty KT, Frederick DT, Jané-Valbuena J, Yoon CH, Rozenblatt-Rosen O, Shalek AK, Regev A, Garraway LA. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science. 2016 Apr 8;352(6282):189-96

[ref+4] Tung, P. Y., J. D. Blischak, C. J. Hsiao, D. A. Knowles, J. E. Burnett, J. K. Pritchard, and Y. Gilad. 2017. “Batch effects and the effective design of single-cell gene expression studies.” Sci. Rep. 7 (January): 39921.

[ref+5] Robinson, M.D., Oshlack, A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol 11, R25 (2010). https://doi.org/10.1186/gb-2010-11-3-r25

[ref+6] Chen, Y., A. T. Lun, and G. K. Smyth. 2016. “From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.” F1000Res 5: 1438.

[ref+7] Alexey Sergushichev. An algorithm for fast preranked gene set   enrichment analysis using cumulative statistic calculation. bioRxiv (2016), doi:10.1101/060012

[ref+8] Arthur Liberzon, Aravind Subramanian, Reid Pinchback, Helga Thorvaldsdóttir, Pablo Tamayo, Jill P. Mesirov, Molecular signatures database (MSigDB) 3.0, Bioinformatics, Volume 27, Issue 12, 15 June 2011, Pages 1739–1740, https://doi.org/10.1093/bioinformatics/btr260

[ref+9] Holland, C.H., Tanevski, J., Perales-Patón, J. et al. Robustness and applicability of transcription factor and pathway analysis tools on single-cell RNA-seq data. Genome Biol 21, 36 (2020). https://doi.org/10.1186/s13059-020-1949-z

[ref+10] Yu G, Wang L-G, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.
