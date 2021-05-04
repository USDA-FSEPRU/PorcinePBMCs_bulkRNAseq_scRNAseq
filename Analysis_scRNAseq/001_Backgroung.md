---
Title:  "Porcine immune cell atlas single cell RNA sequencing analysis pipeline"
Author: "Sathesh K. Sivasankaran"
Date:   "Feb 15, 2021"
---

# Porcine Immune Single Cell Atlas

## Summary
Pigs are an important protein source supporting global food security and are valuable as a human biomedical model. To examine immune cell gene expression and annotate genome function in this important species, single-cell RNA-sequencing (scRNA-seq) of PBMC populations was performed to more highly resolve PBMC population transcriptomes. Across seven PBMC samples partitioned and sequenced, 28,810 cells distributed across 36 clusters were classified into 13 general cell types including plasmacytoid dendritic cells (DC), conventional DCs, monocytes, >10 B cell clusters including plasmablasts, conventional CD4 and CD8 alpha-beta T cells, NK cells, and 4 clusters of gamma-delta T cells. This resource will be invaluable for annotation of pig genes controlling immunogenetic traits, as well as further study of, and development of new reagents for, porcine immunology.

## Sample information
| Internal | Publication |
|----------|-------------|
| CT230OCT | B1          |
| CT21NOV  | B2          |
| A        | C1          |
| B        | C2          |
| C        | C3          |
| PBMC1    | D1          |
| PBMC2    | D2          |

## Materials and Methods - scRNA-seq Data Analysis
### scRNA-seq read alignment and gene quantification
Quality of raw reads were checked by using the software FASTQC (Andrews 2010). Reads 2 (R2) were corrected for errors using Rcorrector (Song and Florea., 2015) and 3’ polyA tails longer than 10 bases were trimmed using custom Perl scripts. After trimming, R2 longer than 25 bases were re-paired using BBMap (v38.36; Bushnell 2014). Sus scrofa genome Sscrofa 11.1 and annotation GTF (v11.1.97) from Ensemble were used to build the reference genome index (Yates et al., 2020). The annotation GTF file was modified to include both gene symbol (if available) and Ensembl ID as a gene reference (e.g. GZMA_ENSSSCG00000016903) using custom Perl scripts. Processed paired end reads were aligned and gene expression count matrices generated using the CellRanger (v4.0; 10X Genomics) with default parameters. Only reads that were confidently mapped (MAPQ=255), non-PCR duplicates with valid barcodes and unique molecular identifiers (UMIs) were counted to generate gene expression count matrix. Reads with the same cell barcodes, the same UMIs, and mapped to same gene features were collapsed into a single read. Subsequent downstream analysis including quality control, dimensionality reduction, cell clustering and differential gene expression was performed within R environment (v3.6.1; https://www.r-project.org/).

### Ambient RNA removal
Ambient RNA in scRNA-seq samples were cleaned using R package SoupX (v1.4.5; Young and Behjati 2018). Briefly, raw count matrices, gene features and t-SNE visualization coordinates from CellRanger output were imported to SoupX. Ambient RNA fractions were estimated independently for each sample using the automated estimation function autoEstCont with default parameters. Then calculated background contamination level was used to correct the gene count matrix and outputted in CellRanger file format using the R package DropletUtils (v1.6.1; Lun et al., 2019). Newly generated files with modified gene counts were used for further analysis.
Non-expressed gene and poor-quality cell filtering
Modified gene count matrices from seven samples were imported to R package, Seurat (v3.2.1; Stuart et al., 2019) and processed using custom R scripts. Gene counts, Cell barcodes and UMIs were corrected and filtered as described in (Zheng et al., 2017). Briefly, genes without any reads detected across all cells and samples were removed from the dataset. Data was checked for cell barcode duplication and removed the from the dataset if any. Cells with greater than 10% of total reads attributed to mitochondrial genes, number of genes detected less than 500 and less than 1000 unique transcripts (UMIs) detected in within each sample were removed from the dataset. Finally, filtered count matrix for each individual sample was generated and outputted into CellRanger file format using DropletUtils R package for cell doublet removal.

### Doublet removal
progressed cells and genes from each sample were analyzed using the Python package, Scrublet (v0.2; Wolock et al., 2019), to remove highly probable neotypic doublets. For each sample, synthetic doublets were generated using default parameters from the filtered count matrices. A Scrublet object was created using an expected doublet rate of 0.07. Doublet scores were calculated, and a doublet score threshold of 0.25 was set based on a bimodal distribution observed between embedded and neotypic doublets. Lists of doublet probability scores were generated and imported into R, where cells with corresponding doublet probability scores grater than 0.25 were removed from the dataset. Filtered count matrix for each sample was outputted into CellRanger format using DropletUtils.

### Sample integration, clustering and differential gene expression
Seurat was used for data transformation, integration, dimensionality reduction and cell clustering for the dataset. In brief, the gene counts of cells from each sample that had survived all quality control steps was loaded into Seurat and transformed individually using SCTransform function. To integrate cells from different samples, integration gene features were selected using SelectIntegrationFeatures function, identified features were verified using PrepSCTIntegration and anchored between Seurat objects using FindIntegrationAnchors function with default parameters. Then the dataset was integrated with IntegrateData function using default parameters. The expression level of highly variable genes in the cells was scaled and centered along each gene, and the principal component analysis was conducted. We then assessed the number of PCs to be included in all further analyses by calculating the cumulative standard deviations accounted for each PC to identify the ‘knee’ point PC number (PC shows less than 0.1% degree of variation with successive PC). The selected significant PCs were used for two-dimensional t-distributed stochastic neighbor embedding (t-SNE) and uniform manifold approximation and projection (UMAP) coordinates generation for visualization.
We used FindCluster function to identify cell clusters with the cluster resolution value of 1.85. To profile gene expression across cell cluster and calculate differentially expressed (DE) genes between cluster pairs, filtered count matrices were log normalized using NormalizeData function with default parameters and log-normalized data was scaled using ScaleData. The functions FeaturePlot and DotPlot in Seurat was utilized to generate expression profile of individual genes across cell clusters. DE genes between each pairwise cluster combination was performed using Seurat function FindMarkers. Genes that expressed at least 20% of cells within one of the cluster pair, with a |logFC| > 0.25, and an adjusted p-value < 0.05 were considered as DE genes.

### Cluster subsetting
To visualize gene expression in specific subsets of cell clusters, specific clusters were isolated from the original dataset using subset, and DietSeurat was performed to remove scaled data. RNA counts were rescaled for further visualization of relative gene expression within only the new subset using ScaleData. UMAP and t-SNE coordinates of only the specific clusters included in subset analysis were reconstructed. To leave the original cluster identities intact, cells were not re-integrated or re-clustered.

### Data availability
The raw fastq files were uploaded to https://www.ebi.ac.uk/ena/browser/view/PRJEB43267.

### Reference
Andrews S (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc.

Song L and Florea L (2015). Rcorrector: efficient and accurate error correction for Illumina RNA-seq reads. Gigascience. 4:48. doi: 10.1186/s13742-015-0089-y.

Bushnell B (2014) BBMap: a fast, accurate, splice-aware aligner. In: 9th Annual Genomics of Energy & Environment Meeting, Walnut Creek, CA. Available online at: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/.

Yates AD et al., (2020), Ensembl 2020. Nucleic Acids Res. 48:D682-D688. doi: 10.1093/nar/gkz966.

Young MD, Behjati S (2020). SoupX removes ambient RNA contamination from droplet based single-cell RNA sequencing data. Gigascience. 9(12):giaa151. doi: 10.1093/gigascience/giaa151.

Lun ATL, Riesenfeld S, Andrews T, Dao T, Gomes T, participants in the 1st Human Cell Atlas Jamboree, Marioni JC (2019). EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. Genome Biol. 20, 63. doi: 10.1186/s13059-019-1662-y.

Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM, Hao Y, Stoeckius M, Smibert P, Satija R (2019). Comprehensive Integration of Single-Cell Data. Cell. 177:1888-1902.e21. doi: 10.1016/j.cell.2019.05.031.

Zheng GX, Terry JM, and Belgrader P (2017). Massively parallel digital transcriptional profiling of single cells. Nat Commun 8. doi: 10.1038/ncomms14049.

Wolock SL, Lopez R and Klein AM (2019). Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. Cell Syst. 8:281-291.e9. doi: 10.1016/j.cels.2018.11.005
