## Analysis_scRNAseq

This directory contains scripts used to perform scRNA-seq analysis from seven porcine PBMC samples.

Selected processed output files/data from analyses described in our scripts can be found in the Ag Data Commons at https://data.nal.usda.gov/dataset/data-reference-transcriptomics-porcine-peripheral-immune-cells-created-through-bulk-and-single-cell-rna-sequencing. Raw data used for analysis can be found at https://www.ebi.ac.uk/ena/browser/view/PRJEB43826.

Materials in this repository include:

| File Name | File Description |
| ----------- | ------------------ |
| 001_Background | Background information for initial processing of data |
| 002_Alignment | Alignment of sequencing data to pig genome |
| 003_AmbientRNAcorrection | Estimation and removal of ambient RNA contamination |
| 004_QCfiltering | Removal of poor quality cells and non-expressed genes |
| 005_DoubletRemoval | Removal of highly probable doublets |
| 006_IntegrationAndClustering | Normalization of reads, integration of samples, and clustering of cells |
| 007_ClusterCharacterization | Further characterizing cell clusters to identify general cell types described in previous porcine literature |
| 008_GeneSetEnrichmentAnalysis | Performing gene set enrichment analyses based on gene set signatures for particular sorted cell populations used for bulk RNA-seq analyses; bulk RNA-seq gene sets are compared to single-cell gene profiles to determine enrichment for each individual cell in our scRNA-seq dataset |
| 009_CD4subset | Revamped analysis of only CD4+ T cells from scRNA-seq data |
| 010_GDsubset | Revamped analysis of only CD4+ T cells from scRNA-seq data |
| 011_AgDataCommonsDeposition | Creating data files to deposit on Ag Data Commons |
| RandomForestModeling | Training of a random forest model to cell cluster identities |
| ReferenceBasedLabelTransfer | Cell label prediction and mapping to an external dataset of human PBMCs |
| MitoGenes | List of mitochondrial genes used to calculate percent mitochondrial reads in cells |
| UnfilteredGeneInfo | List of all genes in the v97 genome annotation that were considered in analyses |

Analyses were performed collaboratively by Jayne Wiarda, Sathesh Sivasankaran, and Lance Daharsh. Some scripts for analyses can also be found originally in personal repositories: https://github.com/jwiarda/scRNAseq_PBMCs and https://github.com/satheshsiva/PBMC_scRNAseq.
