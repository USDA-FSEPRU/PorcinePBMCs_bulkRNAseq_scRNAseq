Reference Based Label Transfer
================
### Background

	We used a publically available human reference CITE-seq dataset of 162,000 PBMCs measured with 228 antibodies to transfer cell type annotations onto our newly generated pig scRNA-seq PBMC dataset (Hao et al.). Due to the cross-species comparison, we distilled the human reference and pig query datasets to only include 1:1 orthologous genes. Orthologous genes were identified by querying Ensemble BioMart for pig and human orthologous genes. Subsequently pig gene IDs were updated to match the orthologous human gene ID. To map the pig dataset to the human reference dataset, we employed similar methods that were used to create the human reference dataset. The Iowa State University High Performance Computing Nova Cluster was used to load the hdf5 and R modules into the programming environment. Within R, the Seurat, SeuratDisk, ggplot2, and patchwork libraries were loaded. The human reference dataset containing only orthologous genes was partitioned by sample and separately normalized using SCTransform. SelectIntegrationFeatures was ran on the list of partitioned samples to identify 3000 features. PrepSCTIntegration was ran on the list of partitioned samples before features were scaled and used for PCA. Integration anchors were found using Day 0 samples as a reference. The identified anchors were then used to integrate the partitioned human samples. The integrated human dataset was used to create a supervised PCA (sPCA). The sPCA was created by first finding variable features in the integrated dataset for PCA. FindNeighbors was ran using the first 20 dimensions of the PCA. Finally, the sPCA was created and the first 50 dimensions were calculated to reduce the time needed for multiple object mapping.
	The pig dataset containing only orthologous genes was partitioned by sample and separately normalized using SCTransform. Anchors were found between the human reference dataset and each individual pig sample. These anchors were used to calculate mapping scores for each cell. These mapping scores provided a 0-1 confidence value of how well an individual cell, or an average of all cells in a cluster, were represented by the human reference dataset.
	A feature of the human reference CITE-seq dataset was that cell types had been well annotated using weighted RNA and protein modalities. After mapping our pig dataset to the human reference dataset, we were able to transfer cell type annotations from the human reference dataset to each of the cells in the pig dataset. We predicted three increasingly specific levels of cell type annotations from the human reference dataset, along with corresponding prediction scores. Prediction scores provided a 0-1 confidence value for an individual cell type prediction, based on how many nearby human cells shared the same cell type annotation that was predicted. Finally, we took the predicted Level 2 cell annotations and projected back onto our original UMAP of the pig data.
	In order to identify cells from the pig dataset that were not well represented by the human reference dataset, we integrated the two datasets and performed a de novo visualization. In the de novo visualization, cell types that are unique or more highly represented in the pig dataset should remain separated and appear pink in the figure.

### Load required software packages

``` r
module load hdf5/1.10.5-openmpi3-fx7c4wu
module load r

R

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
```

### Set size of globals
``` r
options(future.globals.maxSize = 10000 * 4000^2)
```

### Load Suerat reference object created with only orthologous genes, partition and perform SCTransform individually

``` r
reference <- readRDS("Human_Reference_Only_Orthos.rds”)

reference.list <- SplitObject(reference, split.by = "SampleID")

for (i in 1:length(reference.list)) {
reference.list[[i]] <- SCTransform(reference.list[[i]], verbose = FALSE)
}
```

### Follow initial steps of SCT integration

``` r
reference.features <- SelectIntegrationFeatures(object.list = reference.list, nfeatures = 3000)

reference.list <- PrepSCTIntegration(object.list = reference.list, anchor.features = reference.features, verbose = FALSE)
```

### Follow reciprocal PCA workflow to integrate data, use “Day 0” samples as reference

``` r
reference.list <- lapply(X=reference.list, FUN=function(x) {
+ x <- ScaleData(x, features = reference.features, verbose = FALSE)
+ x <- RunPCA(x, features = reference.features, verbose = FALSE)})
anchors <- FindIntegrationAnchors(object = reference.list, reference = c(8,9,11,12,14,16,18,23), reduction = "rpca", dims = 1:50)

reference.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

saveRDS(reference.integrated, "Reference_Integrated.rds")
```
### Create spca

``` r
reference.integrated@active.assay <- “SCT”

reference.integrated <- FindVariableFeatures(reference.integrated, selection.method = "vst", nfeatures = 3000)

reference.integrated <- RunPCA(reference.integrated, verbose = FALSE)
```

### Check dimensions (1:20 dims chosen)
``` r
ElbowPlot(reference.integrated)

reference.integrated <- FindNeighbors(reference.integrated, dims = 1:20)
reference.integrated <- RunUMAP(reference.integrated, dims = 1:20)
reference.integrated <- RunSPCA(reference.integrated, assay = "SCT", graph = "integrated_snn")
```

### Calculate first 50 dims to save time with multiple object mapping
``` r
reference.integrated <- FindNeighbors(
    object = reference.integrated,
    reduction = "spca",
    dims = 1:50,
    graph.name = "spca.annoy.neighbors",
    k.param = 50,
    cache.index = TRUE,
    return.neighbor = TRUE,
)

saveRDS(reference.integrated, "Reference_Integrated_SPCA.rds")
SaveAnnoyIndex(object = reference.integrated[[“spca.annoy.neighbors"]], file = "reftmp.idx")
```

### If you want to load the annoy index (not saved in the rds object)

``` r
reference.integrated[[“spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = reference.integrated[[“spca.annoy.neighbors"]], file = "reftmp.idx")
```

### Load Suerat object created with only orthologous genes, partition and perform SCTransform individually

``` r
query <- readRDS("Pig_Query_Only_Orthos.rds")

query.list <- SplitObject(query, split.by = "Sample")

for (i in 1:length(query.list)) {
query.list[[i]] <- SCTransform(query.list[[i]], verbose = FALSE)
}
```

### Find anchors between each pig sample and the human reference dataset. Use the anchors found to calculate a mapping score for each cell. Then transfer the three levels of cell type annotation from the human reference dataset to each cell.
``` r
for (i in 1:length(query.list)) {
	query1 <- query.list[[i]]

	anchors <- FindTransferAnchors(reference = reference.integrated, query = query1, 	normalization.method="SCT", reference.reduction="spca", dims=1:50)
	#1 Found 5770 anchors
	#2 Found 4942 anchors
	#3 Found 5421 anchors
	#4 Found 4775 anchors
	#5 Found 6514 anchors
	#6 Found 8726 anchors
	#7 Found 9971 anchors

	query1 <- AddMetaData(
    		object = query1,
    		metadata = MappingScore(anchors = anchors),
    		col.name = "mapping.score"
	)

	query1 <- MapQuery(
   		 anchorset = anchors,
    		query = query1,
    		reference = reference,
    		refdata = list(
    			celltype.l1 = “celltype.l1",
    			celltype.l2 = "celltype.l2”,
    			celltype.l3 = “celltype.l3”,
    		),
    		reference.reduction = "spca",
	)

	query.list[[i]] <- query1
}
```

### Merge the human and pig datasets to create a de novo visualization
``` r
reference.integrated$id <- 'reference'
query.list$id <- 'query'
refquery <- merge(reference, query.list)
refquery[["spca"]] <- merge(reference[["spca"]], query.list[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)
```

### View session information

``` r
sessionInfo()
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server 7.9 (Maipo)

Matrix products: default
BLAS:   /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/r-3.6.3-sxv6dw3qtl4ih3sncnebmvbgmdngyqwz/rlib/R/lib/libRblas.so
LAPACK: /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/r-3.6.3-sxv6dw3qtl4ih3sncnebmvbgmdngyqwz/rlib/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] patchwork_1.1.1       ggplot2_3.3.3         SeuratDisk_0.0.0.9013
[4] Seurat_3.9.9.9010

loaded via a namespace (and not attached):
  [1] nlme_3.1-151           matrixStats_0.57.0     bit64_4.0.5
  [4] RcppAnnoy_0.0.18       RColorBrewer_1.1-2     httr_1.4.2
  [7] sctransform_0.3.2.9001 tools_3.6.3            R6_2.5.0
 [10] irlba_2.3.3            rpart_4.1-15           KernSmooth_2.23-18
 [13] uwot_0.1.10            lazyeval_0.2.2         mgcv_1.8-33
 [16] colorspace_2.0-0       withr_2.3.0            tidyselect_1.1.0
 [19] gridExtra_2.3          bit_4.0.4              compiler_3.6.3
 [22] cli_2.2.0              hdf5r_1.3.3            plotly_4.9.2.1
 [25] scales_1.1.1           lmtest_0.9-38          spatstat.data_1.7-0
 [28] ggridges_0.5.2         pbapply_1.4-3          spatstat_1.64-1
 [31] goftest_1.2-2          stringr_1.4.0          digest_0.6.27
 [34] spatstat.utils_1.17-0  pkgconfig_2.0.3        htmltools_0.5.0
 [37] parallelly_1.23.0      fastmap_1.0.1          htmlwidgets_1.5.3
 [40] rlang_0.4.10           shiny_1.5.0            generics_0.1.0
 [43] zoo_1.8-8              jsonlite_1.7.2         ica_1.0-2
 [46] dplyr_1.0.2            magrittr_2.0.1         Matrix_1.3-0
 [49] fansi_0.4.1            Rcpp_1.0.5             munsell_0.5.0
 [52] abind_1.4-5            reticulate_1.18        lifecycle_0.2.0
 [55] stringi_1.5.3          MASS_7.3-53            Rtsne_0.15
 [58] plyr_1.8.6             grid_3.6.3             parallel_3.6.3
 [61] listenv_0.8.0          promises_1.1.1         ggrepel_0.9.0
 [64] crayon_1.3.4           miniUI_0.1.1.1         deldir_0.2-3
 [67] lattice_0.20-41        cowplot_1.1.0          splines_3.6.3
 [70] tensor_1.5             pillar_1.4.7           igraph_1.2.6
 [73] future.apply_1.7.0     reshape2_1.4.4         codetools_0.2-18
 [76] leiden_0.3.6           glue_1.4.2             data.table_1.13.4
 [79] png_0.1-7              vctrs_0.3.6            httpuv_1.5.4
 [82] gtable_0.3.0           RANN_2.6.1             purrr_0.3.4
 [85] polyclip_1.10-0        tidyr_1.1.2            assertthat_0.2.1
 [88] future_1.21.0          rsvd_1.0.3             mime_0.9
 [91] xtable_1.8-4           later_1.1.0.1          survival_3.2-7
 [94] viridisLite_0.3.0      tibble_3.0.4           cluster_2.1.0
 [97] globals_0.14.0         fitdistrplus_1.1-3     ellipsis_0.3.1
[100] ROCR_1.0-11
```

### References
Hao, Y., Hao, S., Andersen-Nissen, E., Mauck, W. M., Zheng, S., Butler, A., ... & Satija, R. (2020). Integrated analysis of multimodal single-cell data. bioRxiv.
