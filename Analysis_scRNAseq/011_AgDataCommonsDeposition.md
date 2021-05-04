Creating Data Files For Ag Data Commons Deposition
================

### Load required software packages

Refer to sessionInfo() at the bottom of the page for the R and package
versions used.

``` r
library(Seurat)
library(SeuratDisk)
library(DropletUtils)
library(dplyr)
```

### Import relevant data

Read in our Seurat object from .h5Seurat file format:

``` r
pbmc <- LoadH5Seurat('/home/Jayne.Wiarda/PBMCscRNAseq/PBMC7_AllCells.h5Seurat')
```

### Extract relevant metadata & count matrix:

Create CellRanger-like output files of our Seurat object, including cell
barcode IDs (barcodes.tsv.gz), gene/feature IDs (features.tsv.gz), and
raw count matrix\* (matrix.mtx.gz):

  - The ‘raw’ count matrix is actually gene counts following ambient RNA
    removal with the R package, SoupX. During ambient RNA removal, we
    specified to calculate non-integer count estimations, so most gene
    counts are actually non-integer values in this matrix but should
    still be treated as raw/unnormalized data that requires further
    normalization/transformation.

<!-- end list -->

``` r
write10xCounts(x = pbmc@assays$RNA@counts, path = "/home/Jayne.Wiarda/PBMCscRNAseq/PBMC7_AllCells", version = "3") 
```

Extract metadata:

``` r
meta <- pbmc@meta.data
colnames(meta)
```

    ##  [1] "nCount_SCT"              "nFeature_SCT"           
    ##  [3] "orig.ident"              "nCount_RNA"             
    ##  [5] "nFeature_RNA"            "barcode"                
    ##  [7] "SampleID"                "Loupe"                  
    ##  [9] "UmiSums"                 "GenesDetected"          
    ## [11] "prcntTop"                "prcntMito"              
    ## [13] "DuplicatedBarcodes"      "PassViability"          
    ## [15] "PassGenesDet"            "PassLibSize"            
    ## [17] "PassBarcodeFreq"         "PassAll"                
    ## [19] "Scrublet"                "PassScrub"              
    ## [21] "integrated_snn_res.1.85" "seurat_clusters"        
    ## [23] "PaperIDs"                "phyloorder"             
    ## [25] "neworder"                "celltypes"

``` r
meta <- meta %>% select(nCount_RNA, nFeature_RNA, Loupe, prcntMito, 
                        Scrublet, seurat_clusters, PaperIDs, celltypes) # keep only the most pertinent categories
write.csv(meta, '/home/Jayne.Wiarda/PBMCscRNAseq/PBMC7_AllCells_meta.csv')
```

Extract UMAP, t-SNE, & PCA coordinates:

``` r
pca = Embeddings(pbmc[["pca"]])
pca <- as.data.frame(pca)
pca$CellID <- rownames(pca)
write.csv(pca, '/home/Jayne.Wiarda/PBMCscRNAseq/PBMC7_AllCells_PCAcoord.csv')

tsne = Embeddings(pbmc[["tsne"]])
tsne <- as.data.frame(tsne)
tsne$CellID <- rownames(tsne)
write.csv(tsne, '/home/Jayne.Wiarda/PBMCscRNAseq/PBMC7_AllCells_tSNEcoord.csv')

umap = Embeddings(pbmc[["umap"]])
umap <- as.data.frame(umap)
umap$CellID <- rownames(umap)
write.csv(umap, '/home/Jayne.Wiarda/PBMCscRNAseq/PBMC7_AllCells_UMAPcoord.csv')
```

I will also export the UMAP & t-SNE coordinates for data subsets created
for only CD4+ ab T cells (clusters 0, 3, 4, 28) and only gd T cells
(clusters 6, 21, 24, 31). I will not export these separately using
write10xCounts or the full metadata since that can be obtained from the
output files of the PBMC7\_AllCells dataset.

For CD4 T cells only:

``` r
CD4only <- LoadH5Seurat('/home/Jayne.Wiarda/PBMCscRNAseq/PBMC7_CD4only.h5Seurat')

tsne = Embeddings(CD4only[["tsne"]])
tsne <- as.data.frame(tsne)
tsne$CellID <- rownames(tsne)
write.csv(tsne, '/home/Jayne.Wiarda/PBMCscRNAseq/PBMC7_CD4only_tSNEcoord.csv')

umap = Embeddings(CD4only[["umap"]])
umap <- as.data.frame(umap)
umap$CellID <- rownames(umap)
write.csv(umap, '/home/Jayne.Wiarda/PBMCscRNAseq/PBMC7_CD4only_UMAPcoord.csv')
```

For gd T cells only:

``` r
GDonly <- LoadH5Seurat('/home/Jayne.Wiarda/PBMCscRNAseq/PBMC7_GDonly.h5Seurat')

tsne = Embeddings(GDonly[["tsne"]])
tsne <- as.data.frame(tsne)
tsne$CellID <- rownames(tsne)
write.csv(tsne, '/home/Jayne.Wiarda/PBMCscRNAseq/PBMC7_GDonly_tSNEcoord.csv')

umap = Embeddings(GDonly[["umap"]])
umap <- as.data.frame(umap)
umap$CellID <- rownames(umap)
write.csv(umap, '/home/Jayne.Wiarda/PBMCscRNAseq/PBMC7_GDonly_UMAPcoord.csv')
```

### View session information

``` r
sessionInfo()
```

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Red Hat Enterprise Linux
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.3.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] dplyr_1.0.2                 DropletUtils_1.8.0         
    ##  [3] SingleCellExperiment_1.10.1 SummarizedExperiment_1.18.2
    ##  [5] DelayedArray_0.14.1         matrixStats_0.57.0         
    ##  [7] Biobase_2.48.0              GenomicRanges_1.40.0       
    ##  [9] GenomeInfoDb_1.24.2         IRanges_2.22.2             
    ## [11] S4Vectors_0.26.1            BiocGenerics_0.34.0        
    ## [13] SeuratDisk_0.0.0.9013       Seurat_3.2.2               
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.15             colorspace_2.0-0       deldir_0.2-3          
    ##   [4] ellipsis_0.3.1         ggridges_0.5.2         XVector_0.28.0        
    ##   [7] spatstat.data_1.5-2    leiden_0.3.6           listenv_0.8.0         
    ##  [10] ggrepel_0.9.1          bit64_4.0.5            fansi_0.4.1           
    ##  [13] R.methodsS3_1.8.1      codetools_0.2-16       splines_4.0.2         
    ##  [16] knitr_1.30             polyclip_1.10-0        jsonlite_1.7.2        
    ##  [19] ica_1.0-2              cluster_2.1.0          R.oo_1.24.0           
    ##  [22] png_0.1-7              uwot_0.1.9             HDF5Array_1.16.1      
    ##  [25] shiny_1.5.0            sctransform_0.3.1      compiler_4.0.2        
    ##  [28] httr_1.4.2             dqrng_0.2.1            assertthat_0.2.1      
    ##  [31] Matrix_1.2-18          fastmap_1.0.1          lazyeval_0.2.2        
    ##  [34] limma_3.44.3           cli_2.2.0              later_1.1.0.1         
    ##  [37] htmltools_0.5.0        tools_4.0.2            rsvd_1.0.3            
    ##  [40] igraph_1.2.6           GenomeInfoDbData_1.2.3 gtable_0.3.0          
    ##  [43] glue_1.4.2             RANN_2.6.1             reshape2_1.4.4        
    ##  [46] Rcpp_1.0.5             spatstat_1.64-1        vctrs_0.3.5           
    ##  [49] nlme_3.1-148           lmtest_0.9-38          xfun_0.19             
    ##  [52] stringr_1.4.0          globals_0.14.0         mime_0.9              
    ##  [55] miniUI_0.1.1.1         lifecycle_0.2.0        irlba_2.3.3           
    ##  [58] goftest_1.2-2          future_1.21.0          edgeR_3.30.3          
    ##  [61] zlibbioc_1.34.0        MASS_7.3-51.6          zoo_1.8-8             
    ##  [64] scales_1.1.1           promises_1.1.1         spatstat.utils_1.17-0 
    ##  [67] rhdf5_2.32.4           RColorBrewer_1.1-2     yaml_2.2.1            
    ##  [70] reticulate_1.18        pbapply_1.4-3          gridExtra_2.3         
    ##  [73] ggplot2_3.3.2          rpart_4.1-15           stringi_1.5.3         
    ##  [76] BiocParallel_1.22.0    bitops_1.0-6           rlang_0.4.9           
    ##  [79] pkgconfig_2.0.3        evaluate_0.14          lattice_0.20-41       
    ##  [82] Rhdf5lib_1.10.1        ROCR_1.0-11            purrr_0.3.4           
    ##  [85] tensor_1.5             patchwork_1.1.0        htmlwidgets_1.5.3     
    ##  [88] cowplot_1.1.0          bit_4.0.4              tidyselect_1.1.0      
    ##  [91] parallelly_1.21.0      RcppAnnoy_0.0.17       plyr_1.8.6            
    ##  [94] magrittr_2.0.1         R6_2.5.0               generics_0.1.0        
    ##  [97] pillar_1.4.7           withr_2.3.0            mgcv_1.8-31           
    ## [100] fitdistrplus_1.1-3     RCurl_1.98-1.2         survival_3.2-7        
    ## [103] abind_1.4-5            tibble_3.0.4           future.apply_1.6.0    
    ## [106] crayon_1.3.4           hdf5r_1.3.3            KernSmooth_2.23-17    
    ## [109] plotly_4.9.2.1         rmarkdown_2.7          locfit_1.5-9.4        
    ## [112] grid_4.0.2             data.table_1.13.4      digest_0.6.27         
    ## [115] xtable_1.8-4           tidyr_1.1.2            httpuv_1.5.4          
    ## [118] R.utils_2.10.1         munsell_0.5.0          viridisLite_0.3.0
