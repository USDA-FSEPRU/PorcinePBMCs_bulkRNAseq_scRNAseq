Cell type-enriched and cell type-specific gene identification
================
### Background
The normalized count table was used for differential gene expression (DGE) analysis with DESeq2 by setting the size factor for each sample to 1. A generalized linear model was fitted for each gene in the count table, with negative binomial response and log link function of the effect of cell types and pig subjects. nbinomWaldTest() was used to estimate and test the significance of regression coefficients with the following explicit parameter settings: betaPrior=FALSE,maxit=50000,useOptim=TRUE,useT=FALSE,useQR=TRUE. Cell type-enriched genes and cell type-specific genes were identified using the results function separately. A gene was labeled as cell-type enriched if the expression level (averaged across replicates) in one cell type was at least 2x higher than the average across all remaining cell types and adjusted p-value <0.05. A gene was labeled as cell type-specific if the averaged expression level in one cell type was at least 2x higher in pairwise comparison to the average in each other cell type and adjusted p-value <0.05 (Benjamini and Hochberg, 1995). Heatmaps to display specific genes were created as mentioned above.

### Install and load R packages
``` r
BiocManager::install("genefilter")
BiocManager::install("yarn")
BiocManager::install(c("DESeq2", "sva"))

library("DESeq2")
library("edgeR")
library("yarn")
library(pca3d)
library(sva)
library(limma)
library(rgl)
library("pheatmap")
library("RColorBrewer")
library(ggplot2)
library(dplyr)
library("genefilter")
library("AnnotationDbi")
library("WriteXLS")
library("gplots")
library("Biobase")
```

### STep 1: Combining count tables from paired-end and single end mapping on HPC
``` r
count <- read.delim("1.1.V97.Salmon.counts.17.samples.txt",
                    header = TRUE, as.is = TRUE)
load(file = "../Robjects/metadata.formatted.RData")
meta_fmt <- meta_fmt[-18, ]
meta_fmt <- meta_fmt[order(meta_fmt$CellType), ]
count <- count[, rownames(meta_fmt)]
all(rownames(meta_fmt) == colnames(count))

#### plot raw count per sample
pdf("1.0 Number of reads assigened to gene features.pdf",
    width = 6, height = 4)
par(mar=c(9.1,4.1,4.1,2.1))
barplot(colSums(count)/1000000, width = 0.8,
        xlim = c(0,24), ylim = c(0,50),
        ylab ="# Reads assigned to genes (M)",
        col = "gray",
        names.arg = meta_fmt$SampleID, las = 2)
dev.off()
```

### YARN
``` r
yarn_norm <- TRUE

if (yarn_norm)
{
    meta <- new("AnnotatedDataFrame", data = meta_fmt)
    exprset_imm <- new("ExpressionSet", exprs = as.matrix(count), phenoData = meta)

    exprset_imm_filtered <- filterLowGenes(exprset_imm, groups = "CellType",
                                           threshold = 4, minSamples = 1)

    exprset_imm_filtered <- normalizeTissueAware(obj = exprset_imm_filtered,
                                                 group = "CellType", normalizationMethod = "qsmooth")

    ##
    meta_fmt$CellType <- factor(meta_fmt$CellType,
                                levels = sort(unique(meta_fmt$CellType)))
    col <- rainbow(length(unique(meta_fmt$CellType)))[as.numeric(meta_fmt$CellType)]

    ## plot density distribution before and after "qsmooth" normalization"
    pdf(file = "1.1.Density plots showing gene expression before and after qsmooth normalization.pdf", height = 5, width = 12)
    op <- par(mfrow = c(1, 2))
    plotDensity(exprset_imm_filtered, groups = "CellType",
                normalized = FALSE,
                main = "Before Normalization",
                legendPos = "right", ylab = "Density",
                xlab = "log2(raw expression)")
    plotDensity(exprset_imm_filtered, groups = "CellType",
                normalized= TRUE, main="After Normalization",
                legendPos = "right", ylab = "Density",
                xlab = "Normalized expression, log2" )

    dev.off()
    par(op)

    pdf(file = "1.3. MDS plots showing sample distance.pdf", width = 8, height = 6)
    par(mar = c(5.1,4.1,4.1,7.1))
    col <- rainbow(length(unique(meta_fmt$CellType)))[as.numeric(meta_fmt$CellType)]

    plotCMDS(exprset_imm_filtered, pch = 16, cex = 2, normalized = TRUE,
             col = col)

    legend("bottomright",
           legend = sort(unique(meta_fmt$CellType)), pch = 16,
           col = rainbow(length(unique(meta_fmt$CellType)))[as.numeric(sort(unique(meta_fmt$CellType)))])

    dev.off()

    ### output normalized expression matrix
    normalized_expr_mat <- assayData(exprset_imm_filtered)[["normalizedMatrix"]]
    dim(normalized_expr_mat) ## 11171 genes
    count_back <- round(2^normalized_expr_mat - 1)
    counts <- count_back   ## 11270 genes
} else {
    #### filter count table to remove the lowly expressed genes
    cpms = cpm(count)
    keep = rowSums(cpms > 1) >= 2
    counts = count[keep,]  ## 15592 (mincount >= 4); 12805 (cpms > 1)
}

meta_fmt <- meta_fmt[colnames(counts), ]
## sanity check
all(rownames(meta_fmt) == colnames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta_fmt,
                              design = ~ 0+ CellType + PigID)
if (yarn_norm)
{
    "sizeFactors"(dds) <- rep(1, 17)
} else {
    dds <- estimateSizeFactors(dds)
}

dds <- estimateDispersions(dds)
### exploratory analysis
rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))

#### Heatmap showing sample distances
#### samples SP4 and SP8 are very different from SP5-7 of the same treatment
distancePlot <- function(sampleDists, sampleNames)
{
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- sampleNames
    colnames(sampleDistMatrix) <- sampleNames
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors)
    sampleDistMatrix
}

if (yarn_norm)
{
    suffix <- "-YARN.norm"
} else {
    suffix <- "-DEseq2.norm"
}

pdf(paste0("Figure 1.4 Heatmap showing original sample distances", suffix, ".pdf"),
    width = 6, height = 5)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                 sampleNames = rld$SampleID)
dev.off()

## PCA plot showing sample relationship between and within groups
pdf(paste0("Figure 1.5.PCA and MDS plots showing sample relationahip", suffix, ".pdf"),
    width=5, height=6)
plotPCA(rld, intgroup = c("CellType"))
dev.off()
```

### TPM data
``` r
TPM <- read.delim("1.0.V97.Salmon.abundance.TPM.17.samples.txt",
                  header = TRUE, as.is = TRUE)
TPM <- TPM[, rownames(meta_fmt)]
all(rownames(meta_fmt) == colnames(TPM))

## filtering TPM
TPM <- TPM[rownames(counts), ]

count <- TPM
## scaling TPM using YARN

yarn_norm <- TRUE

if (yarn_norm)
{
    meta <- new("AnnotatedDataFrame", data = meta_fmt)
    exprset_imm <- new("ExpressionSet", exprs = as.matrix(count + 1), phenoData = meta)

    exprset_imm_filtered <- filterLowGenes(exprset_imm, groups = "CellType",
                                           threshold = 1, minSamples = 1)

    exprset_imm_filtered <- normalizeTissueAware(obj = exprset_imm_filtered,
                                                 group = "CellType", normalizationMethod = "qsmooth")

    ##
    meta_fmt$CellType <- factor(meta_fmt$CellType,
                                levels = sort(unique(meta_fmt$CellType)))
    col <- rainbow(length(unique(meta_fmt$CellType)))[as.numeric(meta_fmt$CellType)]

    ## plot density distribution before and after "qsmooth" normalization"
    pdf(file = "1.1.Density plots showing gene expression in TPM scale before and after qsmooth normalization.pdf", height = 5, width = 12)
    op <- par(mfrow = c(1, 2))
    plotDensity(exprset_imm_filtered, groups = "CellType",
                normalized = FALSE,
                main = "Before Normalization",
                legendPos = "right", ylab = "Density",
                xlab = "log2(raw expression)")
    plotDensity(exprset_imm_filtered, groups = "CellType",
                normalized= TRUE, main="After Normalization",
                legendPos = "right", ylab = "Density",
                xlab = "Normalized expression, log2" )

    dev.off()
    par(op)

    pdf(file = "1.3.MDS plots showing sample distance in TPM scale.pdf", width = 8, height = 6)
    par(mar = c(5.1,4.1,4.1,7.1))
    col <- rainbow(length(unique(meta_fmt$CellType)))[as.numeric(meta_fmt$CellType)]

    plotCMDS(exprset_imm_filtered, pch = 16, cex = 2, normalized = TRUE,
             col = col)

    legend("bottomright",
           legend = sort(unique(meta_fmt$CellType)), pch = 16,
           col = rainbow(length(unique(meta_fmt$CellType)))[as.numeric(sort(unique(meta_fmt$CellType)))])

    dev.off()

    ### output normalized expression matrix
    normalized_expr_mat <- assayData(exprset_imm_filtered)[["normalizedMatrix"]]
    dim(normalized_expr_mat) ## 11171 genes
    count_back <- 2^normalized_expr_mat - 2  ## 11270 genes
} else {
    #### filter count table to remove the lowly expressed genes
    cpms = cpm(count)
    keep = rowSums(cpms > 1) >= 2
    counts = count[keep,]  ## 15592 (mincount >= 4); 12805 (cpms > 1)
}

meta_fmt <- meta_fmt[colnames(counts), ]
## sanity check
all(rownames(meta_fmt) == colnames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta_fmt,
                              design = ~ 0+ CellType + PigID)
if (yarn_norm)
{
    "sizeFactors"(dds) <- rep(1, 17)
} else {
    dds <- estimateSizeFactors(dds)
}

dds <- estimateDispersions(dds)
### exploratory analysis
rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
```

### Heatmap showing sample distances
#### samples SP4 and SP8 are very different from SP5-7 of the same treatment
``` r
distancePlot <- function(sampleDists, sampleNames)
{
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- sampleNames
    colnames(sampleDistMatrix) <- sampleNames
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors)
    sampleDistMatrix
}

if (yarn_norm)
{
    suffix <- "-YARN.norm"
} else {
    suffix <- "-DEseq2.norm"
}

pdf(paste0("Figure 1.4 Heatmap showing original sample distances in TPM scale", suffix, ".pdf"),
    width = 6, height = 5)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                 sampleNames = rld$SampleID)
dev.off()

## PCA plot showing sample relationship between and within groups
pdf(paste0("Figure 1.5.PCA and MDS plots showing sample relationahip in TPM scale", suffix, ".pdf"),
    width=5, height=6)
plotPCA(rld, intgroup = c("CellType"))
dev.off()

pdf(paste0("Figure 1.6.Boxplot showing normalized TPM", suffix, ".pdf"),
    width=6, height = 8)
count_back_0 <- count_back
par(mar=c(9.1,4.1,4.1,2.1))
colnames(count_back_0) <- meta_fmt$SampleID

boxplot(log2(count_back_0 + 1), las = 2, ylab = "log2(TPM + 1)")
dev.off()

write.table(count_back,
            paste0("Salmon.TPM.normalized.expression.table",
                   suffix, ".txt"),
            sep = "\t",
            quote = FALSE, row.names = TRUE)
```

### Preprocessing
``` r
normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))
write.table(normCounts, paste0("expression.count.table", suffix, ".txt"),
            sep = "\t",
            quote = FALSE, row.names = TRUE)


models <- model.matrix(~ 0 + CellType + PigID, data = meta_fmt)

contrast.matrix <- matrix(c(1, rep(-1/8, 8), 0,  ## CD21-
                            -1/8, 1, rep(-1/8, 7), 0,      ## CD21+
                            rep(-1/8, 2), 1, rep(-1/8, 6), 0,     ## CD8-CD4+
                            rep(-1/8, 3), 1, rep(-1/8, 5), 0,    ## CD8+CD4-
                            rep(-1/8, 4), 1, rep(-1/8, 4), 0,    ## CD8+CD4+
                            rep(-1/8, 5), 1, rep(-1/8, 3), 0,    ## GD-TCR
                            rep(-1/8, 6), 1, rep(-1/8, 2), 0,   ## Monocytes
                            rep(-1/8, 7), 1, -1/8, 0,  ## Neutrophils
                            rep(-1/8, 8), 1, 0    ## NK
                          ), nrow = 9, byrow = TRUE)
rownames(contrast.matrix) <- c("CD21-", "CD21+", "CD8-CD4+",
                                "CD8+CD4-", "CD8+CD4+", "GD-TCR",
                                "Monocytes", "Neutrophils",
                                "NK")



dds <- nbinomWaldTest(dds, modelMatrix = NULL,betaPrior = FALSE,
                      maxit = 50000, useOptim = TRUE, quiet = FALSE,
                      useT = FALSE, useQR = TRUE)
```


### Test models with 4 sva
``` r
pig_gene_ID_maps <- read.delim("pig_genes.Ensembl.92.txt", stringsAsFactors = FALSE)

output_DESeq <- function(i, dds, contrast.matrix, threshold)
{
    print(i)
    res <- results(dds, alpha = 0.01, contrast=contrast.matrix[i,],
                   lfcThreshold = log2(threshold), altHypothesis="greaterAbs")
    res<- res[order(res$padj),]
    res <- as.data.frame(res)
    res <- merge(res, pig_gene_ID_maps, by.x="row.names", by.y = "EnsemblID", all.x = T)
    colnames(res)[1] <- "EnsemblID"
    res
}

DESeq_out <- lapply(1:nrow(contrast.matrix), output_DESeq, dds = dds,
                    contrast.matrix = contrast.matrix, threshold = 2)
save.image(file = paste0("sorted.immune.cells", suffix, ".RData"))
file <- file.path("./", paste0("Sorted.immunecells.subtype-enriched.genes.w.pigID.effect", suffix, ".xlsx"))

WriteXLS(x = DESeq_out, ExcelFileName = file,
         row.names = FALSE, SheetNames = rownames(contrast.matrix))
```
### Specificity analysis
``` r
load("YARN normalization minCount=4 minSample=2/sorted.immune.cells-YARN.norm.RData")
models <- model.matrix(~ 0 + CellType+ PigID, data = meta_fmt)
contrast.matrix <- matrix(c(1, -1, rep(0, 8),  ## CD21-
                            1, 0, -1, rep(0, 7),
                            1, 0, 0, -1, rep(0, 6),
                            1, 0, 0, 0, -1, rep(0, 5),
                            1, 0, 0, 0, 0, -1, rep(0, 4),
                            1, 0, 0, 0, 0, 0, -1, rep(0, 3),
                            1, 0, 0, 0, 0, 0, 0, -1, rep(0, 2),
                            1, 0, 0, 0, 0, 0, 0, 0, -1, rep(0, 1),

                           -1, 1, rep(0, 8),
                            0, 1, -1, rep(0, 7),
                            0, 1, 0, -1, rep(0, 6),
                            0, 1, 0, 0, -1, rep(0, 5),
                            0, 1, 0, 0, 0, -1, rep(0, 4),
                            0, 1, 0, 0, 0, 0, -1, rep(0, 3),
                            0, 1, 0, 0, 0, 0, 0, -1, rep(0, 2),
                            0, 1, 0, 0, 0, 0, 0, 0, -1, rep(0, 1),

                           -1,  0, 1, rep(0, 7),
                            0, -1, 1, rep(0, 7),
                            0, 0,  1, -1, rep(0, 6),
                            0, 0,  1, 0, -1, rep(0, 5),
                            0, 0,  1, 0, 0, -1, rep(0, 4),
                            0, 0,  1, 0, 0, 0, -1, rep(0, 3),
                            0, 0,  1, 0, 0, 0, 0, -1, rep(0, 2),
                            0, 0,  1, 0, 0, 0, 0, 0, -1, rep(0, 1),

                            -1, 0, 0, 1, rep(0, 6),
                            0, -1, 0, 1, rep(0, 6),
                            0, 0, -1, 1, rep(0, 6),
                            0, 0,  0, 1, -1, rep(0, 5),
                            0, 0,  0, 1,  0, -1, rep(0, 4),
                            0, 0,  0, 1,  0, 0, -1, rep(0, 3),
                            0, 0,  0, 1,  0, 0, 0, -1, rep(0, 2),
                            0, 0,  0, 1,  0, 0, 0, 0, -1, rep(0, 1),

                            -1, 0, 0, 0, 1, rep(0, 5),
                            0, -1, 0, 0, 1, rep(0, 5),
                            0, 0, -1, 0, 1, rep(0, 5),
                            0, 0,  0,-1, 1, rep(0, 5),
                            0, 0,  0, 0, 1, -1, rep(0, 4),
                            0, 0,  0, 0, 1, 0, -1, rep(0, 3),
                            0, 0,  0, 0, 1, 0, 0, -1, rep(0, 2),
                            0, 0,  0, 0, 1, 0, 0, 0, -1, rep(0, 1),

                           -1, 0, 0, 0, 0, 1, rep(0, 4),
                           0, -1, 0, 0, 0, 1, rep(0, 4),
                           0, 0, -1, 0, 0, 1, rep(0, 4),
                           0, 0,  0,-1, 0, 1, rep(0, 4),
                           0, 0,  0, 0,-1, 1, rep(0, 4),
                           0, 0,  0, 0, 0, 1, -1, rep(0, 3),
                           0, 0,  0, 0, 0, 1, 0, -1, rep(0, 2),
                           0, 0,  0, 0, 0, 1, 0, 0, -1, rep(0, 1),

                           -1, 0, 0, 0, 0, 0, 1, rep(0, 3),
                           0, -1, 0, 0, 0, 0, 1, rep(0, 3),
                           0, 0, -1, 0, 0, 0, 1, rep(0, 3),
                           0, 0,  0,-1, 0, 0, 1, rep(0, 3),
                           0, 0,  0, 0,-1, 0, 1, rep(0, 3),
                           0, 0,  0, 0, 0,-1, 1, rep(0, 3),
                           0, 0,  0, 0, 0, 0, 1, -1, rep(0, 2),
                           0, 0,  0, 0, 0, 0, 1, 0, -1, rep(0, 1),


                           -1, 0, 0, 0, 0, 0, 0, 1, rep(0, 2),
                           0, -1, 0, 0, 0, 0, 0, 1, rep(0, 2),
                           0, 0, -1, 0, 0, 0, 0, 1, rep(0, 2),
                           0, 0,  0,-1, 0, 0, 0, 1, rep(0, 2),
                           0, 0,  0, 0,-1, 0, 0, 1, rep(0, 2),
                           0, 0,  0, 0, 0,-1, 0, 1, rep(0, 2),
                           0, 0,  0, 0, 0, 0,-1, 1, rep(0, 2),
                           0, 0,  0, 0, 0, 0, 0, 1, -1, rep(0, 1),

                           -1, 0, 0, 0, 0, 0, 0, 0, 1, rep(0, 1),
                           0, -1, 0, 0, 0, 0, 0, 0, 1, rep(0, 1),
                           0, 0, -1, 0, 0, 0, 0, 0, 1, rep(0, 1),
                           0, 0,  0,-1, 0, 0, 0, 0, 1, rep(0, 1),
                           0, 0,  0, 0,-1, 0, 0, 0, 1, rep(0, 1),
                           0, 0,  0, 0, 0,-1, 0, 0, 1, rep(0, 1),
                           0, 0,  0, 0, 0, 0,-1, 0, 1, rep(0, 1),
                           0, 0,  0, 0, 0, 0, 0,-1, 1, rep(0, 1)

), nrow = 72, byrow = TRUE)

cell_types <- c("CD21-", "CD21+", "CD8-CD4+",
                "CD8+CD4-", "CD8+CD4+", "GD-TCR",
                "Monocytes", "Neutrophils","NK")
contrast_names <- expand.grid(cell_types, cell_types, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
contrast_names <- contrast_names[ contrast_names[, 1] != contrast_names[, 2], ]

rownames(contrast.matrix) <- paste(contrast_names[, 2],
                                    contrast_names[, 1], sep = "-")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta_fmt,
                              design = ~ 0+ CellType + PigID)
if (yarn_norm)
{
    "sizeFactors"(dds) <- rep(1, 18)
} else {
    dds <- estimateSizeFactors(dds)
}

dds <- estimateDispersions(dds)

dds <- nbinomWaldTest(dds, modelMatrix = NULL,betaPrior = FALSE,
                      maxit = 50000, useOptim = TRUE, quiet = FALSE,
                      useT = FALSE, useQR = TRUE)
```

### Test models with 4 sva
``` r
pig_gene_ID_maps <- read.delim("pig_genes.Ensembl.92.txt", stringsAsFactors = FALSE)

output_DESeq <- function(i, dds, contrast.matrix, threshold)
{
    print(i)
    res <- results(dds, alpha = 0.01, contrast=contrast.matrix[i,],
                   lfcThreshold = log2(threshold), altHypothesis="greaterAbs")
    res<- res[order(res$padj),]
    res <- as.data.frame(res)
    res <- res[sort(rownames(res)), ]
    colnames(res) <- paste(colnames(res), rownames(contrast.matrix)[i], sep = "_")
    res
}


for (j in 1:8)
{
    DESeq_out <- lapply((7*j - 6):(7*j), output_DESeq, dds = dds,
                        contrast.matrix = contrast.matrix,
                        threshold = 2)

    cmbined <- do.call(cbind, DESeq_out)
    cmbined <- cmbined[, -c(seq(7, 48, by = 6))]
    cmbined <- merge(cmbined, pig_gene_ID_maps, by.x = "row.names",
                     by.y = "EnsemblID", all.x = TRUE)
    colnames(cmbined)[1] <- "EnsemblID"

    file <- file.path("./", paste0("Sorted.immunecells.subtype-enriched.genes.w.pigID.effect.", gsub("([^-]+).+", "\\1", rownames(contrast.matrix)[7*j], perl = TRUE), ".txt"))

    write.table(x = cmbined, file = file, sep = "\t", quote = FALSE, row.names = FALSE)

}
```

#### Find cell-type specific genes
``` r
files <- dir(".", "*effects")
file.names <- gsub(".+?effects.(.+?).txt", "\\1", files)

data_lsit <- lapply(files, function(.x){
    temp <- read.delim(.x, header = TRUE, as.is = TRUE, check.names = FALSE)
    temp <- temp[rowSums(temp[, seq(7, 42, by = 5)] <= 0.05) == 7 &
                      rowSums(temp[, seq(3, 42, by = 5)] >= 1) == 7, ]
    print(dim(temp))
    temp
})

file <- file.path("./", paste0("Sorted.immunecells.subtype-specific.genes.w.pigID.effect-YARN.mincount=4.minsample=2.xlsx"))

WriteXLS(x = data_lsit, ExcelFileName = file,
         row.names = FALSE, SheetNames = file.names)
```

### Heatmap of union of DEGs
``` r
load("DESeq2 normalization cpm.gt.1 minSample=2/sorted.immune.cells-DEseq2.norm.RData")
common_degs <- unique(do.call(c,lapply(DESeq_out, function(.x){
    ids <- .x$EnsemblID[.x$padj <= 0.05 & .x$log2FoldChange >= 2 ]
    ids <- ids[!is.na(ids)]
})))  ## YARN (2611)  DESeq2 (3173)


common_pig_degs <- pig_gene_ID_maps[pig_gene_ID_maps$EnsemblID %in% common_degs, c(1,2)]
rownames(common_pig_degs) <- common_pig_degs[, 1]
common_pig_degs <- with(common_pig_degs, common_pig_degs[common_degs, ])
common_pig_degs$Symbol <- ifelse(common_pig_degs$Symbol == "", common_pig_degs$EnsemblID, common_pig_degs$Symbol)
common_pig_degs$Symbol[duplicated(common_pig_degs$Symbol)] #GZMA

## avoid duplicates
common_pig_degs$Symbol <- ifelse(common_pig_degs$Symbol == "GZMA", common_pig_degs$EnsemblID, common_pig_degs$Symbol)

if (yarn_norm)
{
    adj_exp <- read.delim("YARN normalization minCount=4 minSample=2/expression.count.table-YARN.norm.txt")
} else {
    adj_exp <- read.delim("DESeq2 normalization cpm.gt.1 minSample=2/expression.count.table-DEseq2.norm.txt")
}

adj_exp <- voom(adj_exp, design = model.matrix(~ CellType + PigID, data = meta_fmt), plot = FALSE)

adj_exp <- adj_exp$E

all(colnames(adj_exp) == rownames(meta_fmt))
adj_exp_deg <- adj_exp[rownames(adj_exp) %in% common_degs, ]
adj_exp_deg <- adj_exp_deg[common_degs, ]

all(rownames(common_pig_degs) == rownames(adj_exp_deg))
rownames(adj_exp_deg) <- common_pig_degs$Symbol


annotation_col <- data.frame(CellType = meta_fmt$CellType)
rownames(annotation_col) <- rownames(meta_fmt)
annotation_colors <- list(CellType = c("CD21N"  = "#FF0000FF",
                                     "CD21P"    = "#FFAA00FF",
                                     "CD8NCD4P" = "#AAFF00FF",
                                     "CD8PCD4N" = "#00FF00FF",
                                     "CD8PCD4P" = "#00FFAAFF",
                                     "GDT"      = "#00AAFFFF",
                                     "Monocytes" = "#0000FFFF",
                                     "Neutrophils" = "#AA00FFFF",
                                     "NK"          = "#FF00AAFF"))


pdf("1.7 Heatmap showing cell type enriched genes.pdf",
    height = 55, width = 6)
heat <- pheatmap(as.matrix(adj_exp_deg),
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 clustering_method = "ward.D2",
                 clustering_distance_rows = "correlation",
                 clustering_distance_cols = "euclidean",
                 scale = "row",
                 #cellheight = 0.3,
                 border_color= NA,
                 annotation_col = annotation_col,
                 annotation_colors = annotation_colors,
                 color =  rev(redgreen(75)),
                 show_colnames = FALSE,
                 show_rownames = TRUE,
                 fontsize = 5,
                 fontsize_row = 1.2,
                 fontsize_col = 5)
print(heat)
dev.off()
```

### References
Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol 57(1), 289-300. doi: 10.2307/2346101.
