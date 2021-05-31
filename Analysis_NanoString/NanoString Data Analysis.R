library("NanoStringDiff")
library("metaMA")
library("calibrate")
library("sva")
library("limma")
library("lme4")
library("optimx")
library("lmerTest")
library("dplyr")
library("DESeq2")
library("MASS")
library("pheatmap")
library("gplots")
library("RColorBrewer")
library("grDevices")
library("ggplot2")

## using sva to reduce hidden variations
get.sva <- function (expr.data=NULL, meta.data=NULL)
{
    mod <-  model.matrix(~ 0 + Cell.Population + RIN, data=meta.data)
    mod0 <- model.matrix(~ RIN, data=meta.data)
    
    num.sva <-svaseq(as.matrix(expr.data), mod, mod0)$n.sv
    sv <- svaseq(as.matrix(expr.data), mod, mod0,n.sv=num.sva)$sv
    
    colnames(sv)<- paste0("sv",1:num.sva)
    
    meta.data.sva <-cbind(meta.data,sv ) 
    
    meta.data.sva
}

## adjust expression for hidden variations for EDA plots
## log2cpm
get.adj <- function(expr.data=NULL, design=NULL,  meta.data=NULL)
{
    
    v <- voom(expr.data, design=design)
    fit <- lmFit(v, design=design)
    adj.exp <- v$E -fit$coefficients[, 9:ncol(design)] %*% 
        t(design[ ,9:ncol(design)])
    
    list(raw = v$E, adj = adj.exp)
}


#### Heatmap showing sample distances
#### samples SP4 and SP8 are very different from SP5-7 of the same treatment
distancePlot <- function(sampleDists, sampleNames)
{
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- sampleNames
    colnames(sampleDistMatrix) <- sampleNames
    colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors)
    sampleDistMatrix 
}

## Heatmap
plot.NanoString.heatmap <- function(norm_expr, metadata, filename)
{
    row_cor <- cor(t(as.matrix(norm_expr)), method = "spearman")
    row_distance <- as.dist(1 - row_cor)
    col_distance <- dist(t(norm_expr))
    
    ## heatmap annotation settings
    annotation_col <- data.frame(Population = 
                                     factor(gsub("\\s+$","", 
                                                 metadata$Cell.Population),
                                levels =c("B cell", "CD4- CD8+", "CD4+ CD8-", 
                                          "CD4+ CD8+", "Gd T cell", "Monocyte",
                                          "NK cell", "PBMC"), 
                                labels = c("B cell", "CD4-CD8+", "CD4+CD8-", 
                                           "CD4+CD8+", "Gd T cell", "Monocyte", 
                                           "NK cell", "PBMC")))
    
    rownames(annotation_col) <- colnames(norm_expr)
    
    colors <- c("red", "orange", "yellow","brown", "blue", 
                "cyan", "green", "purple")
    names(colors) <- c("B cell", "CD4-CD8+", "CD4+CD8-",
                       "CD4+CD8+", "Gd T cell", "Monocyte", 
                       "NK cell", "PBMC")
    annotation_colors <- list(Population = colors)
    
    ## plot heatmap
    pdf(filename, width = 5, height = 11.5)
    pheatmap(as.matrix(norm_expr), 
             cluster_rows = T,
             cluster_cols = T,
             cellwidth = 5, cellheight = 3.3,
             clustering_method="ward.D2",
             clustering_distance_rows = row_distance,
             clustering_distance_cols = col_distance,
             scale="row",
             border_color= NA,
             annotation_col = annotation_col,
             annotation_colors = annotation_colors,
             treeheight_row = 30,
             color = greenred(75),
             show_colnames= FALSE,
             show_rownames= TRUE,
             fontsize = 8,
             fontsize_row= 3, 
             fontsize_col = 3)
    dev.off()
}

setwd("C:\\Users\\Haibo\\Desktop\\Tuggle_lab-help\\Nanostring\\NanoString")
expr_data <- read.delim("Crystal Naostring raw data 21518 KW.csv", 
                        sep = ",", header = TRUE, as.is = TRUE, 
                        check.names = FALSE)
colnames(expr_data) <- gsub("\\s+(_\\d+)", "\\1", colnames(expr_data))
metadata <- read.delim("Crystal Nanostring metadata.txt",
                       header = TRUE, as.is =TRUE)
metadata <- within(metadata, {
    PigID <- factor(PigID)
    Cell.Population <- factor(gsub("\\s+$","", Cell.Population))
})
metadata <- with(metadata, metadata[order(Cell.Population, PigID), ])

expr_data <- expr_data[, c(colnames(expr_data)[1:3], metadata$Nanostring.ID)]
colnames(expr_data) <- c(colnames(expr_data)[1:3], metadata$SampleID)
write.table(expr_data, file = "Crystal Naostring formatted raw data.csv",
            sep = ",", row.names = FALSE, quote = FALSE)

## Normalization of NanoString raw data using NanoStringDiff
path = file.path("C:\\Users\\Haibo\\Desktop\\Tuggle_lab-help\\Nanostring\\NanoString", "Crystal Naostring formatted raw data.csv")
NanoStringData <- createNanoStringSetFromCsv(path = path, 
                                             header=TRUE, 
                                             designs = metadata)
NanoStringData <- estNormalizationFactors(NanoStringData)
## get normalized data for further analysis
norm_expr <- NanoStringDataNormalization(path = path, 
                                         header=TRUE, 
                                         designs = metadata)[[2]]

## remove ISG20 which is not expressed in any sample
norm_expr <- norm_expr[rownames(norm_expr) != "ISG20", ]
norm_expr <- log(norm_expr+1)

plot.NanoString.heatmap(norm_expr = norm_expr, 
                        metadata = metadata, 
                        filename = 
             paste("Heatmap showing unadjusted 230 gene expression",
                   "in diff blood cell populations.pdf"))

## Surrogate variable analysis because of large hidden variations asw shown by
## the abobe heatmap
meta.sva <- get.sva(expr.data= norm_expr, meta.data=metadata)
model <- model.matrix(~ 0 + Cell.Population + RIN + sv1 + sv2 +sv3, 
                      data = meta.sva)
## adjust for sva
adj.exp <- get.adj(expr.data=norm_expr, design=model, 
                   meta.data=meta.sva)

## sample distance-based clustering
pdf("Distance between samples with or without adjusting hidden variations.pdf",
    height =7, width =7.5)
sampleDists_adj <- dist(t(adj.exp$adj))
distancePlot(sampleDists_adj, gsub("B cell ", "B cell",
                                   gsub("20180209_USDA \\d+ \\d+_\\d+ |.RCC","",
                                        meta.sva$Nanostring.ID)))
sampleDists_raw <- dist(t(adj.exp$raw))
distancePlot(sampleDists_raw, gsub("B cell ", "B cell",
                                   gsub("20180209_USDA \\d+ \\d+_\\d+ |.RCC","",
                                        meta.sva$Nanostring.ID)))
dev.off()

## heatmap showing gene expression before SVA-adjusting
plot.NanoString.heatmap(norm_expr = adj.exp$raw, 
                        metadata, 
                        filename = "Raw Nanostring gene exp-scaled heatmap.pdf")

## heatmap showing gene expression adjusted for 3 surrogate variables
plot.NanoString.heatmap(norm_expr = adj.exp$adj, 
                        metadata, 
            filename = "SVA adjusted Nanostring gene exp-scaled heatmap.pdf")

## using voom adjusted data to fit linear mixed effect model
adj_exp <- as.data.frame(adj.exp$raw)
adj_exp <- as.data.frame(t(adj_exp))
adj_exp <- cbind(adj_exp, meta.sva)

contrasts <- matrix(c(0, -1, rep(0, 10), ## B cell vs. others
                      0, 0, -1, rep(0, 9),
                      0, 0, 0, -1, rep(0, 8),
                      0, 0, 0, 0, -1, rep(0, 7),
                      0, 0, 0, 0, 0, -1, rep(0, 6),
                      0, 0, 0, 0, 0, 0, -1, rep(0, 5),
                      0, 1, -1, rep(0, 9),  ## CD4-CD8+
                      0, 1, 0, -1, rep(0, 8),
                      0, 1, 0, 0, -1, rep(0, 7),
                      0, 1, 0, 0, 0, -1, rep(0, 6),
                      0, 1, 0, 0, 0, 0, -1, rep(0, 5),
                      0, 0, 1, -1, rep(0, 8),  ## CD4+CD8-
                      0, 0, 1, 0, -1, rep(0, 7),
                      0, 0, 1, 0, 0, -1, rep(0, 6),
                      0, 0, 1, 0, 0, 0, -1, rep(0, 5),
                      0, 0, 0, 1, -1, rep(0, 7),  ## CD4+CD8+
                      0, 0, 0, 1, 0, -1, rep(0, 6),
                      0, 0, 0, 1, 0, 0, -1, rep(0, 5),
                      0, 0, 0, 0, 1, -1, rep(0, 6),  ## GD T cells
                      0, 0, 0, 0, 1, 0, -1, rep(0, 5),
                      0, 0, 0, 0, 0, 1, -1, rep(0, 5)), ## Monocytes
                    ncol = 12, byrow = TRUE)

rownames(contrasts) <- c("B - CD4-CD8+", "B - CD4+CD8-", "B - CD4+CD8+",
                         "B - gd.T", "B - Mo", "B - NK", "CD4-CD8+ - CD4+CD8-",
                         "CD4-CD8+ - CD4+CD8+", "CD4-CD8+ - gd.T", 
                         "CD4-CD8+ - Mo", "CD4-CD8+ - NK", 
                         "CD4+CD8- - CD4+CD8+", "CD4+CD8- - gd.T", 
                         "CD4+CD8- - Mo", "CD4+CD8- - NK", "CD4+CD8+ - gd.T", 
                         "CD4+CD8- - Mo", "CD4+CD8- - NK", "gd.T - Mo", 
                         "gd.T - NK", "Mo - NK" )

test_results <- vector("list", 229)
test_lsmeans <- vector("list", 229)
coeff <- vector("list", "229")
for (i in 1:(ncol(adj_exp )-9))
{
    
    try({
        lmm <- lmer(adj_exp[, i] ~ 1 + Cell.Population + RIN + sv1 + 
                        sv2 + sv3 + (1|PigID),
                    data = adj_exp, REML = TRUE,
                    control = lmerControl(optimizer ='optimx', 
                                          optCtrl=list(maxit=20000, 
                                                       method = "nlminb")),
                    start = NULL)
        
        coeff[[i]] <- coefficients(summary(lmm))[,1, drop = FALSE]
        lmm <- as_lmerModLmerTest(lmm, tol = 1e-08)
        lmmtest <- do.call("rbind", lapply(1: nrow(contrasts), function(.x) {
                      contest1D(lmm, L = contrasts[.x, ],
                      rhs = 0, confint = TRUE,
                      level = 0.95, ddf = "Kenward-Roger")}))
        test_results[[i]] <- cbind(data.frame(contrast = rownames(contrasts),
                                    Gene = rep(colnames(adj_exp)[i], 21)), 
                                   lmmtest)

        ls_means <- as.data.frame(ls_means(lmm, which = NULL, level = 0.95,
                 ddf = c("Kenward-Roger"), pairwise = FALSE))
        ls_means <-  cbind(data.frame(Gene = rep(colnames(adj_exp)[i], 8)),
                           ls_means)
        test_lsmeans[[i]] <- ls_means[, -c(2)]
    })
}

coeffs <- do.call('cbind', coeff)
adjusted_express <- as.matrix(adj_exp[, 1:(ncol(adj_exp )-9)]) - 
    as.matrix(meta.sva[, c(5,7:9)]) %*% coeffs[9:12, ]
head(adjusted_express)
rownames(adjusted_express) <- gsub("\\s+(_\\d+)", "\\1", 
                                   meta.sva$Nanostring.ID)

final_contrast_test <- do.call('rbind', test_results)
final_lsmeans <- do.call('rbind', test_lsmeans)
colnames(final_lsmeans)[2] <- "Cell.population"

### plot lsmean +- SEM
## TPM is the same CPM
pdf("2-24-2019 Genes.lsmeans.expression.pdf", width = 12, height = 15)
for (i in seq(1, length(unique(final_lsmeans$Gene)), 60))
{
    temp <- final_lsmeans[final_lsmeans$Gene %in% 
                              unique(final_lsmeans$Gene)[i:(i+59)], ]
    p <- ggplot(temp, aes(x= Cell.population, y= Estimate, 
                          color = Cell.population)) +
        geom_point()+ geom_errorbar(aes(ymin= Estimate - `Std. Error`,
                                        ymax= Estimate + `Std. Error`),
                                    size=0.3, width=.2) + 
        facet_wrap(~ Gene, nrow = 10, ncol = 6, scales ="free_y") +
        ylab("log2(TPM)") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
        labs(color = "Cell population") + xlab(NULL)
    print(p)
}
dev.off()
