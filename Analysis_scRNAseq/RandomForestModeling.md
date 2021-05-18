Random Forest Modeling
================
### Background
	Using the Iowa State University High Performance Computing Nova Cluster, the GNU gcc compiler and R modules were loaded into the environment. The R packages caret (Kuhn et al.) and ranger (Wright et al.) were used to create the random forest models. The ranger package is a more efficient C++ implementation of the randomForest package (Brieman et al.), resulting in faster run times for larger datasets.

	A normalized count matrix was used as the input data for the random forest models. Each cell was labeled by its previously defined cluster. The Seurat subset function was used to subset the data based on which clusters were included in each model. Three different types of models were created. The first were pairwise models where the training data included only cells from two different clusters (ex. Clusters 0 & 3). The second were cell type specific models where the training data included only cells from clusters that were all a certain cell type (ex. CD4 T cells – Clusters 0, 3, 4, & 28). The third type of model compared cells from all clusters from a certain cell type (ex. CD4 T cells – Clusters 0, 3, 4, & 28) to the cells from all other clusters.

	Each model was ran as classification for the cluster identity of each cell. The number of trees created was 500, the target node size was 1, the number of variables was 14,386, and the number of variables to sample at each split (Mtry) was 119. Each tree in the model is grown from a bootstrap resampling process which provides an out-of-bag (OOB) error that provides an efficient and reasonable approximation of the test error. Variable importance was assigned by measuring node impurity (Impurity) and using permutations (Permutation) to calculate variable importance. Features that reduce error in predictive accuracy are ranked as more important. A high error rate in the model would suggest that cells from the groups being compared are more similar to each other, whereas a low error rate in the model would suggest the cells from each cluster are unique. Variable importance was used to find genes or sets of genes that can be used to identify certain types of cells or discriminate groups of cells from one another.

### Load required software packages

``` r
module load gcc
module load r

R

library(caret)
library(ranger)
```

### Load normalized count data

Read in our normalized count matrix with Cluster IDs output from Seurat object

In this example the matrix is 28810 x 14387, with 28810 cells in rows and 14386 genes in columns. The first column is the Cluster ID for each cell.

``` r
data <- read.csv(“data.csv”)
data <- as.matrix(data)
dim(data)
[1] 28810 14387
```

### Subset data based on which clusters to include in training data (ex. all CD4 T cell clusters)

Subset to only CD4 T cell clusters:

In this example the matrix is 5082 x 14387, with 5082 cells in rows and 14386 genes in columns. The first column is the Cluster ID for each cell.

``` r
subset_data <- subset(data, rownames(data) == “0” | rownames(data) == “3” | rownames(data) == “4” | rownames(data) == “28”)
dim(subset_data)
[1] 5082 14387
```

### Train the random forest model on cluster ID

Run model using permutation as the importance variable

``` r
set.seed(100)
rf <- ranger(cluster ~., num.trees = 500, importance = “permutation”, data = subset_data)
rf
Ranger result

Call:
 ranger(cluster ~ ., num.trees = 500, importance = "permutation", data = subset_data)

Type:                             Classification
Number of trees:                  500
Sample size:                      5082
Number of independent variables:  14386
Mtry:                             119
Target node size:                 1
Variable importance mode:         permutation
Splitrule:                        gini
OOB prediction error:             19.05 %

importance <- rf$variable.importance
sorted <- sort(importance, decreasing=TRUE)
write.csv(sorted, "RF_Permutation.csv")
```

Run model using impurity as the importance variable

``` r
set.seed(100)
rf <- ranger(cluster ~., num.trees = 500, importance = “impurity”, data = subset_data)
rf
Ranger result

Call:
 ranger(cluster ~ ., num.trees = 500, importance = "impurity", data = subset_data)

Type:                             Classification
Number of trees:                  500
Sample size:                      5082
Number of independent variables:  14386
Mtry:                             119
Target node size:                 1
Variable importance mode:         impurity
Splitrule:                        gini
OOB prediction error:             19.05 %
importance <- rf$variable.importance
sorted <- sort(importance, decreasing=TRUE)
write.csv(sorted, "RF_Impurity.csv")
```

### View session information

``` r
sessionInfo()
R version 4.0.4 (2021-02-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server 7.9 (Maipo)

Matrix products: default
BLAS/LAPACK: /opt/rit/spack-app/linux-rhel7-x86_64/gcc-10.2.0/openblas-0.3.13-yr6yp5gfuo7zduwfvqi47ypiywy5aikc/lib/libopenblas-r0.3.13.so

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
[1] ranger_0.12.1   caret_6.0-86    ggplot2_3.3.3   lattice_0.20-41

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.6           pillar_1.5.1         compiler_4.0.4
 [4] gower_0.2.2          plyr_1.8.6           tools_4.0.4
 [7] class_7.3-18         iterators_1.0.13     rpart_4.1-15
[10] ipred_0.9-11         lubridate_1.7.10     lifecycle_1.0.0
[13] tibble_3.1.0         gtable_0.3.0         nlme_3.1-152
[16] pkgconfig_2.0.3      rlang_0.4.10         Matrix_1.3-2
[19] foreach_1.5.1        DBI_1.1.1            prodlim_2019.11.13
[22] stringr_1.4.0        withr_2.4.1          dplyr_1.0.5
[25] pROC_1.17.0.1        generics_0.1.0       vctrs_0.3.6
[28] recipes_0.1.15       stats4_4.0.4         nnet_7.3-15
[31] grid_4.0.4           tidyselect_1.1.0     data.table_1.14.0
[34] glue_1.4.2           R6_2.5.0             fansi_0.4.2
[37] survival_3.2-10      lava_1.6.9           reshape2_1.4.4
[40] purrr_0.3.4          magrittr_2.0.1       ModelMetrics_1.2.2.2
[43] splines_4.0.4        MASS_7.3-53.1        scales_1.1.1
[46] codetools_0.2-18     ellipsis_0.3.1       assertthat_0.2.1
[49] timeDate_3043.102    colorspace_2.0-0     utf8_1.2.1
[52] stringi_1.5.3        munsell_0.5.0        crayon_1.4.1
```

### References
Kuhn, Max. "Building predictive models in R using the caret package." J Stat Softw 28.5 (2008): 1-26.

Wright, Marvin N., and Andreas Ziegler. "ranger: A fast implementation of random forests for high dimensional data in C++ and R." arXiv preprint arXiv:1508.04409 (2015).

Breiman, Leo. "Random forests." Machine learning 45.1 (2001): 5-32.
