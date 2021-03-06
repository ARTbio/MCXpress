Abstract
--------

Single-cell RNA-sequencing allows unbiased transcriptome profiling of hundreds to thousands of individual cells, enabling the analysis of genes expression at the cellular level. By uncovering cell heterogeneity in a given cell type or tissue, this technology leads to the potential discovery of rare subpopulations of cells which could be responsible of the onset and progression of specific diseases. Here we developed a fast and comprehensible hand on tool for the analysis of RNA-seq data using an approach based on Multiple Correspondence Analysis (MCA). MCA is a dimensionality reduction technique that allows the representation of both individuals (cells) and the variables (genes) within the same Euclidean space, thus allowing the simultaneous identification of subpopulations of cells and their gene signatures. MCA was coupled with pre-ranked gene set enrichment analysis; a powerful analytical method for characterizing differentially expressed genes/pathways within each groups of individuals. This singular combination allows a joint comparison of gene expression and pathway enrichment across all the groups, overcoming the limitations of a pair-wise based analysis used by standard differential expression approaches. An extensive visualization tools integrated in a shiny interface, allowing the user to intuitively inspect results and easily obtain a functional interpretation, complement the package. By integrating a singular dimensionality reduction approach, a robust analytic method and a sophisticated visualization tool in an easy to use framework, the MCXpress R package can enhance greatly the interpretation of single cell but also bulk RNA-sequencing data analysis.

Installation of MCXpress
------------------------

MCXpress is still under development but can be installed via github using the devtools library.

``` r
library(devtools)
install_github("Cortalak/MCXpress", ref = "master")
library(MCXpress)
```

Initialisation of MCXpress object
---------------------------------

In order to performs the analysis with MCXpress, the user needs to provide a gene expression matrix. The matrix (TPM, FPKM, count) must have as rownames the genes and as column names the cells/samples names. It should be preferrentially log transformed. first. The MCXpress object is first initialised using the `Initialise_MCXpress` function. It will be the main input for all the remaining MCXpress function and all the analysis results will be stored inside as list. Note that the genes name and cell/samples name must not contain any duplicates. We use here the example dataset of MCXpress derived from GSE64553 in GEO.

``` r
#GSE64553 an example matrix of log transformed single cell expression 
your_analysis <- Initialise_MCXpress(GSE64553)
#your_analysis$ExpressionMatrix to call your matrix
dim(your_analysis$ExpressionMatrix)
```

    ## [1] 36454    35

If not done previously, it is possible to filter the matrix by removing poorly expressed genes or selecting the most variable genes. By calling these function it will update the Expression Matrix inside the object.

``` r
your_analysis <- filter_outlier(your_analysis, percentage = 0.1, threshold = 3)
your_analysis <- select_most_variable_genes(X = your_analysis, ngenes = 10000)
#your_analysis$ExpressionMatrix to call your filtered matrix
dim(your_analysis$ExpressionMatrix)
```

    ## [1] 10000    35

Disjunctive Matrix Creation
---------------------------

In order to analyse the expression matrix with MCA, a transformation into a disjunctive matrix is required. The function `discretisation_01` performs a simple scaling of the expression value from 0 to 1. It is also possible to discretise the matrix using `discretisation_bsplines`.

``` r
your_analysis <- discretisation_01(your_analysis, scaled = FALSE)
#your_analysis <- Discretisation_Bsplines(your_analysis) For bsplines discretisation method
#your_analysis$Disjunctive_Matrix to call the Disjunctive Matrix
```

Multiple Corespondence Analysis
-------------------------------

Multiple Corespondece Analysis is performed on the created Disjunctive Matrix.

``` r
#Perform MCA you can choose the number of Axis to retain with the parameter Dim
your_analysis <- MCA(your_analysis)
```

``` r
#All attributes of MCA
your_analysis$MCA %>% attributes
```

    ## $names
    ## [1] "cells_principal"          "cells_standard"          
    ## [3] "genes_standard"           "genes_principal"         
    ## [5] "eigen_value"              "explained_eigen_variance"
    ## [7] "Methods"                  "Axis_Gene_Cor"           
    ## [9] "plot"                    
    ## 
    ## $class
    ## [1] "MCA"

``` r
#Plot first two Axes of MCA
your_analysis$MCA$plot
#Visualise your data interactively with your_analysis$Shiny
```

![](../dev/README_files/figure-markdown_github/MCA_2-1.svg) \#\# Clustering To distinguish the different subtypes of cell in the data clustering can be performed.

``` r
#Performing K-means Clustering
your_analysis <- cluster_kmeans(your_analysis, k = 6, dim=5)
```

    ## Performing Kmeans Clustering with  6  Cluster
    ##  Calculating Centroids 
    ## 
    ##  Calculating Distance between Cluster Centroids 
    ## 
    ##  Calculating Distance betwenn Cluster and Genes

``` r
#List of attributes for Clustering Section
your_analysis$cluster %>% attributes
```

    ## $names
    ## [1] "labels"                 "dim"                   
    ## [3] "nClusters"              "coord_centroids"       
    ## [5] "cluster_distances"      "gene_cluster_distances"
    ## [7] "closest_cluster"        "plot1"                 
    ## [9] "plot2"                 
    ## 
    ## $class
    ## [1] "Cluster"

``` r
#Finding Clustering Results
your_analysis$cluster$labels
```

    ## # A tibble: 35 x 2
    ##     Cluster               Sample
    ##       <chr>                <chr>
    ##  1 Cluster4   HFF_PD22_Rotenone1
    ##  2 Cluster4   HFF_PD22_Rotenone2
    ##  3 Cluster4   HFF_PD22_Rotenone3
    ##  4 Cluster4 HFF_PD22_woRotenone1
    ##  5 Cluster4 HFF_PD22_woRotenone2
    ##  6 Cluster4 HFF_PD22_woRotenone3
    ##  7 Cluster5   HFF_PD26_Rotenone1
    ##  8 Cluster5   HFF_PD26_Rotenone2
    ##  9 Cluster5   HFF_PD26_Rotenone3
    ## 10 Cluster1 HFF_PD26_woRotenone1
    ## # ... with 25 more rows

``` r
#Plot Cluster in GGplot2
your_analysis$cluster$plot1
```

![](../dev/README_files/figure-markdown_github/Clustering-1.svg) It is possible to visualise the most important genes for each cluster in the form of a heatmap. MCXpress will calculate the closest genes for each cluster.

``` r
Heatmap_Cluster(your_analysis, n = 3, plotly = F) 
```

    ## Joining, by = "Cells"

![](../dev/README_files/figure-markdown_github/Biomarkers-1.svg)

GSEA
----

Performs Geneset enrichment analysis using MCA distance metric. You need to supply a gmt file downloaded from BroadGSEA or your very own list of geneset.

``` r
your_analysis <- GSEA(your_analysis, GMTfile = Hallmark, nperm = 1000)
```

``` r
#Top 10 genes specific to cluster1
your_analysis$GSEA$Ranking$Cluster1 %>% head(10)
```

    ## RP11-184E9.2   AC090954.5       DIRAS2     COL6A4P1         ALX3 
    ##    1.0000000    0.5857108    0.3120977    0.2934250    0.2853957 
    ##       COL6A6       BMPR1B    LINC00599     ADAMTS19        BANK1 
    ##    0.2815028    0.2742676    0.2729280    0.2545550    0.2196198

``` r
#Top 10 genes specific to cluster2
your_analysis$GSEA$Ranking$Cluster2 %>% head(10)
```

    ##       FAM101A RP11-384C12.2   RP5-908D6.1  RP1-170D19.3         CXCL6 
    ##     1.0000000     0.9948866     0.9336523     0.9167362     0.9031346 
    ##         ASCL1  RP11-271M1.1       PRSS30P         KRT81      KRTAP3-1 
    ##     0.8787518     0.8628496     0.8433485     0.8321772     0.8086064

``` r
#Heatmap of GSEA results for each cluster
your_analysis %>% GSEA_Heatmap_Cluster(plotly=F, margin = c(15,20), pval = 0.25, es = 0.3)
```

![](../dev/README_files/figure-markdown_github/GSEA_Heatmap-1.svg)
