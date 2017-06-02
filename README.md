Abstract
--------

Single-cell RNA-sequencing allows unbiased transcriptome profiling of
hundreds to thousands of individual cells, enabling the analysis of
genes expression at the cellular level. By uncovering cell heterogeneity
in a given cell type or tissue, this technology leads to the potential
discovery of rare subpopulations of cells which could be responsible of
the onset and progression of specific diseases. Here we developed a fast
and comprehensible hand on tool for the analysis of RNA-seq data using
an approach based on Multiple Correspondence Analysis (MCA). MCA is a
dimensionality reduction technique that allows the representation of
both individuals (cells) and the variables (genes) within the same
Euclidean space, thus allowing the simultaneous identification of
subpopulations of cells and their gene signatures. MCA was coupled with
pre-ranked gene set enrichment analysis; a powerful analytical method
for characterizing differentially expressed genes/pathways within each
groups of individuals. This singular combination allows a joint
comparison of gene expression and pathway enrichment across all the
groups, overcoming the limitations of a pair-wise based analysis used by
standard differential expression approaches. An extensive visualization
tools integrated in a shiny interface, allowing the user to intuitively
inspect results and easily obtain a functional interpretation,
complement the package. By integrating a singular dimensionality
reduction approach, a robust analytic method and a sophisticated
visualization tool in an easy to use framework, the MCXpress R package
can enhance greatly the interpretation of single cell but also bulk
RNA-sequencing data analysis.

Installation of MCXpress
------------------------

MCXpress is still under development but can be installed via github
using the devtools library.

    library(devtools)
    install_github("cbl-imagine/MCXpress", ref = "master")
    library(MCXpress)

Initialisation of MCXpress object
---------------------------------

In order to performs the analysis with MCXpress, the user needs to
provide a gene expression matrix. The matrix (TPM, FPKM, count) must
have as rownames the genes and as column names the cells/samples names.
It should be preferrentially log transformed. first. The MCXpress object
is first initialised using the `Initialise_MCXpress` function. It will
be the main input for all the remaining MCXpress function and all the
analysis results will be stored inside as list. Note that the genes name
and cell/samples name must not contain any duplicates. We use here the
example dataset of MCXpress derived from GSE64553 in GEO.

    #GSE64553 an example matrix of log transformed single cell expression 
    your_analysis <- Initialise_MCXpress(GSE64553)
    #your_analysis$ExpressionMatrix to call your matrix
    dim(your_analysis$ExpressionMatrix)

    ## [1] 36454    35

If not done previously, it is possible to filter the matrix by removing
poorly expressed genes or selecting the most variable genes. By calling
these function it will update the Expression Matrix inside the object.

    your_analysis <- filter_outlier(your_analysis, percentage = 0.1, threshold = 3)
    your_analysis <- select_most_variable_genes(X = your_analysis, ngenes = 10000)
    #your_analysis$ExpressionMatrix to call your filtered matrix
    dim(your_analysis$ExpressionMatrix)

    ## [1] 10000    35

Disjunctive Matrix Creation
---------------------------

In order to analyse the expression matrix with MCA, a transformation
into a disjunctive matrix is required. The function `discretisation_01`
performs a simple scaling of the expression value from 0 to 1. It is
also possible to discretise the matrix using `discretisation_bsplines`.

    your_analysis <- discretisation_01(your_analysis, scaled = FALSE)
    #your_analysis <- Discretisation_Bsplines(your_analysis) For bsplines discretisation method
    #your_analysis$Disjunctive_Matrix to call the Disjunctive Matrix

Multiple Corespondence Analysis
-------------------------------

Multiple Corespondece Analysis is performed on the created Disjunctive
Matrix.

    #Perform MCA you can choose the number of Axis to retain with the parameter Dim
    your_analysis <- MCA(your_analysis, Dim = 3)

    ## Peforming MCA...
    ## ===========================================================================
    ## Calculating Cell to Cell Distance...
    ## MCA is finished

    #All attributes of MCA
    your_analysis$MCA %>% attributes

    ## $names
    ##  [1] "cells_principal"          "cells_standard"          
    ##  [3] "genes_standard"           "genes_principal"         
    ##  [5] "Cell2Cell_Distance"       "eigen_value"             
    ##  [7] "explained_eigen_variance" "Methods"                 
    ##  [9] "Axis_Gene_Cor"            "plot"                    
    ## 
    ## $class
    ## [1] "MCA_Object"

    #Plot first two Axes of MCA
    your_analysis$MCA$plot

![](C:\Users\Akira\AppData\Local\Temp\RtmpcxtMOu\preview-283059411c56.dir\MCXpress_vignette_files/figure-markdown_strict/unnamed-chunk-6-1.png)

    #Visualise your data interactively with your_analysis$Shiny

Clustering
----------

To distinguish the different subtypes of cell in the data clustering is
performed.

    #Performing K-means Clustering
    your_analysis <- cluster_kmeans(your_analysis, k = 6)

    ## Performing Kmeans Clustering with  6  Cluster
    ##  Calculating Centroids 
    ## 
    ##  Calculating Distance between Cluster Centroids 
    ## 
    ##  Calculating Distance betwenn Cluster and Genes

    #List of attributes for Clustering Section
    your_analysis$cluster %>% attributes

    ## $names
    ## [1] "labels"                 "nClusters"             
    ## [3] "coord_centroids"        "cluster_distances"     
    ## [5] "gene_cluster_distances" "closest_cluster"       
    ## [7] "plot1"                  "plot2"                 
    ## 
    ## $class
    ## [1] "Cluster_Object"

    #Finding Clustering Results
    your_analysis$cluster$labels

    ## # A tibble: 35 x 2
    ##                  Sample  Cluster
    ##                   <chr>    <chr>
    ##  1   HFF_PD22_Rotenone1 Cluster2
    ##  2   HFF_PD22_Rotenone2 Cluster2
    ##  3   HFF_PD22_Rotenone3 Cluster2
    ##  4 HFF_PD22_woRotenone1 Cluster2
    ##  5 HFF_PD22_woRotenone2 Cluster2
    ##  6 HFF_PD22_woRotenone3 Cluster2
    ##  7   HFF_PD26_Rotenone1 Cluster6
    ##  8   HFF_PD26_Rotenone2 Cluster6
    ##  9   HFF_PD26_Rotenone3 Cluster6
    ## 10 HFF_PD26_woRotenone1 Cluster4
    ## # ... with 25 more rows

    #Plot Cluster in GGplot2
    your_analysis$cluster$plot1

![](C:\Users\Akira\AppData\Local\Temp\RtmpcxtMOu\preview-283059411c56.dir\MCXpress_vignette_files/figure-markdown_strict/unnamed-chunk-7-1.png)

GSEA
----

Performs Geneset enrichment analysis using MCA distance metric.

    your_analysis <- GSEA(your_analysis, GMTfile = reactome_gmtfile, nperm = 1000)

    ## Calculating ranking of genes correlation for each axis 
    ## Beginning enrichment analysis for axis 
    ## processing: Axis1 
    ## processing: Axis2 
    ## 
    ## Calculating ranking of genes for each clusters 
    ## Beginning enrichment analysis for clusters
    ## processing: Cluster1 
    ## processing: Cluster2 
    ## processing: Cluster3 
    ## processing: Cluster4 
    ## processing: Cluster5 
    ## processing: Cluster6 
    ## processing: Origin 
    ## Enrichment Analysis Completed

    #List of attributes for GSEA Section
    your_analysis$GSEA %>% attributes

    ## $names
    ## [1] "GSEA_Results_Axis" "RankingAxis"       "GSEA_Results"     
    ## [4] "Ranking"           "GMTfile"           "Pathways"         
    ## [7] "AllRanking"        "gseaParam"        
    ## 
    ## $class
    ## [1] "GSEA_Object"

    #Top 10 genes specific to cluster1
    your_analysis$GSEA$Ranking$Cluster1 %>% tail(10)

    ##         RUNX3         MEOX2       SLC24A3       FAM101A       PRSS30P 
    ##          9991          9992          9993          9994          9995 
    ## RP11-224O19.2         ASCL1       SPATA22      KRTAP3-1 CTD-2231H16.1 
    ##          9996          9997          9998          9999         10000

    #Top 10 genes specific to cluster2
    your_analysis$GSEA$Ranking$Cluster2 %>% tail(10)

    ##  ASS1P11   GPRC5C   RNF175   SLC7A2    SFRP2     MLC1    AJAP1     NEFL 
    ##     9991     9992     9993     9994     9995     9996     9997     9998 
    ## TMEM132D    SMOC2 
    ##     9999    10000
