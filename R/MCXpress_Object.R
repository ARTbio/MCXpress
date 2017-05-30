print.MCXpress_object <- function(obj, ...){
   cat("\n$ExpressionMatrix\tyour input expression matrix")
if(obj$Disjunctive_Matrix %>% is.null() %>% not()){
  cat("\n$DisjunctiveMatrix\tgenerated Burt Matrix")
  if(obj$MCA %>% is.null() %>% not()){
    cat("\n$MCA\t\tMCA results")
      if(obj$cluster %>% is.null() %>% not()){
        cat("\n$cluster\t\tClustering Results")
          if(obj$GSEA %>% is.null() %>% not()){
            cat("\n$GSEA\tGSEA Results")
  }
  }
  }
}
}


print.MCA_Object <- function(obj, ...) {
  cat(obj$Methods, 'Dimension Reduction Results', "\n", "\n")
  NAME <-
    c(
      "",
      "$Cells_Standard",
      "$Cells_Principal",
      "$Genes_Standard",
      "$Genes_Principal",
      "$Cell2Cell_Distance",
      "$Axis_Gene_Correlation",
      "$Graph",
      "$Eigen_Value",
      "$Wilcoxon",
      "$Shiny"
    )
  DESCRIPTION <-
    c(
      "",
      "Cells raw coordinate",
      "Cells normalised coordinate",
      "Genes raw coordinate",
      "Genes normalised coordinate",
      "Distance Between Individuals in the Euclidean Space",
      "Pearson Correlation between cells coordinate and Gene Expression",
      "Plot of Cell space",
      "eigenvalue",
      "Wilcoxon Test p-value for the number of axis to keep",
      "Interactive Visualisation of the MCA Results"
    )
  tibble(NAME, DESCRIPTION) %>%  print.data.frame(row.names = F, right = F)
}

print.Cluster_Object <- function(obj, ...) {
  cat('Clustering Results', "\n", "\n")
  NAME <-
    c(
      "",
      "$nClusters",
      "$Cluster_Quali",
      "$Closest_Cluster",
      "$Coord_Centroids",
      "$Graph1",
      "$Graph2"
    )
  DESCRIPTION <-
    c(
      "",
      "Number of Cluster",
      "A vector indicating the cluster for each samples",
      "Table indicating for each genes the closest cluster centroids",
      "Coordinates of the Cluster Centroids",
      "Cluster Plot in Cell Space",
      "Centroids of cluster plotted in Genespace"
    )
  tibble(NAME, DESCRIPTION) %>%  print.data.frame(row.names = F, right = F)
}

print.GSEA_Object <- function(obj, ...) {
  cat('Gene Set Enrichment Analysis Results', "\n", "\n")
  NAME <- c("", "$Ranking", "$GSEA_Results","$Ranking_Axis", "$GSEA_Results_Axis","$GMTfile", "$gseaParam")
  DESCRIPTION <-
    c(
      "",
      "Gene Ranking for each cluster",
      "fgsea package Gene Set Enrichment Analysis Results for each cluster",
      "Gene Ranking for each axis",
      "fgsea package Gene Set Enrichment Analysis Results for each axis",
      "GMTfile used for GSEA",
      "GSEA Parameter"
    )
  tibble(NAME, DESCRIPTION) %>%  print.data.frame(row.names = F, right = F)
}

