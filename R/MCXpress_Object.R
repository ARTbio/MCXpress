print.Dim_Red_Object <- function(obj, ...) {
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
      "Pearson Correlation between Individuals Coordinate on Axis and Gene Expression",
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
      "$Graph2",
      "$Shiny"
    )
  DESCRIPTION <-
    c(
      "",
      "Number of Cluster",
      "A vector indicating the cluster for each samples",
      "Table indicating for each genes the cluster that is expressing the most",
      "Coordinates of the Cluster Centroids",
      "Cluster Plot in Cell Space",
      "Centroids of cluster plotted in Genespace",
      "Interactive Plot"
    )
  tibble(NAME, DESCRIPTION) %>%  print.data.frame(row.names = F, right = F)
}

print.FA_Object <- function(obj, ...) {
  cat('Functionnal Analysis Results', "\n", "\n")
  NAME <- c("", "$Ranking", "$GSEA_Results", "$Shiny")
  DESCRIPTION <-
    c(
      "",
      "Table with Gene Ranking for each cluster",
      "fgsea package Gene Set Enrichment Analysis Results for each cluster",
      "Interactive Plot"
    )
  tibble(NAME, DESCRIPTION) %>%  print.data.frame(row.names = F, right = F)
}
