
#' Initialisation of the MCXpress Object
#'
#' Create a MCXpressObject from an Expression Matrix
#'
#' @param X DESCRIPTION.
#' @param min_reads DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' @export
Initialise_MCXpress <- function(X, min_reads = NULL) {
  if (X %>% is.matrix) {
    if ((X %>% rownames %>% is.null) | (X %>% colnames %>%
      is.null)) {
      errormessage <- "Gene name should be the rownames of the matrix and Sample name the column name"
      stop(errormessage)
    } else {
      MCXpress <- list()
      X <- as.matrix((X[apply(X, 1, var) > 0,]))
      X <- X[str_length(X %>% rownames) > 0,]
      X <- X %>% subset(X %>% rownames %>% duplicated %>%  not)
      colnames(X)<-gsub(pattern = " ", replacement ="_", x = colnames(X))
      colnames(X)<-gsub(pattern = "-", replacement ="_", x = colnames(X))
      colnames(X)<-gsub(pattern = "\\.", replacement ="_", x = colnames(X))
      if (min_reads %>% is.numeric()) {
        X <- X[rowSums(X) > min_reads, ]
      }
      MCXpress$ExpressionMatrix <- X
      class(MCXpress) <- "MCXpress_object"
      return(MCXpress)
    }
  } else {
    errormessage <- "The input is not a matrix"
    stop(errormessage)
  }
}

#' Interactive plot enrichment score of a pathway
#'
#' Plotly version of the plotEnrichment function of the fgsea package
#'
#' @param pathway A single geneset
#' @param stats Ranking of the genes
#' @param gseaParam GSEA parameter value
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
plotlyEnrichment <- function(pathway, stats, gseaParam = 0) {
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  name <- rnk[(rnk %>% names) %in% pathway] %>% sort %>%  names
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
    returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  plot_ly(data = toPlot) %>% add_lines(x = ~x, y = ~y, mode = "lines",
    line = list(color = "rgb(154, 240, 24)", width = 2),
    name = "Enrichment") %>% add_lines(x = ~x, y = ~min(bottoms),
    line = list(color = "rgb(110, 193, 248)", width = 2,
      dash = "dash"), name = "Lower limit", hoverinfo = "text",
    text = ~round(min(bottoms), digits = 4)) %>% add_lines(x = ~x,
    y = ~max(tops), line = list(color = "rgb(250, 150, 10)",
      width = 2, dash = "dash"), name = "Upper limit",
    hoverinfo = "text", text = ~round(max(tops), digits = 4)) %>% add_segments(x =~pathway, xend=~pathway ,y =~diff/2 ,yend =~-diff/2, line=list(color = "rgb(0, 10, 10)",width=1), text=~paste0(name,
      "</br>", "Rank:", pathway), showlegend = FALSE, hoverinfo="text")
}


