Initialise_MCXpress <- function(X, min_reads = NULL){
  if (X %>% is.matrix){
    if((X %>% rownames %>% is.null)|(X %>%  colnames %>% is.null)){errormessage<-"Gene name should be the rownames of the matrix and Sample name the column name"
    stop(errormessage)} else{
      MCXpress <- list()
      X <- as.matrix(unique(X[apply(X, 1, var) > 0, ]))
      X <- X[str_length(X %>% rownames) > 0, ]
      if(min_reads %>%  is.numeric()){
      X <- X[rowSums(X) > min_reads, ]}
      MCXpress$ExpressionMatrix<-X
      class(MCXpress)<- "MCXpress_object"
      return(MCXpress)}}else{errormessage<-"The input is not a matrix"
      stop(errormessage)}
}

plotlyEnrichment<-function(pathway,stats, gseaParam=1)
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  name<-statsAdj[as.vector(na.omit(match(pathway, names(statsAdj))))] %>% sort %>% names %>% rev
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

  data_markers <- data.frame(Genes=rep(name, each=100),x=rep(pathway, 100) %>% sort, y=rep(seq(from=diff/2, to=-diff/2, length.out = 100),pathway %>% length))
  plot_ly(data = toPlot) %>%  add_lines( x=~x, y=~y, mode="lines", line = list(color = 'rgb(154, 240, 24)', width = 2), name= "Enrichment") %>%
    add_lines(x=~x ,y = ~min(bottoms), line = list(color = 'rgb(110, 193, 248)', width = 2, dash = 'dash'), name= "Lower limit", hoverinfo="text", text=~round(min(bottoms), digits = 4)) %>%
    add_lines(x=~x ,y = ~max(tops), line = list(color = 'rgb(250, 150, 10)', width = 2, dash = 'dash'), name="Upper limit",  hoverinfo="text", text=~round(max(tops), digits = 4)) %>%
    add_markers(data=data_markers , x=~x ,y =~y, marker=list(color = 'rgb(0, 10, 10)', size=1), name="Genes", hoverinfo="text", text=~paste0(Genes,'</br>',"Rank:", x), showlegend=FALSE)
}

Reactome_GridPlot <-
  function(X, Info_File, Hierarchy_File, GMTfile, p = 0.25) {
    GMT <- GMTfile %>%  attributes() %>%  as_tibble
    colnames(GMT) <- 'pathway'
    Category <- Reactome_Category_Generator(Info_File, Hierarchy_File)
    GMT_Category <- inner_join(GMT, Category, by = 'pathway')
    Reactome_Table <-
      inner_join(X$Functionnal_Analysis$GSEA_Results, GMT_Category, by = 'pathway')
    Reactome_Table_filtered <-
      Reactome_Table %>%  filter(!is.na(padj)) %>%  filter(padj < p) %>%  mutate(Padj =
                                                                                   format(padj, scientific = TRUE, digits = 2))
    Reactome_Grid <- NULL
    Reactome_Grid <-
      ggplot(Reactome_Table_filtered,
             aes(
               Category,
               ES,
               colour = NES,
               text = pathway,
               text2 = Padj
             )) + geom_jitter(alpha = 0.7, size = (round(log(
               Reactome_Table_filtered$size
             ))) / 3) + facet_grid(Cluster ~ .,
                                   scales = "fixed",
                                   shrink = FALSE,
                                   space = "free_x") + theme_light() + theme(panel.grid.major.y = element_line(colour =
                                                                                                                 "gray")) + theme(axis.text.x = element_text(
                                                                                                                   face = "bold",
                                                                                                                   angle = 60,
                                                                                                                   hjust = 1,
                                                                                                                   size = 8
                                                                                                                 )) + scale_colour_gradient2(
                                                                                                                   low = "blue",
                                                                                                                   mid = "black",
                                                                                                                   high = "red",
                                                                                                                   midpoint = 0
                                                                                                                 )
    Reactome_Grid <-
      Reactome_Grid + geom_vline(xintercept = seq(1.5, length(
        unique(Reactome_Table_filtered$Category)
      ), 1), color = "gray")
    Reactome_Grid <- Reactome_Grid
    return(Reactome_Grid)
  }
