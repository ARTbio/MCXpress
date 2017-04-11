
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param X DESCRIPTION.
#' @param min_reads DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
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
  name <- statsAdj[as.vector(na.omit(match(pathway, names(statsAdj))))] %>%
    sort %>% names %>% rev
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

  data_markers <- data.frame(Genes = rep(name, each = 100),
    x = rep(pathway, 100) %>% sort, y = rep(seq(from = diff/2,
      to = -diff/2, length.out = 100), pathway %>% length))
  plot_ly(data = toPlot) %>% add_lines(x = ~x, y = ~y, mode = "lines",
    line = list(color = "rgb(154, 240, 24)", width = 2),
    name = "Enrichment") %>% add_lines(x = ~x, y = ~min(bottoms),
    line = list(color = "rgb(110, 193, 248)", width = 2,
      dash = "dash"), name = "Lower limit", hoverinfo = "text",
    text = ~round(min(bottoms), digits = 4)) %>% add_lines(x = ~x,
    y = ~max(tops), line = list(color = "rgb(250, 150, 10)",
      width = 2, dash = "dash"), name = "Upper limit",
    hoverinfo = "text", text = ~round(max(tops), digits = 4)) %>%
    add_markers(data = data_markers, x = ~x, y = ~y, marker = list(color = "rgb(0, 10, 10)",
      size = 1), name = "Genes", hoverinfo = "text", text = ~paste0(Genes,
      "</br>", "Rank:", x), showlegend = FALSE)
}


plot_immune <- function(X, GMTfile, quali, p = 0.05, thres =0.5){
  GMT <- GMTfile %>% attributes() %>% as_tibble
  colnames(GMT) <- "pathway"
  Category <- quali
  GMT_Category <- inner_join(GMT, Category, by = "pathway")
  Reactome_Table <- inner_join(X$Functionnal_Analysis$GSEA_Results %>%
    map2(.y = X$Functionnal_Analysis$GSEA_Results %>%
      names, .f = function(x, y) {
      x %>% mutate(Cluster = y)
    }) %>% bind_rows(), GMT_Category, by = "pathway")
  Reactome_Table_filtered <- Reactome_Table %>% filter(!is.na(padj)) %>%
    filter(padj < p) %>%  filter(ES > thres) %>% mutate(Padj = format(padj, scientific = TRUE,
    digits = 2))


  # Reactome_Table_filtered$Cluster <- Reactome_Table_filtered$Cluster %>%
  #   factor(levels = c("STHSC", "LTHSC", "MPP", "CMP", "MEP",
  #     "GMP", "LMPP"))

  Reactome_Grid <- NULL
  Reactome_Grid <- ggplot(Reactome_Table_filtered, aes(Category,
    ES, colour = ES, text1 = pathway, text2 = Padj, text3 = NES,
    text4 = Category)) + geom_jitter(alpha = 0.7, size = (round(log(Reactome_Table_filtered$size)))/3) +
    facet_grid(Cluster ~ Subcategory, scales = "fixed", shrink = FALSE, drop = FALSE) + theme_light() + theme(panel.grid.major.y = element_line(colour = "gray")) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold",
      angle = 60, hjust = 1, size = 8)) + scale_colour_gradient2(low = "blue",
    mid = "black", high = "red")
  if(Reactome_Table_filtered$Subcategory %>% unique %>% length>1){
  Reactome_Grid <- Reactome_Grid + geom_vline(xintercept = seq(1.5,
    length(unique(Reactome_Table_filtered$Subcategory)), 1),
    color = "gray")}
  Reactome_Grid %<>% ggplotly(tooltip = c("x","text1"))%>% layout(margin=list(
    b=100,
    t=75, l=75, r=75, pad=0
))
  return(Reactome_Grid)
}



#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param X DESCRIPTION.
#' @param Info_File DESCRIPTION.
#' @param Hierarchy_File DESCRIPTION.
#' @param GMTfile DESCRIPTION.
#' @param p adjusted p-value threshold.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
Reactome_GridPlot <- function(X, Info_File, Hierarchy_File, GMTfile,
  p = 0.05, width=1000, height=1000) {
  GMT <- GMTfile %>% attributes() %>% as_tibble
  colnames(GMT) <- "pathway"
  Category <- Reactome_Category_Generator2(Info_File, Hierarchy_File)
  GMT_Category <- inner_join(GMT, Category, by = "pathway")
  Reactome_Table <- inner_join(X$Functionnal_Analysis$GSEA_Results %>%
    map2(.y = X$Functionnal_Analysis$GSEA_Results %>%
      names, .f = function(x, y) {
      x %>% mutate(Cluster = y)
    }) %>% bind_rows(), GMT_Category, by = "pathway")
  Reactome_Table_filtered <- Reactome_Table %>% filter(!is.na(padj)) %>%
    filter(padj < p) %>% mutate(Padj = format(padj, scientific = TRUE,
    digits = 2))

Reactome_Category_Generator(Info_File = file2 , Hierarchy_File = file1)
  # Reactome_Table_filtered$Cluster <- Reactome_Table_filtered$Cluster %>%
  #   factor(levels = c("STHSC", "LTHSC", "MPP", "CMP", "MEP",
  #     "GMP", "LMPP"))

  Reactome_Grid <- NULL
  Reactome_Grid <- ggplot(Reactome_Table_filtered, aes(Category,
    NES, colour = ES, text1 = pathway, text2 = Padj, text3 = NES,
    text4 = Category)) + geom_jitter(alpha = 0.7, size = (round(log(Reactome_Table_filtered$size)))/3) +
    facet_grid(Cluster ~ Category., scales = "fixed", shrink = FALSE,
      space = "free_x") + theme_light() + theme(panel.grid.major.y = element_line(colour = "gray")) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(face = "bold",
      angle = 60, hjust = 1, size = 8)) + scale_colour_gradient2(low = "blue",
    mid = "black", high = "red")
  Reactome_Grid <- Reactome_Grid + geom_vline(xintercept = seq(1.5,
    length(unique(Reactome_Table_filtered$Category)), 1),
    color = "gray")
  Reactome_Grid %<>% ggplotly(tooltip = c("text1", "text2", "text3"), autosize = F, width = width, height = height) %>% layout(margin=list(
    b=250,
    t=25, l=50, r=50, pad=0
))
  return(Reactome_Grid)
}

Reactome_Category_Generator <- function(Info_File, Hierarchy_File) {
  colnames(Info_File)[1] = "value"
  MAS <- Hierarchy_File[grep("MMU", Hierarchy_File$X1), ]
  A <- !(Hierarchy_File$X1 %in% Hierarchy_File$X2)
  B <- Hierarchy_File$X1[A]
  C <- B[grep("MMU", B)] %>% unique %>% as_tibble
  D <- left_join(C, Info_File, by = "value") %>% select(1:2)
  X <- MAS$X2 %>% unique
  W <- vector(length = X %>% length)
  for (i in 1:(X %>% length)) {
    Y <- X[i]
    Flag <- Y %in% MAS$X2
    while (Flag == TRUE) {
      Y <- MAS$X1[grepl(Y, MAS$X2) %>% which %>%  head(1)]
      Flag <- Y %in% MAS$X2
    }
    W[i] <- Y
  }
  Q <- W %>% as_tibble()
  colnames(Q) <- "value"
  Q <- Q %>% add_column(X)
  XC <- left_join(Q, D, by="value") %>% select(-value) %>% arrange(X2)
  colnames(XC)[1] <- "value"
  Standard <- bind_rows(XC, D)
  Final <- left_join(Standard, Info_File, by = "value") %>%
    select(-value)
  colnames(Final) <- c("Category", "pathway", "species")
  return(Final)
}



# file1<-read_tsv("C:/Users/Akira/Documents/Reactome/ReactomePathwaysRelation.txt", col_names = FALSE)
# file2<-read_tsv("C:/Users/Akira/Documents/Reactome/ReactomePathways.txt", col_names = FALSE)
# file1
# Info_File<- file2
# Hierarchy_File<-file1
# Reactome_Category_Generator2 <- function(Info_File,
#     Hierarchy_File) {
#     colnames(Info_File)[1] = "value"
#     MAS <- Hierarchy_File[grep("MMU", Hierarchy_File$X1),
#         ]
#     A <- !(Hierarchy_File$X1 %in% Hierarchy_File$X2)
#     B <- Hierarchy_File$X1[A]
#     C <- B[grep("MMU", B)] %>% unique %>% as_tibble
#     D <- left_join(C, Info_File, by = "value") %>%
#         select(1:2)
#     X <- MAS$X2 %>% unique
#     W <- vector(length = X %>% length)
#     for (i in 1:(X %>% length)) {
#         Y <- X[i]
#         Flag <- Y %in% MAS$X2
#         while (Flag == TRUE) {
#             Y <- MAS$X1[grep(Y, MAS$X2)]
#             Flag <- Y %in% MAS$X2
#         }
#         W[i] <- Y
#     }
#     Q <- W %>% as_tibble()
#     colnames(Q) <- "value"
#     Q <- Q %>% add_column(X)
#     XC <- left_join(Q, D) %>% select(-value) %>% arrange(X2)
#     colnames(XC)[1] <- "value"
#     Standard <- bind_rows(XC, D)
#     Final <- left_join(Standard, Info_File, by = "value") %>%
#         select(-value)
#     colnames(Final) <- c("Category", "pathway", "species")
#     return(Final)
# }
