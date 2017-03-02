Functional_Analysis_GSEA <- function(X, GMTfile, nperm = 1000, maxSize = 1000, minSize = 0,
                                     nproc = 1, mode = "all") {
  df <- X$Dim_Red$Axis_Gene_Cor %>% separate(col = Genes, into = c("Genes", "bin"),
                                             sep = "-bin", convert = TRUE) %>% filter(bin == 1) %>% select(-bin)

  axis_rank <- df %>% select(-Genes) %>% purrr::map(abs) %>% purrr::map(multiply_by,
                                                                        -1) %>% purrr::map(set_names, nm = df$Genes) %>% purrr::map(rank, ties.method = "first") %>%
    purrr::map(sort)

  axis_gsea <- axis_rank %>% purrr::map(.f = function(x) {
    fgsea(pathways = GMTfile, stats = x, nperm = nperm, maxSize = maxSize, minSize = minSize,
          nproc = nproc, BPPARAM = SerialParam(), gseaParam = 1)
  }) %>% purrr::map(as_tibble)

  df2 <- X$cluster$Gene_Cluster_Distance %>% separate(col = Genes, into = c("Genes",
                                                                            "bin"), sep = "-bin", convert = TRUE) %>% filter(bin == 1) %>% select(-bin)

  cluster_rank <- df2 %>% select(-Genes) %>% purrr::map(abs) %>% purrr::map(multiply_by,
                                                                            -1) %>% purrr::map(set_names, nm = df$Genes) %>% purrr::map(rank, ties.method = "first") %>%
    purrr::map(sort)

  cluster_gsea <- cluster_rank %>% purrr::map(.f = function(x) {
    fgsea(pathways = GMTfile, stats = x, nperm = nperm, maxSize = maxSize, minSize = minSize,
          nproc = nproc, BPPARAM = SerialParam(), gseaParam = 1)
  }) %>% purrr::map(as_tibble)

  X$Functionnal_Analysis$GSEA_Results_Axis <- axis_gsea
  X$Functionnal_Analysis$RankingAxis <- axis_rank
  X$Functionnal_Analysis$GSEA_Results <- cluster_gsea
  X$Functionnal_Analysis$Ranking <- cluster_rank
  X$Functionnal_Analysis$GMTfile <- GMTfile
  X$Functionnal_Analysis$Pathways <- axis_gsea$Axis1$pathway
  X$Functionnal_Analysis$AllRanking <- append(X$Functionnal_Analysis$RankingAxis,
                                              X$Functionnal_Analysis$Ranking)
  X$Functionnal_Analysis$Shiny <- Create_Shiny_Functionnal_Analysis(X)
  cat(paste0("DONE\n"))
  class(X$Functionnal_Analysis) <- "FA.object"
  return(X)
}

Reactome_Category_Generator <- function(Info_File, Hierarchy_File) {
  colnames(Info_File)[1] = "value"
  MAS <- Hierarchy_File[grep("HSA", Hierarchy_File$X1), ]
  A <- !(Hierarchy_File$X1 %in% Hierarchy_File$X2)
  B <- Hierarchy_File$X1[A]
  C <- B[grep("HSA", B)] %>% unique %>% as_tibble
  D <- left_join(C, Info_File, by = "value") %>% select(1:2)
  X <- MAS$X2 %>% unique
  W <- vector(length = X %>% length)
  for (i in 1:(X %>% length)) {
    Y <- X[i]
    Flag <- Y %in% MAS$X2
    while (Flag == TRUE) {
      Y <- MAS$X1[grep(Y, MAS$X2)]
      Flag <- Y %in% MAS$X2
    }
    W[i] <- Y
  }
  Q <- W %>% as_tibble()
  colnames(Q) <- "value"
  Q <- Q %>% add_column(X)
  XC <- left_join(Q, D) %>% select(-value) %>% arrange(X2)
  colnames(XC)[1] <- "value"
  Standard <- bind_rows(XC, D)
  Final <- left_join(Standard, Info_File, by = "value") %>% select(-value)
  colnames(Final) <- c("Category", "pathway", "species")
  return(Final)
}
