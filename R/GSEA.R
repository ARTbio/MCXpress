Functional_Analysis_GSEA <-
  function(X,
           GMTfile,
           nperm = 1000,
           maxSize = 1000,
           minSize = 0,
           nproc = 1) {
    PCResults<- tibble()
    All_PCRanking<- tibble()
    for(i in 1:(min(X$Dim_Red$Eigen_Value%>% length,10)) ){
      cat(paste0("Processing PC", i, '...\n'))
      PCRanking<- X$Dim_Red$PC_Gene_Cor %>%
        filter(Component==paste0("PC",i)) %>%
        separate(col= Genes, into = c("Genes","bin"),sep = "-bin", convert = TRUE) %>%
        filter(bin==1) %>%
        mutate(AbsCor = Cor %>%  abs) %>%
        arrange(Cor) %>%
        mutate(Ranking=rank(Cor))
      GSEA_PCRank<-PCRanking$Ranking %>% as.matrix %>%  as.vector() %>%   set_names(value = PCRanking$Genes)
      GSEAResults <-
        fgsea(
          pathways = GMTfile,
          stats = GSEA_PCRank,
          nperm = nperm,
          maxSize = maxSize,
          minSize = minSize,
          nproc = nproc,
          BPPARAM = SerialParam(),
          gseaParam = 1
        )
      GSEAResults %<>%  as_tibble()%>% mutate(PC = rep(paste0("PC",i), GSEAResults %>%  nrow()))
      All_PCRanking <- bind_rows(All_PCRanking, PCRanking)
      PCResults <- bind_rows(PCResults, GSEAResults)
    }
    X$Functionnal_Analysis$GSEA_Results_PC<- PCResults[, c(ncol(PCResults), 1:(ncol(PCResults) - 1))]
    X$Functionnal_Analysis$GSEA_Results_PC[c("ES", "NES", "padj","pval")]<- X$Functionnal_Analysis$GSEA_Results_PC[c("ES", "NES", "padj","pval")] %>% format(scientific = TRUE, digits = 2)
    X$Functionnal_Analysis$RankingPC<- All_PCRanking %>%  select(-bin)
    cat("Beginning Gene Set Enrichment Analysis For Cluster...\n")
    Results <- tibble()
    All_Ranking <- tibble()
    Gene_All_Clusters_Distance <-
      X$cluster$Gene_Cluster_Distance %>% separate(col = Genes,
                                                   sep = '-bin',
                                                   into = c('Genes', 'bin')) %>% filter(bin == 1) %>% dplyr::select(-bin) %>%  gather(contains("Cluster"), key = Cluster, value =
                                                                                                                                        Distance)
    ID <- Gene_All_Clusters_Distance$Cluster %>% unique()
    for (i in 1:X$cluster$nClusters) {
      cat(paste0("Processing Cluster", i, '...\n'))
      Ranking <- NULL
      Gene_Cluster_Distance <-
        Gene_All_Clusters_Distance %>%  filter(Cluster == ID[i]) %>%  arrange(desc(Distance))
      Ranking_Table <-
        Gene_Cluster_Distance %>% select(Cluster, Genes, Distance) %>% mutate(Ranking =
                                                                                rank(Distance, ties.method = "random"))
      Ranking<-Ranking_Table$Ranking %>% rev
      names(Ranking) <- Ranking_Table$Genes
      GSEAResults <-
        fgsea(
          GMTfile,
          Ranking,
          nperm = nperm,
          maxSize = maxSize,
          minSize = minSize,
          nproc = nproc,
          BPPARAM = SerialParam(),
          gseaParam = 1
        )
      GSEAResults <-
        GSEAResults %>%  as_tibble()  %>% mutate(Cluster = rep(ID[i], GSEAResults %>%  nrow()))
      All_Ranking <- bind_rows(All_Ranking, Ranking_Table)
      Results <- bind_rows(Results, GSEAResults)
    }

    X$Functionnal_Analysis$GSEA_Results <-
      Results[, c(ncol(Results), 1:(ncol(Results) - 1))]
    X$Functionnal_Analysis$GSEA_Results[c("ES", "NES", "padj","pval")]<- X$Functionnal_Analysis$GSEA_Results[c("ES", "NES", "padj","pval")] %>% format(scientific = TRUE, digits = 2)
    X$Functionnal_Analysis$Ranking  <-  All_Ranking %>% mutate(Ranking=Ranking %>% rev)
    X$Functionnal_Analysis$Distance_Var <-
      X$cluster$Gene_Cluster_Distance %>% gather(-Genes, key = "Cluster", value = "Distance") %>% group_by(Genes) %>%  summarise(Dis_Var =
                                                                                                                                   var(Distance)) %>% mutate(Rank_Var = -Dis_Var %>% rank)

    X$Functionnal_Analysis$Shiny  <-Create_Shiny_Functionnal_Analysis(X)
    A<-X$Functionnal_Analysis$Ranking %>%  select(-Distance)
    B<-X$Functionnal_Analysis$RankingPC %>%  select(-Cor, -AbsCor) %>%  extract(,c(2,1,3))
    colnames(A)[1]<-"Group"
    colnames(B)[1]<-"Group"
    Grouped<-bind_rows(A,B)
    X$Functionnal_Analysis$Grouped<-Grouped
    X$Functionnal_Analysis$Info <-
      paste0("number of permutation:",
             nperm,
             "maxSize:",
             maxSize,
             "minSize:",
             minSize)
    cat(paste0('DONE\n'))
    class(X$Functionnal_Analysis) <- "FA_Object"
    return(X)
  }

Reactome_Category_Generator <- function(Info_File, Hierarchy_File) {
  colnames(Info_File)[1] = 'value'
  MAS <-  Hierarchy_File[grep("HSA", Hierarchy_File$X1), ]
  A   <-  !(Hierarchy_File$X1 %in% Hierarchy_File$X2)
  B   <-  Hierarchy_File$X1[A]
  C   <-  B[grep("HSA", B)] %>%  unique %>%  as_tibble
  D   <-  left_join(C, Info_File, by = "value") %>%  select(1:2)
  X   <-  MAS$X2 %>%  unique
  W   <-  vector(length = X %>%  length)
  for (i in 1:(X %>% length)) {
    Y <- X[i]
    Flag <- Y %in% MAS$X2
    while (Flag == TRUE) {
      Y <- MAS$X1[grep(Y, MAS$X2)]
      Flag <- Y %in% MAS$X2
    }
    W[i] <- Y
  }
  Q <- W %>%  as_tibble()
  colnames(Q) <- 'value'
  Q <- Q %>%  add_column(X)
  XC <- left_join(Q, D) %>% select(-value) %>% arrange(X2)
  colnames(XC)[1] <- 'value'
  Standard <- bind_rows(XC, D)
  Final <-
    left_join(Standard, Info_File, by = 'value') %>%  select(-value)
  colnames(Final) <- c('Category', 'pathway', 'species')
  return(Final)
}


