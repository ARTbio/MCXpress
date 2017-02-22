Functional_Analysis_GSEA <-
  function(X,
           GMTfile,
           nperm = 1000,
           maxSize = 1000,
           minSize = 0,
           nproc = 1,
           mode = "all"){
    AxisResults<- tibble()
    All_AxisRanking<- tibble()
    for(i in 1:(min(X$Dim_Red$Cells_Principal %>% ncol ,5))){
      cat(paste0("Processing Axis", i, '...\n'))
      #Create Ranking to feed in the FGSEA
      AxisRanking<- X$Dim_Red$Axis_Gene_Cor %>%
        filter(Axis==paste0("Axis",i)) %>%
        separate(col= Genes, into = c("Genes","bin"),sep = "-bin", convert = TRUE) %>%
        filter(bin==1) %>%
        mutate(AbsCor = Cor %>%  abs) %>%
        arrange(desc(AbsCor)) %>%
        mutate(Ranking=rank(-AbsCor))
      GSEA_AxisRank<-AxisRanking$Ranking %>% as.matrix %>%  as.vector() %>%   set_names(value = AxisRanking$Genes)
      #End Ranking Creation for Axis
      GSEAResults <-
        fgsea(
          pathways = GMTfile,
          stats = GSEA_AxisRank,
          nperm = nperm,
          maxSize = maxSize,
          minSize = minSize,
          nproc = nproc,
          BPPARAM = SerialParam(),
          gseaParam = 1
        )
      GSEAResults %<>%
        as_tibble() %>%
        mutate(Axis = rep(paste0("Axis",i), GSEAResults %>%  nrow()))
      All_AxisRanking %<>% bind_rows(AxisRanking)
      AxisResults %<>% bind_rows(GSEAResults)
    }
    X$Functionnal_Analysis$GSEA_Results_Axis<- AxisResults[, c(ncol(AxisResults), 1:(ncol(AxisResults) - 1))]
    X$Functionnal_Analysis$GSEA_Results_Axis[c("ES", "NES", "padj","pval")]<- X$Functionnal_Analysis$GSEA_Results_Axis[c("ES", "NES", "padj","pval")] %>% format(scientific = TRUE, digits = 2)
    X$Functionnal_Analysis$RankingAxis<- All_AxisRanking %>%  select(-bin)

    cat("Beginning Gene Set Enrichment Analysis For Cluster...\n")
    Results <- tibble()
    All_Ranking <- tibble()

    if(mode=="all")
    {Gene_All_Clusters_Distance<-X$cluster$Gene_Cluster_Distance %>%
      select(Genes,contains("Cluster")) %>%
      separate(Genes, into=c("Genes", "bin"), sep="-bin") %>%
      group_by(Genes)%>%
      summarise_at(vars(contains("Cluster")), min) %>%
      gather(contains("Cluster"), key = Cluster, value = Distance)}
    if(mode=="bin1"){
    Gene_All_Clusters_Distance <-
      X$cluster$Gene_Cluster_Distance %>%
      separate(col = Genes, sep = '-bin', into = c('Genes', 'bin')) %>%
      filter(bin == 1) %>%
      select(-bin) %>%
      gather(contains("Cluster"), key = Cluster, value =Distance)}

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


    A<-X$Functionnal_Analysis$Ranking %>%  select(-Distance)
    B<-X$Functionnal_Analysis$RankingAxis %>%  select(-Cor, -AbsCor) %>%  select(Axis, everything())
    colnames(A)[1]<-"Group"
    colnames(B)[1]<-"Group"
    X$Functionnal_Analysis$Grouped<-bind_rows(A,B)
    X$Functionnal_Analysis$Shiny<-Create_Shiny_Functionnal_Analysis(X)
    cat(paste0('DONE\n'))
    class(X$Functionnal_Analysis) <- "FA.object"
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


