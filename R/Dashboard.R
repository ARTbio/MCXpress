#' Create Dashboard
#'
#' @param X MCXpress object
#' @return Shiny object
create_dashboard1 <- function(X)
{
  dr_axis <-
    X$MCA$cells_principal %>% select(contains("Axis")) %>% colnames
  ui <-
    dashboardPage(skin = sample(c("red", "blue", "yellow", "green"), size = 1),
      dashboardHeader(title = h1("MCXpress")),
      dashboardSidebar(sidebarMenu(
        menuItem("MCA",
                 tabName = "mca", icon = icon("arrows"))
      )),
      dashboardBody(CSS, tabItems(TabMCA(X,
                                         dr_axis)))
    )
  server <- function(input, output, clientData, session)
  {
    options(warn = -1)
    MCA_axis_name <-
      X$MCA$cells_principal %>% select(contains("Axis")) %>% colnames
    output$CellSpaceGeneCor <- renderPlotly({
      if (input$MCA_CS_AC_Type == "Principal")
      {
        axis_cor <-
          X$MCA$cells_principal %>% rownames_to_column(var = "Sample") %>%
          inner_join(
            X$ExpressionMatrix[input$MCA_CS_AC_Gene, ] %>% data.frame() %>%
              tibble::rownames_to_column() %>% set_colnames(c("Sample", "Expression")),
            by = "Sample"
          )
        p <-
          plot_ly(data = axis_cor,
                  x = ~ axis_cor[[input$MCA_CS_AC_Axis_x]],
                  y = ~ axis_cor[[input$MCA_CS_AC_Axis_y]]) %>% add_markers(
                    text = ~ Sample,
                    hoverinfo = "text",
                    alpha = input$MCA_CS_AC_Alpha,
                    color = ~ Expression,
                    marker = list(size = input$MCA_CS_AC_Size)
                  ) %>% layout(
                    xaxis = list(title = input$MCA_CS_AC_Axis_x),
                    yaxis = list(title = input$MCA_CS_AC_Axis_y)
                  )
        p
      } else
      {
        axis_cor <-
          X$MCA$cells_standard %>% rownames_to_column(var = "Sample") %>%
          inner_join(
            X$ExpressionMatrix[input$MCA_CS_AC_Gene, ] %>% data.frame() %>%
              tibble::rownames_to_column() %>% set_colnames(c("Sample", "Expression")),
            by = "Sample"
          )
        p <-
          plot_ly(data = axis_cor,
                  x = ~ axis_cor[[input$MCA_CS_AC_Axis_x]],
                  y = ~ axis_cor[[input$MCA_CS_AC_Axis_y]]) %>% add_markers(
                    text = ~ Sample,
                    hoverinfo = "text",
                    alpha = input$MCA_CS_AC_Alpha,
                    color = ~ Expression,
                    marker = list(size = input$MCA_CS_AC_Size)
                  ) %>% layout(
                    xaxis = list(title = input$MCA_CS_Axis_x),
                    yaxis = list(title = input$MCA_CS_Axis_y)
                  )
        p
      }
    })
    output$CellSpace2D <- renderPlotly({
      if (input$MCA_CS_2D_Type == "Principal")
      {
        d3 <- X$MCA$cells_principal %>% rownames_to_column(var = "Sample")
        p <-
          plot_ly(
            data = d3,
            x = ~ d3[[input$MCA_CS_Axis_x]],
            y = ~ d3[[input$MCA_CS_Axis_y]],
            type = "scatter",
            mode = "markers",
            text = ~ Sample,
            hoverinfo = "text",
            alpha = input$Alpha,
            marker = list(size = input$Size)
          ) %>% layout(
            xaxis = list(title = input$MCA_CS_Axis_x),
            yaxis = list(title = input$MCA_CS_Axis_y)
          )
        p
      } else
      {
        d3 <- X$MCA$cells_standard %>% rownames_to_column(var = "Sample")
        p <-
          plot_ly(
            data = d3,
            x = ~ d3[[input$MCA_CS_Axis_x]],
            y = ~ d3[[input$MCA_CS_Axis_y]],
            type = "scatter",
            mode = "markers",
            text = ~ Sample,
            hoverinfo = "text",
            alpha = input$Alpha,
            marker = list(size = input$Size)
          ) %>% layout(
            xaxis = list(title = input$MCA_CS_Axis_x),
            yaxis = list(title = input$MCA_CS_Axis_y)
          )
        p
      }
    })

    output$CellSpace3D <-
      renderPlotly(if (input$MCA_CS_3D_Type == "Standard")
      {
        plot_ly(
          X$MCA$cells_standard,
          mode = "markers",
          text = ~ paste(
            rownames(X$MCA$cells_standard),
            "</br>",
            input$MCA_CS_Axis1_3D,
            ": ",
            X$MCA$cells_standard[[input$MCA_CS_Axis1_3D]] %>%
              signif(digits = 4),
            "</br>",
            input$Axis2_3D,
            ": ",
            X$MCA$cells_standard[[input$MCA_CS_Axis2_3D]] %>%
              signif(digits = 4),
            "</br>",
            input$Axis3_3D,
            ": ",
            X$MCA$cells_standard[[input$MCA_CS_Axis3_3D]] %>%
              signif(digits = 4)
          ),
          x = ~ X$MCA$cells_standard[[input$MCA_CS_Axis1_3D]],
          y = ~ X$MCA$cells_standard[[input$MCA_CS_Axis2_3D]],
          z = ~ X$MCA$cells_standard[[input$MCA_CS_Axis3_3D]],
          hoverinfo = "text",
          marker = list(
            opacity = input$MCA_CS_Alpha_3D,
            size = input$MCA_CS_Size_3D,
            line = list(color = "rgba(0, 0, 0, .8)",
                        width = 2)
          )
        ) %>% add_markers() %>% layout(
          autosize = FALSE,
          scene = list(
            xaxis = list(title = input$MCA_CS_Axis1_3D),
            yaxis = list(title = input$MCA_CS_Axis2_3D),
            zaxis = list(title = input$MCA_CS_Axis3_3D)
          )
        )
      } else
      {
        plot_ly(
          X$MCA$cells_standard,
          mode = "markers",
          text = ~ paste(
            rownames(X$MCA$cells_principal),
            "</br>",
            input$MCA_CS_Axis1_3D,
            ": ",
            X$MCA$cells_principal[[input$MCA_CS_Axis1_3D]] %>%
              signif(digits = 4),
            "</br>",
            input$MCA_CS_Axis2_3D,
            ": ",
            X$MCA$cells_principal[[input$MCA_CS_Axis2_3D]] %>%
              signif(digits = 4),
            "</br>",
            input$MCA_CS_Axis3_3D,
            ": ",
            X$MCA$cells_principal[[input$MCA_CS_Axis3_3D]] %>%
              signif(digits = 4)
          ),
          x = ~ X$MCA$cells_principal[[input$MCA_CS_Axis1_3D]],
          y = ~ X$MCA$cells_principal[[input$MCA_CS_Axis2_3D]],
          z = ~ X$MCA$cells_principal[[input$MCA_CS_Axis3_3D]],
          hoverinfo = "text",
          marker = list(
            opacity = input$MCA_CS_Alpha_3D,
            size = input$MCA_CS_Size_3D,
            line = list(color = "rgba(0, 0, 0, .8)",
                        width = 2)
          )
        ) %>% add_markers() %>% layout(
          margin = list(
            b = 50,
            t = 50,
            l = 50,
            r = 100,
            pad = 0
          ),
          autosize = FALSE,
          scene = list(
            xaxis = list(title = input$MCA_CS_Axis1_3D),
            yaxis = list(title = input$MCA_CS_Axis2_3D),
            zaxis = list(title = input$MCA_CS_Axis3_3D)
          )
        )
      })
    output$TableGeneCor <-
      DT::renderDataTable(X$MCA$Axis_Gene_Cor %>% extract(1:6))
    output$GeneSpace <- renderPlotly({
      if (input$Type_Gene == "Principal")
      {
        d3 <- X$MCA$genes_principal %>% rownames_to_column(var = "Genes")
      } else
      {
        d3 <- X$MCA$genes_standard %>% rownames_to_column(var = "Genes")
      }
      p <-
        plot_ly(
          data = d3,
          x = ~ d3[[input$Axis1_Gene]],
          y = ~ d3[[input$Axis2_Gene]],
          type = "scatter",
          mode = "markers",
          text = ~ Genes,
          hoverinfo = "text",
          alpha = input$Alpha_Gene,
          marker = list(size = input$Size_Gene)
        ) %>%
        layout(
          xaxis = list(title = input$Axis1_Gene),
          yaxis = list(title = input$Axis2_Gene)
        )
      p
    })
    output$Eigen <- renderPlotly({
      Shiny_Eigen <- X$MCA$explained_eigen_variance
      Shiny_Eigen <- Shiny_Eigen[1:input$amount_adjust, ]
      plot_ly(Shiny_Eigen) %>% add_bars(
        x = ~ Axis,
        y = ~ Shiny_Eigen[[input$Mode]],
        hoverinfo = "text",
        text = ~ paste0(
          "Axis: ",
          Axis,
          "</br>Explained Variance: ",
          Explained_Variance %>% round(3),
          "</br>Cumulated Explained Variance: ",
          Cumulative %>% round(3)
        )
      ) %>% layout(
        yaxis = list(title = input$Mode),
        margin = list(
          b = 100,
          t = 25,
          l = 50,
          r = 50,
          pad = 0
        )
      )
    })
  }
  return(shinyApp(ui, server))
}

#' Create Dashboard
#'
#' @param X MCXpress object
#' @return Shiny object
create_dashboard2 <- function(X)
{
  dr_axis <-
    X$MCA$cells_principal %>% select(contains("Axis")) %>% colnames
  ui <-
    dashboardPage(skin = sample(c("red", "blue", "yellow", "green"), size = 1),
      dashboardHeader(title = h1("MCXpress")),
      dashboardSidebar(sidebarMenu(
        menuItem("MCA",
                 tabName = "mca", icon = icon("arrows")),
        menuItem(
          "Clustering",
          tabName = "clus",
          icon = icon("object-group")
        )
      )),
      dashboardBody(CSS, tabItems(TabMCA(X, dr_axis),
                                  TabClus(X, dr_axis)))
    )

  server <- function(input, output, clientData, session)
  {
    MCA_axis_name <-
      X$MCA$cells_principal %>% select(contains("Axis")) %>% colnames

    output$CellSpaceGeneCor <- renderPlotly({
      if (input$MCA_CS_AC_Type == "Principal")
      {
        axis_cor <-
          X$MCA$cells_principal %>% rownames_to_column(var = "Sample") %>%
          inner_join(
            X$ExpressionMatrix[input$MCA_CS_AC_Gene, ] %>% data.frame() %>%
              tibble::rownames_to_column() %>% set_colnames(c("Sample", "Expression")),
            by = "Sample"
          )
        p <-
          plot_ly(data = axis_cor,
                  x = ~ axis_cor[[input$MCA_CS_AC_Axis_x]],
                  y = ~ axis_cor[[input$MCA_CS_AC_Axis_y]]) %>% add_markers(
                    text = ~ Sample,
                    hoverinfo = "text",
                    alpha = input$MCA_CS_AC_Alpha,
                    color = ~ Expression,
                    marker = list(size = input$MCA_CS_AC_Size)
                  ) %>% layout(
                    xaxis = list(title = input$MCA_CS_AC_Axis_x),
                    yaxis = list(title = input$MCA_CS_AC_Axis_y)
                  )
        p
      } else
      {
        axis_cor <-
          X$MCA$cells_standard %>% rownames_to_column(var = "Sample") %>%
          inner_join(
            X$ExpressionMatrix[input$MCA_CS_AC_Gene, ] %>% data.frame() %>%
              tibble::rownames_to_column() %>% set_colnames(c("Sample", "Expression")),
            by = "Sample"
          )
        p <-
          plot_ly(data = axis_cor,
                  x = ~ axis_cor[[input$MCA_CS_AC_Axis_x]],
                  y = ~ axis_cor[[input$MCA_CS_AC_Axis_y]]) %>% add_markers(
                    text = ~ Sample,
                    hoverinfo = "text",
                    alpha = input$MCA_CS_AC_Alpha,
                    color = ~ Expression,
                    marker = list(size = input$MCA_CS_AC_Size)
                  ) %>% layout(
                    xaxis = list(title = input$MCA_CS_Axis_x),
                    yaxis = list(title = input$MCA_CS_Axis_y)
                  )
        p
      }
    })
    output$CellSpace2D <- renderPlotly({
      if (input$MCA_CS_2D_Type == "Principal")
      {
        d3 <- X$MCA$cells_principal %>% rownames_to_column(var = "Sample")
        p <-
          plot_ly(
            data = d3,
            x = ~ d3[[input$MCA_CS_Axis_x]],
            y = ~ d3[[input$MCA_CS_Axis_y]],
            type = "scatter",
            mode = "markers",
            text = ~ Sample,
            hoverinfo = "text",
            alpha = input$Alpha,
            marker = list(size = input$Size)
          ) %>% layout(
            xaxis = list(title = input$MCA_CS_Axis_x),
            yaxis = list(title = input$MCA_CS_Axis_y)
          )
        p
      } else
      {
        d3 <- X$MCA$cells_standard %>% rownames_to_column(var = "Sample")
        p <-
          plot_ly(
            data = d3,
            x = ~ d3[[input$MCA_CS_Axis_x]],
            y = ~ d3[[input$MCA_CS_Axis_y]],
            type = "scatter",
            mode = "markers",
            text = ~ Sample,
            hoverinfo = "text",
            alpha = input$Alpha,
            marker = list(size = input$Size)
          ) %>% layout(
            xaxis = list(title = input$MCA_CS_Axis_x),
            yaxis = list(title = input$MCA_CS_Axis_y)
          )
        p
      }
    })

    output$CellSpace3D <-
      renderPlotly(if (input$MCA_CS_3D_Type == "Standard")
      {
        plot_ly(
          X$MCA$cells_standard,
          mode = "markers",
          text = ~ paste(
            rownames(X$MCA$cells_standard),
            "</br>",
            input$MCA_CS_Axis1_3D,
            ": ",
            X$MCA$cells_standard[[input$MCA_CS_Axis1_3D]] %>%
              signif(digits = 4),
            "</br>",
            input$Axis2_3D,
            ": ",
            X$MCA$cells_standard[[input$MCA_CS_Axis2_3D]] %>%
              signif(digits = 4),
            "</br>",
            input$Axis3_3D,
            ": ",
            X$MCA$cells_standard[[input$MCA_CS_Axis3_3D]] %>%
              signif(digits = 4)
          ),
          x = ~ X$MCA$cells_standard[[input$MCA_CS_Axis1_3D]],
          y = ~ X$MCA$cells_standard[[input$MCA_CS_Axis2_3D]],
          z = ~ X$MCA$cells_standard[[input$MCA_CS_Axis3_3D]],
          hoverinfo = "text",
          marker = list(
            opacity = input$MCA_CS_Alpha_3D,
            size = input$MCA_CS_Size_3D,
            line = list(color = "rgba(0, 0, 0, .8)",
                        width = 2)
          )
        ) %>% add_markers() %>% layout(
          autosize = FALSE,
          scene = list(
            xaxis = list(title = input$MCA_CS_Axis1_3D),
            yaxis = list(title = input$MCA_CS_Axis2_3D),
            zaxis = list(title = input$MCA_CS_Axis3_3D)
          )
        )
      } else
      {
        plot_ly(
          X$MCA$cells_standard,
          mode = "markers",
          text = ~ paste(
            rownames(X$MCA$cells_principal),
            "</br>",
            input$MCA_CS_Axis1_3D,
            ": ",
            X$MCA$cells_principal[[input$MCA_CS_Axis1_3D]] %>%
              signif(digits = 4),
            "</br>",
            input$MCA_CS_Axis2_3D,
            ": ",
            X$MCA$cells_principal[[input$MCA_CS_Axis2_3D]] %>%
              signif(digits = 4),
            "</br>",
            input$MCA_CS_Axis3_3D,
            ": ",
            X$MCA$cells_principal[[input$MCA_CS_Axis3_3D]] %>%
              signif(digits = 4)
          ),
          x = ~ X$MCA$cells_principal[[input$MCA_CS_Axis1_3D]],
          y = ~ X$MCA$cells_principal[[input$MCA_CS_Axis2_3D]],
          z = ~ X$MCA$cells_principal[[input$MCA_CS_Axis3_3D]],
          hoverinfo = "text",
          marker = list(
            opacity = input$MCA_CS_Alpha_3D,
            size = input$MCA_CS_Size_3D,
            line = list(color = "rgba(0, 0, 0, .8)",
                        width = 2)
          )
        ) %>% add_markers() %>% layout(
          margin = list(
            b = 50,
            t = 50,
            l = 50,
            r = 100,
            pad = 0
          ),
          autosize = FALSE,
          scene = list(
            xaxis = list(title = input$MCA_CS_Axis1_3D),
            yaxis = list(title = input$MCA_CS_Axis2_3D),
            zaxis = list(title = input$MCA_CS_Axis3_3D)
          )
        )
      })
    output$TableGeneCor <-
      DT::renderDataTable(X$MCA$Axis_Gene_Cor %>% extract(1:6))
    output$GeneSpace <- renderPlotly({
      if (input$Type_Gene == "Principal")
      {
        d3 <- X$MCA$genes_principal %>% rownames_to_column(var = "Genes")
      } else
      {
        d3 <- X$MCA$genes_standard %>% rownames_to_column(var = "Genes")
      }
      p <-
        plot_ly(
          data = d3,
          x = ~ d3[[input$Axis1_Gene]],
          y = ~ d3[[input$Axis2_Gene]],
          type = "scatter",
          mode = "markers",
          text = ~ Genes,
          hoverinfo = "text",
          alpha = input$Alpha_Gene,
          marker = list(size = input$Size_Gene)
        ) %>%
        layout(
          xaxis = list(title = input$Axis1_Gene),
          yaxis = list(title = input$Axis2_Gene)
        )
      p
    })
    output$Eigen <- renderPlotly({
      Shiny_Eigen <- X$MCA$explained_eigen_variance
      Shiny_Eigen <- Shiny_Eigen[1:input$amount_adjust, ]
      plot_ly(Shiny_Eigen) %>% add_bars(
        x = ~ Axis,
        y = ~ Shiny_Eigen[[input$Mode]],
        hoverinfo = "text",
        text = ~ paste0(
          "Axis: ",
          Axis,
          "</br>Explained Variance: ",
          Explained_Variance %>% round(3),
          "</br>Cumulated Explained Variance: ",
          Cumulative %>% round(3)
        )
      ) %>% layout(
        yaxis = list(title = input$Mode),
        margin = list(
          b = 100,
          t = 25,
          l = 50,
          r = 50,
          pad = 0
        )
      )
    })


    # ____________________________________________________________________________
    # Clus server ####

    options(warn = -1)

    Boxplot <-
      X$ExpressionMatrix %>% data.frame %>% rownames_to_column(var = "Genes") %>%
      set_colnames(c("Genes", X$ExpressionMatrix %>% colnames)) %>% gather("Sample",
                                                                           "Expression", -Genes) %>% as_tibble %>% dplyr::arrange(Sample) %>% inner_join(X$cluster$labels,

                                                                                                                                                         by = "Sample") %>% as_tibble
    output$CellSpace_Clus <- renderPlotly({
      d3 <-
        X$MCA$cells_principal %>% rownames_to_column(var = "Sample") %>%
        inner_join(X$cluster$labels, by = "Sample")
      d4 <-
        X$MCA$cells_standard %>% rownames_to_column(var = "Sample") %>%
        inner_join(X$cluster$labels, by = "Sample")

      if (input$Type == "Principal")
      {
        p <-
          plot_ly(
            data = d3,
            x = ~ d3[[input$Axis1_Clus]],
            y = ~ d3[[input$Axis2_Clus]],
            color = ~ Cluster,
            type = "scatter",
            mode = "markers",
            text = ~ Sample,
            hoverinfo = "text",
            alpha = input$Alpha_Clus,
            marker = list(size = input$Size_Clus)
          ) %>%
          layout(
            xaxis = list(title = input$Axis1_Clus),
            yaxis = list(title = input$Axis2_Clus)
          )
        p
      } else
      {
        p <-
          plot_ly(
            data = d4,
            x = ~ d4[[input$Axis1_Clus]],
            y = ~ d4[[input$Axis2_Clus]],
            color = ~ Cluster,
            type = "scatter",
            mode = "markers",
            text = ~ Sample,
            hoverinfo = "text",
            alpha = input$Alpha_Clus,
            marker = list(size = input$Size_Clus)
          ) %>%
          layout(
            xaxis = list(title = input$Axis1_Clus),
            yaxis = list(title = input$Axis2_Clus)
          )
        p
      }
    })

    output$GeneSpace_Clus <- renderPlotly({
      Genes <-
        X$MCA$genes_standard %>% rownames_to_column(var = "Genes") %>%
        select_("Genes", input$Axis1_Gene_Clus, input$Axis2_Gene_Clus) %>%
        set_colnames(c("Genes", "AP1", "AP2"))
      Centroids <-
        X$cluster$coord_centroids %>% select_("Cluster",
                                              input$Axis1_Gene_Clus,
                                              input$Axis2_Gene_Clus) %>% set_colnames(c("Cluster", "AC1", "AC2"))
      p <-
        plot_ly(data = Genes,
                x = ~ AP1,
                y = ~ AP2) %>% add_markers(
                  name = "Genes",
                  text = ~ Genes,
                  hoverinfo = "text",
                  marker = list(
                    size = input$Size_Gene_Clus,
                    color = "black",
                    alpha = input$Alpha_Gene_Clus
                  )
                ) %>% add_markers(
                  data = Centroids,
                  x = ~ AC1,
                  y = ~ AC2,
                  color = ~ Cluster,
                  text = ~ Cluster,
                  hoverinfo = "text",
                  marker = list(size = 10)
                ) %>% layout(
                  xaxis = list(title = input$Axis1_Gene_Clus),
                  yaxis = list(title = input$Axis2_Gene_Clus)
                )
      return(p)
    })

    output$CellSpace3D_Clus <-
      renderPlotly(
        plot_ly(
          X$MCA$cells_principal %>%
            rownames_to_column(var = "Sample") %>% inner_join(X$cluster$labels, by = "Sample"),
          color = ~ Cluster,
          mode = "markers",
          text = ~ paste(Cluster, " ", Sample),
          x = ~ X$MCA$cells_principal[[input$Axis1_3D_Clus]],
          y = ~ X$MCA$cells_principal[[input$Axis2_3D_Clus]],
          z = ~ X$MCA$cells_principal[[input$Axis3_3D_Clus]],
          hoverinfo = "text",
          marker = list(
            opacity = input$Alpha_3D_Clus,
            size = input$Size_3D_Clus,
            line = list(color = "rgba(0, 0, 0, .8)", width = 2)
          )
        ) %>% add_markers() %>%
          layout(
            autosize = FALSE,
            legend = list(x = 100, y = 0.5),
            scene = list(
              xaxis = list(title = input$Axis1_3D_Clus),
              yaxis = list(title = input$Axis2_3D_Clus),
              zaxis = list(title = input$Axis3_3D_Clus)
            )
          )
      )

    output$DTBOX <- DT::renderDataTable({
      DTboxplot <- data_frame()
      Data <-
        X$cluster$gene_cluster_distances %>% gather("Cluster", "Distance", -Genes)
      for (i in (Data$Cluster %>% unique))
      {
        Cluster <- i
        bin1 <-
          Data %>% filter(Cluster == i) %>% arrange_("Distance") %>%
          separate(Genes,
                   into = c("Genes", "bin"),
                   sep = "-bin") %>% filter(bin ==
                                              1) %>% extract2("Genes") %>% head(5) %>% paste(collapse = " ")
        bin2 <-
          Data %>% filter(Cluster == i) %>% arrange_("Distance") %>%
          separate(Genes,
                   into = c("Genes", "bin"),
                   sep = "-bin") %>% filter(bin ==
                                              2) %>% extract2("Genes") %>% head(5) %>% paste(collapse = " ")
        DTboxplot <-
          bind_rows(DTboxplot, data_frame(Cluster, bin1, bin2))
      }
      return(DTboxplot %>% DT::datatable(rownames = FALSE))
    })
    output$Boxplot <-
      (renderPlotly(
        Boxplot %>% filter(Genes %in% input$Genes_Boxplot) %>%
          plot_ly(x = ~ Genes, y = ~ Expression) %>% add_trace(
            type = "box",
            color = ~ Cluster,
            jitter = 0.5,
            pointpos = 0,
            boxpoints = "all",
            boxmean = "sd"
          ) %>% layout(
            boxmode = "group",
            margin = list(
              b = 100,
              t = 25,
              l = 50,
              r = 50,
              pad = 0
            )
          )
      ))

    output$Heatmap_Expression_Clus <-
      renderPlotly(X %>% Heatmap_Cluster(n = input$Num_Top_Genes, plotly = T))
  }
  return(shinyApp(ui, server))
}

#' Create Dashboard
#'
#' @param X MCXpress object
#' @return Shiny object
create_dashboard3 <- function(X)
{
  dr_axis <-
    X$MCA$cells_principal %>% select(contains("Axis")) %>% colnames
  ui <-
    dashboardPage(skin = sample(c("red", "blue", "yellow", "green"), size = 1),
      dashboardHeader(title = h1("MCXpress")),
      dashboardSidebar(sidebarMenu(
        menuItem("MCA",
                 tabName = "mca", icon = icon("arrows")),
        menuItem(
          "Clustering",
          tabName = "clus",
          icon = icon("object-group")
        ),
        menuItem("GSEA", tabName = "gsea", icon = icon("gears"))
      )),
      dashboardBody(CSS, tabItems(
        TabMCA(X, dr_axis), TabClus(X, dr_axis), TabGSEA(X,
                                                         dr_axis)
      ))
    )
  server <- function(input, output, clientData, session)
  {
    MCA_axis_name <-
      X$MCA$cells_principal %>% select(contains("Axis")) %>% colnames

    output$CellSpaceGeneCor <- renderPlotly({
      if (input$MCA_CS_AC_Type == "Principal")
      {
        axis_cor <-
          X$MCA$cells_principal %>% rownames_to_column(var = "Sample") %>%
          inner_join(
            X$ExpressionMatrix[input$MCA_CS_AC_Gene, ] %>% data.frame() %>%
              tibble::rownames_to_column() %>% set_colnames(c("Sample", "Expression")),
            by = "Sample"
          )
        p <-
          plot_ly(data = axis_cor,
                  x = ~ axis_cor[[input$MCA_CS_AC_Axis_x]],
                  y = ~ axis_cor[[input$MCA_CS_AC_Axis_y]]) %>% add_markers(
                    text = ~ Sample,
                    hoverinfo = "text",
                    alpha = input$MCA_CS_AC_Alpha,
                    color = ~ Expression,
                    marker = list(size = input$MCA_CS_AC_Size)
                  ) %>% layout(
                    xaxis = list(title = input$MCA_CS_AC_Axis_x),
                    yaxis = list(title = input$MCA_CS_AC_Axis_y)
                  )
        p
      } else
      {
        axis_cor <-
          X$MCA$cells_standard %>% rownames_to_column(var = "Sample") %>%
          inner_join(
            X$ExpressionMatrix[input$MCA_CS_AC_Gene, ] %>% data.frame() %>%
              tibble::rownames_to_column() %>% set_colnames(c("Sample", "Expression")),
            by = "Sample"
          )
        p <-
          plot_ly(data = axis_cor,
                  x = ~ axis_cor[[input$MCA_CS_AC_Axis_x]],
                  y = ~ axis_cor[[input$MCA_CS_AC_Axis_y]]) %>% add_markers(
                    text = ~ Sample,
                    hoverinfo = "text",
                    alpha = input$MCA_CS_AC_Alpha,
                    color = ~ Expression,
                    marker = list(size = input$MCA_CS_AC_Size)
                  ) %>% layout(
                    xaxis = list(title = input$MCA_CS_Axis_x),
                    yaxis = list(title = input$MCA_CS_Axis_y)
                  )
        p
      }
    })
    output$CellSpace2D <- renderPlotly({
      if (input$MCA_CS_2D_Type == "Principal")
      {
        d3 <- X$MCA$cells_principal %>% rownames_to_column(var = "Sample")
        p <-
          plot_ly(
            data = d3,
            x = ~ d3[[input$MCA_CS_Axis_x]],
            y = ~ d3[[input$MCA_CS_Axis_y]],
            type = "scatter",
            mode = "markers",
            text = ~ Sample,
            hoverinfo = "text",
            alpha = input$Alpha,
            marker = list(size = input$Size)
          ) %>% layout(
            xaxis = list(title = input$MCA_CS_Axis_x),
            yaxis = list(title = input$MCA_CS_Axis_y)
          )
        p
      } else
      {
        d3 <- X$MCA$cells_standard %>% rownames_to_column(var = "Sample")
        p <-
          plot_ly(
            data = d3,
            x = ~ d3[[input$MCA_CS_Axis_x]],
            y = ~ d3[[input$MCA_CS_Axis_y]],
            type = "scatter",
            mode = "markers",
            text = ~ Sample,
            hoverinfo = "text",
            alpha = input$Alpha,
            marker = list(size = input$Size)
          ) %>% layout(
            xaxis = list(title = input$MCA_CS_Axis_x),
            yaxis = list(title = input$MCA_CS_Axis_y)
          )
        p
      }
    })

    output$CellSpace3D <-
      renderPlotly(if (input$MCA_CS_3D_Type == "Standard")
      {
        plot_ly(
          X$MCA$cells_standard,
          mode = "markers",
          text = ~ paste(
            rownames(X$MCA$cells_standard),
            "</br>",
            input$MCA_CS_Axis1_3D,
            ": ",
            X$MCA$cells_standard[[input$MCA_CS_Axis1_3D]] %>%
              signif(digits = 4),
            "</br>",
            input$Axis2_3D,
            ": ",
            X$MCA$cells_standard[[input$MCA_CS_Axis2_3D]] %>%
              signif(digits = 4),
            "</br>",
            input$Axis3_3D,
            ": ",
            X$MCA$cells_standard[[input$MCA_CS_Axis3_3D]] %>%
              signif(digits = 4)
          ),
          x = ~ X$MCA$cells_standard[[input$MCA_CS_Axis1_3D]],
          y = ~ X$MCA$cells_standard[[input$MCA_CS_Axis2_3D]],
          z = ~ X$MCA$cells_standard[[input$MCA_CS_Axis3_3D]],
          hoverinfo = "text",
          marker = list(
            opacity = input$MCA_CS_Alpha_3D,
            size = input$MCA_CS_Size_3D,
            line = list(color = "rgba(0, 0, 0, .8)",
                        width = 2)
          )
        ) %>% add_markers() %>% layout(
          autosize = FALSE,
          scene = list(
            xaxis = list(title = input$MCA_CS_Axis1_3D),
            yaxis = list(title = input$MCA_CS_Axis2_3D),
            zaxis = list(title = input$MCA_CS_Axis3_3D)
          )
        )
      } else
      {
        plot_ly(
          X$MCA$cells_standard,
          mode = "markers",
          text = ~ paste(
            rownames(X$MCA$cells_principal),
            "</br>",
            input$MCA_CS_Axis1_3D,
            ": ",
            X$MCA$cells_principal[[input$MCA_CS_Axis1_3D]] %>%
              signif(digits = 4),
            "</br>",
            input$MCA_CS_Axis2_3D,
            ": ",
            X$MCA$cells_principal[[input$MCA_CS_Axis2_3D]] %>%
              signif(digits = 4),
            "</br>",
            input$MCA_CS_Axis3_3D,
            ": ",
            X$MCA$cells_principal[[input$MCA_CS_Axis3_3D]] %>%
              signif(digits = 4)
          ),
          x = ~ X$MCA$cells_principal[[input$MCA_CS_Axis1_3D]],
          y = ~ X$MCA$cells_principal[[input$MCA_CS_Axis2_3D]],
          z = ~ X$MCA$cells_principal[[input$MCA_CS_Axis3_3D]],
          hoverinfo = "text",
          marker = list(
            opacity = input$MCA_CS_Alpha_3D,
            size = input$MCA_CS_Size_3D,
            line = list(color = "rgba(0, 0, 0, .8)",
                        width = 2)
          )
        ) %>% add_markers() %>% layout(
          margin = list(
            b = 50,
            t = 50,
            l = 50,
            r = 100,
            pad = 0
          ),
          autosize = FALSE,
          scene = list(
            xaxis = list(title = input$MCA_CS_Axis1_3D),
            yaxis = list(title = input$MCA_CS_Axis2_3D),
            zaxis = list(title = input$MCA_CS_Axis3_3D)
          )
        )
      })
    output$TableGeneCor <-
      DT::renderDataTable(X$MCA$Axis_Gene_Cor %>% extract(1:6))
    output$GeneSpace <- renderPlotly({
      if (input$Type_Gene == "Principal")
      {
        d3 <- X$MCA$genes_principal %>% rownames_to_column(var = "Genes")
      } else
      {
        d3 <- X$MCA$genes_standard %>% rownames_to_column(var = "Genes")
      }
      p <-
        plot_ly(
          data = d3,
          x = ~ d3[[input$Axis1_Gene]],
          y = ~ d3[[input$Axis2_Gene]],
          type = "scatter",
          mode = "markers",
          text = ~ Genes,
          hoverinfo = "text",
          alpha = input$Alpha_Gene,
          marker = list(size = input$Size_Gene)
        ) %>%
        layout(
          xaxis = list(title = input$Axis1_Gene),
          yaxis = list(title = input$Axis2_Gene)
        )
      p
    })
    output$Eigen <- renderPlotly({
      Shiny_Eigen <- X$MCA$explained_eigen_variance
      Shiny_Eigen <- Shiny_Eigen[1:input$amount_adjust, ]
      plot_ly(Shiny_Eigen) %>% add_bars(
        x = ~ Axis,
        y = ~ Shiny_Eigen[[input$Mode]],
        hoverinfo = "text",
        text = ~ paste0(
          "Axis: ",
          Axis,
          "</br>Explained Variance: ",
          Explained_Variance %>% round(3),
          "</br>Cumulated Explained Variance: ",
          Cumulative %>% round(3)
        )
      ) %>% layout(
        yaxis = list(title = input$Mode),
        margin = list(
          b = 100,
          t = 25,
          l = 50,
          r = 50,
          pad = 0
        )
      )
    })


    # ____________________________________________________________________________
    # Clus server ####

    options(warn = -1)
    output$CellSpace_Clus <- renderPlotly({
      d3 <-
        X$MCA$cells_principal %>% rownames_to_column(var = "Sample") %>%
        inner_join(X$cluster$labels, by = "Sample")
      d4 <-
        X$MCA$cells_standard %>% rownames_to_column(var = "Sample") %>%
        inner_join(X$cluster$labels, by = "Sample")

      if (input$Type == "Principal")
      {
        p <-
          plot_ly(
            data = d3,
            x = ~ d3[[input$Axis1_Clus]],
            y = ~ d3[[input$Axis2_Clus]],
            color = ~ Cluster,
            type = "scatter",
            mode = "markers",
            text = ~ Sample,
            hoverinfo = "text",
            alpha = input$Alpha_Clus,
            marker = list(size = input$Size_Clus)
          ) %>%
          layout(
            xaxis = list(title = input$Axis1_Clus),
            yaxis = list(title = input$Axis2_Clus)
          )
        p
      } else
      {
        p <-
          plot_ly(
            data = d4,
            x = ~ d4[[input$Axis1_Clus]],
            y = ~ d4[[input$Axis2_Clus]],
            color = ~ Cluster,
            type = "scatter",
            mode = "markers",
            text = ~ Sample,
            hoverinfo = "text",
            alpha = input$Alpha_Clus,
            marker = list(size = input$Size_Clus)
          ) %>%
          layout(
            xaxis = list(title = input$Axis1_Clus),
            yaxis = list(title = input$Axis2_Clus)
          )
        p
      }
    })

    output$GeneSpace_Clus <- renderPlotly({
      Genes <-
        X$MCA$genes_standard %>% rownames_to_column(var = "Genes") %>%
        select_("Genes", input$Axis1_Gene_Clus, input$Axis2_Gene_Clus) %>%
        set_colnames(c("Genes", "AP1", "AP2"))
      Centroids <-
        X$cluster$coord_centroids %>% select_("Cluster",
                                              input$Axis1_Gene_Clus,
                                              input$Axis2_Gene_Clus) %>% set_colnames(c("Cluster", "AC1", "AC2"))
      p <-
        plot_ly(data = Genes,
                x = ~ AP1,
                y = ~ AP2) %>% add_markers(
                  name = "Genes",
                  text = ~ Genes,
                  hoverinfo = "text",
                  marker = list(
                    size = input$Size_Gene_Clus,
                    color = "black",
                    alpha = input$Alpha_Gene_Clus
                  )
                ) %>% add_markers(
                  data = Centroids,
                  x = ~ AC1,
                  y = ~ AC2,
                  color = ~ Cluster,
                  text = ~ Cluster,
                  hoverinfo = "text",
                  marker = list(size = 10)
                ) %>% layout(
                  xaxis = list(title = input$Axis1_Gene_Clus),
                  yaxis = list(title = input$Axis2_Gene_Clus)
                )
      return(p)
    })

    output$CellSpace3D_Clus <-
      renderPlotly(
        plot_ly(
          X$MCA$cells_principal %>%
            rownames_to_column(var = "Sample") %>% inner_join(X$cluster$labels, by = "Sample"),
          color = ~ Cluster,
          mode = "markers",
          text = ~ paste(Cluster, " ", Sample),
          x = ~ X$MCA$cells_principal[[input$Axis1_3D_Clus]],
          y = ~ X$MCA$cells_principal[[input$Axis2_3D_Clus]],
          z = ~ X$MCA$cells_principal[[input$Axis3_3D_Clus]],
          hoverinfo = "text",
          marker = list(
            opacity = input$Alpha_3D_Clus,
            size = input$Size_3D_Clus,
            line = list(color = "rgba(0, 0, 0, .8)", width = 2)
          )
        ) %>% add_markers() %>%
          layout(
            autosize = FALSE,
            legend = list(x = 100, y = 0.5),
            scene = list(
              xaxis = list(title = input$Axis1_3D_Clus),
              yaxis = list(title = input$Axis2_3D_Clus),
              zaxis = list(title = input$Axis3_3D_Clus)
            )
          )
      )

    output$DTBOX <- DT::renderDataTable({
      DTboxplot <- data_frame()
      Data <-
        X$cluster$gene_cluster_distances %>% gather("Cluster", "Distance", -Genes)
      for (i in (Data$Cluster %>% unique))
      {
        Cluster <- i
        bin1 <-
          Data %>% filter(Cluster == i) %>% arrange_("Distance") %>%
          separate(Genes,
                   into = c("Genes", "bin"),
                   sep = "-bin") %>% filter(bin ==
                                              1) %>% extract2("Genes") %>% head(5) %>% paste(collapse = " ")
        bin2 <-
          Data %>% filter(Cluster == i) %>% arrange_("Distance") %>%
          separate(Genes,
                   into = c("Genes", "bin"),
                   sep = "-bin") %>% filter(bin ==
                                              2) %>% extract2("Genes") %>% head(5) %>% paste(collapse = " ")
        DTboxplot <-
          bind_rows(DTboxplot, data_frame(Cluster, bin1, bin2))
      }
      return(DTboxplot %>% DT::datatable(rownames = FALSE))
    })
    Boxplot <-
      X$ExpressionMatrix %>% data.frame %>% rownames_to_column(var = "Genes") %>%
      set_colnames(c("Genes", X$ExpressionMatrix %>% colnames)) %>% gather("Sample",
                                                                           "Expression", -Genes) %>% as_tibble %>% arrange(Sample) %>% inner_join(X$cluster$labels,
                                                                                                                                                  by = "Sample") %>% as_tibble

    output$Heatmap_Expression_Clus <-
      renderPlotly(X %>% Heatmap_Cluster(n = input$Num_Top_Genes, plotly = T))

    output$Boxplot <-
      (renderPlotly(
        Boxplot %>% filter(Genes %in% input$Genes_Boxplot) %>%
          plot_ly(x = ~ Genes, y = ~ Expression) %>% add_trace(
            type = "box",
            color = ~ Cluster,
            jitter = 0.5,
            pointpos = 0,
            boxpoints = "all",
            boxmean = "sd"
          ) %>% layout(
            boxmode = "group",
            margin = list(
              b = 100,
              t = 25,
              l = 50,
              r = 50,
              pad = 0
            )
          )
      ))

    # ____________________________________________________________________________
    # GSEA server ####

    observeEvent(input$Mode_Func_Plot, {
      updateSelectInput(
        session,
        "Choice_Func_Plot",
        label = paste("Choose",
                      input$Mode_Func_Plot),
        choices = {
          switch(
            input$Mode_Func_Plot,
            Cluster = X$GSEA$GSEA_Results %>% names,
            Axis = X$GSEA$GSEA_Results_Axis %>% names
          )
        }
      )
    }, ignoreNULL = TRUE)
    Data <- reactive(X$GSEA$AllRanking[[input$Choice_Func_Plot]])
    output$GSEA <-
      renderPlotly(plotlyEnrichment(X$GSEA$GMTfile[[input$Geneset]],
                                    Data(), gseaParam = X$GSEA$gseaParam))

    # Datatable of Enrichment results
    observeEvent(input$Mode_Func_DT, {
      updateSelectInput(
        session,
        "Choice_Func_DT",
        label = paste("Choose",
                      input$Mode_Func_DT),
        choices = {
          switch(
            input$Mode_Func_DT,
            Cluster = X$GSEA$GSEA_Results %>% names,
            Axis = X$GSEA$GSEA_Results_Axis %>% names
          )
        }
      )
    })
    Table_Enrich1 <- reactive({
      switch(
        input$Mode_Func_DT,
        Cluster = (X$GSEA$GSEA_Results),
        Axis = (X$GSEA$GSEA_Results_Axis)
      )
    })
    Table_Enrich2 <-
      eventReactive(
        input$Choice_Func_DT,
        Table_Enrich1()[[input$Choice_Func_DT]] %>%
          select(-nMoreExtreme, -leadingEdge),
        ignoreNULL = TRUE
      )
    output$Table <-
      DT::renderDataTable(Table_Enrich2() %>% DT::datatable(rownames = FALSE))

    Enrich_Boxplot <-
      X$GSEA$GSEA_Results %>% map2(
        .y = X$GSEA$GSEA_Results %>%
          names,
        .f = function(x, y)
        {
          mutate(.data = x, Cluster = y)
        }
      ) %>% bind_rows

    output$GSEA_Heatmap <- (renderPlotly({
      colvector <- c("Yellow_Red", "Grey", "RG", "Spectral", "Cool_Warm")
      if (input$GSEA_Heatmap_Level == "Cluster")
      {
        func         <-
          list(heatmaply::YlOrRd,
               heatmaply::Greys,
               gplots::redgreen,
               heatmaply::Spectral,
               heatmaply::cool_warm)
        index        <- (colvector == input$GSEA_Heatmap_Color) %>% which
        indexrow     <- (colvector == input$GSEA_Heatmap_Row_Side_Color) %>% which
        colo         <- func[[index]]
        colorow      <- func[[indexrow]]
        GSEA_Heatmap_Cluster(
          X,
          pval    = input$GSEA_Heatmap_Pval,
          es      = input$GSEA_Heatmap_ES,
          nes     = input$GSEA_Heatmap_NES,
          metrics = input$GSEA_Heatmap_Metrics,
          nPath   = input$GSEA_N_Path,
          margin  = c(150, 150),
          plotly  = T,
          cexCol  = input$GSEA_Heatmap_Cex_Col,
          cexRow  = input$GSEA_Heatmap_Cex_Row,
          color   = 50 %>% colo %>%  rev,
          row_color = colorow
        )
      }
      else{
        func      <-
          list(heat.colors,
               cm.colors,
               gplots::redgreen,
               heatmaply::Spectral)
        index        <- (colvector == input$GSEA_Heatmap_Color) %>% which
        indexrow     <- (colvector == input$GSEA_Heatmap_Row_Side_Color) %>% which
        colo         <- func[[index]]
        colorow      <- func[[indexrow]]
        GSEA_Heatmap_SC(
          X,
          pval    = input$GSEA_Heatmap_Pval,
          es      = input$GSEA_Heatmap_ES,
          nes     = input$GSEA_Heatmap_NES,
          metrics = input$GSEA_Heatmap_Metrics,
          nPath   = input$GSEA_N_Path,
          plotly  = T,
          margin  = c(150, 150),
          cexCol  = input$GSEA_Heatmap_Cex_Col,
          cexRow  = input$GSEA_Heatmap_Cex_Row,
          color   = 50 %>% colo %>%  rev,
          row_color = colorow
        )
      }
    }))
  }
  return(shinyApp(ui, server))
}