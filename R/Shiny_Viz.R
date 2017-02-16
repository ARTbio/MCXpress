Create_Shiny_Cluster <- function(X) {
  App <- shinyApp(
    ui =
      navbarPage(theme= shinytheme("cerulean"),
                 "Clustering",
                 tabPanel(
                   "Cell Space",  titlePanel("Cluster in the Cell Space"),
                   wellPanel(fluidRow(
                     column(
                       5,
                       selectInput(
                         "Axis1",
                         label = "Select x Axis",
                         choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                       "Sample") %>% select(contains("Axis")) %>%  colnames,
                         selected = "Axis1"
                       ),
                       selectInput(
                         "Axis2",
                         label = "Select y Axis",
                         choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                       "Sample") %>% select(contains("Axis")) %>%  colnames,
                         selected = "Axis2"
                       ),
                       selectInput(
                         "Type",
                         label = "Type",
                         choices = c("Principal", "Standard"),
                         selected = "Principal"
                       )
                     ),
                     column(
                       5,
                       offset = 1,
                       sliderInput(
                         "Size",
                         label = "Point Size",
                         min = 0,
                         max = 10,
                         value = 5,
                         step = 1
                       ),
                       sliderInput(
                         "Alpha",
                         label = "Transparency",
                         min = 0,
                         max = 1,
                         value = 1,
                         step = 0.1
                       )
                     )
                   )),
                   plotlyOutput("CellSpace", width = "99%", height = "99%")
                 ),
                 tabPanel(
                   title = "Gene Space",  titlePanel("Genespace with Cluster Centroids"),
                   wellPanel(fluidRow(
                     column(
                       5,
                       selectInput(
                         "Axis1_Gene",
                         label = "Select x Axis",
                         choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                       "Sample") %>% select(contains("Axis")) %>%  colnames,
                         selected = "Axis1"
                       ),
                       selectInput(
                         "Axis2_Gene",
                         label = "Select y Axis",
                         choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                       "Sample") %>% select(contains("Axis")) %>%  colnames,
                         selected = "Axis2"
                       )
                     ),
                     column(
                       5,
                       offset = 1,
                       sliderInput(
                         "Size_Gene",
                         label = "Point Size",
                         min = 0,
                         max = 2,
                         value = 1,
                         step = 0.1
                       ),
                       sliderInput(
                         "Alpha_Gene",
                         label = "Transparency",
                         min = 0,
                         max = 1,
                         value = 1,
                         step = 0.1
                       )
                     )
                   )),
                   mainPanel(width=12, plotlyOutput("GeneSpace"))
                 ),
                 tabPanel(
                   "Cell Space 3D",
                   fluidRow(
                     column(
                       5,
                       selectInput(
                         "Axis1_3D",
                         label = "Select x Axis",
                         choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                       "Sample") %>% select(contains("Axis")) %>%  colnames,
                         selected = "Axis1"
                       ),
                       selectInput(
                         "Axis2_3D",
                         label = "Select y Axis",
                         choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                       "Sample") %>% select(contains("Axis")) %>%  colnames,
                         selected = "Axis2"
                       ),
                       selectInput(
                         "Axis3_3D",
                         label = "Select z Axis",
                         choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                       "Sample") %>% select(contains("Axis")) %>%  colnames,
                         selected = "Axis3"
                       )
                     ),
                     column(
                       5,
                       offset = 1,
                       sliderInput(
                         "Size_3D",
                         label = "Point Size",
                         min = 0,
                         max = 10,
                         value = 5,
                         step = 0.1
                       ),
                       sliderInput(
                         "Alpha_3D",
                         label = "Transparency",
                         min = 0,
                         max = 1,
                         value = 1,
                         step = 0.1
                       )
                     )
                   ),mainPanel(width = 12,
                               plotlyOutput("CellSpace3D", width = "100%", height ="100%"
                               ))
                 ),  tabPanel(
                   title = "Boxplot",
                   fluidPage(titlePanel("Genes Expression by Cluster"),
                             fluidRow(column(
                               width = 6,
                               selectInput(
                                 "Genes_Boxplot",
                                 "Choose a Gene:",
                                 choices =
                                   (X$ExpressionMatrix %>%  rownames %>%  sort),
                                 selectize = T, multiple = TRUE, selected = (X$ExpressionMatrix %>%  rownames)[1]
                               )
                             ))),
                   mainPanel(width = 12, h3("Boxplot"), plotlyOutput(
                     "Boxplot", width = "100%", height = "100%"
                   ), h3("Genes Specific to Cluster"), dataTableOutput("DTBOX"))
                 )
      ),
    server = function(input, output){
      options(warn=-1)

      output$DTBOX<-renderDataTable({DTboxplot<-data_frame()
      for(i in (X$cluster$Gene_Cluster_Distance %>% select(-Genes) %>% colnames))
      {
        Cluster<-i
        bin1<-X$cluster$Gene_Cluster_Distance %>%  arrange_(i) %>%  separate(Genes, into=c("Genes", "bin"), sep="-bin") %>%  filter(bin==1) %>% select_("Genes") %>% head(5) %>% as.matrix() %>% as.vector() %>% paste(collapse=" ")
        bin2<-X$cluster$Gene_Cluster_Distance %>%  arrange_(i) %>%  separate(Genes, into=c("Genes", "bin"), sep="-bin") %>%  filter(bin==2) %>% select_("Genes") %>% head(5) %>% as.matrix() %>% as.vector() %>%  paste(collapse=" ")
        DTboxplot<-bind_rows(DTboxplot,data_frame(Cluster,bin1,bin2))
      }
      return(DTboxplot %>% datatable(rownames = FALSE))})

      Boxplot<-X$ExpressionMatrix %>%  data.frame %>%  rownames_to_column(var="Genes") %>%  gather("Sample","Expression",-Genes) %>%  arrange(Sample)  %>%  inner_join(X$cluster$Cluster_Quali, by="Sample")
      output$Boxplot <-
        (renderPlotly(Boxplot %>% filter(Genes %in% input$Genes_Boxplot) %>%  plot_ly(x=~Genes, y=~Expression) %>% add_trace( type="box", split=~Cluster, jitter=0.5, pointpos=0, boxpoints="all") %>%  layout(boxmode="group")))
      output$CellSpace <- renderPlotly({
        if (input$Type == "Principal") {
          d3 <-
            X$Dim_Red$Cells_Principal %>%  rownames_to_column(var = "Sample") %>%  inner_join(X$cluster$Cluster_Quali, by = "Sample")

          p<-plot_ly(data=d3, x=~d3[[input$Axis1]], y=~d3[[input$Axis2]], color =~Cluster, type = "scatter", mode = "markers", text=~Sample, hoverinfo="text",alpha= input$Alpha, marker = list(size = input$Size))%>% layout(xaxis = list(title=input$Axis1), yaxis = list(title=input$Axis2))
          p
        }
        else{
          d3 <-
            X$Dim_Red$Cells_Standard %>%  rownames_to_column(var = "Sample")  %>%  inner_join(X$cluster$Cluster_Quali, by = "Sample")
          p<-plot_ly(data=d3, x=~d3[[input$Axis1]], y=~d3[[input$Axis2]], color =~Cluster, type = "scatter", mode = "markers", text=~Sample, hoverinfo="text",alpha= input$Alpha, marker = list(size = input$Size))%>% layout(xaxis = list(title=input$Axis1), yaxis = list(title=input$Axis2))
          p
        }
      })
      GeneSpace <-
        X$Dim_Red$Genes_Standard %>% rownames_to_column(var = "Genes")
      output$GeneSpace <- renderPlotly({
        p<-plot_ly(data=GeneSpace, x=~GeneSpace[[input$Axis1_Gene]], y=~GeneSpace[[input$Axis2_Gene]], type = "scatter", mode = "markers", text=~Genes, hoverinfo="text",alpha= input$Alpha_Gene, marker = list(size = input$Size_Gene, color='rgba(0,0,0,1)', showlegend=FALSE), name="Genes")%>% layout(xaxis = list(title=input$Axis1_Gene), yaxis = list(title=input$Axis2_Gene)) %>%  add_markers(data=X$cluster$Coord_Centroids,x=~X$cluster$Coord_Centroids[[input$Axis1_Gene]] , y=~X$cluster$Coord_Centroids[[input$Axis2_Gene]],color=~Cluster,colors = "Blues", mode='markers',marker = list(size = 10, color=~Cluster), text=~Cluster)
      })
      output$CellSpace3D <- renderPlotly(
        plot_ly(
          X$Dim_Red$Cells_Standard,
          color =  ~ X$cluster$Cluster_Quali$Cluster,
          mode = 'markers',
          text = ~ paste(
            rownames(X$Dim_Red$Cells_Standard),
            '</br>',
            input$Axis1_3D,
            ': ',
            X$Dim_Red$Cells_Standard[[input$Axis1_3D]] %>%  signif(digits = 4),
            '</br>',
            input$Axis2_3D,
            ': ',
            X$Dim_Red$Cells_Standard[[input$Axis2_3D]] %>%  signif(digits = 4),
            '</br>',
            input$Axis3_3D,
            ': ',
            X$Dim_Red$Cells_Standard[[input$Axis3_3D]] %>%  signif(digits = 4)
          )
          ,
          x = ~ X$Dim_Red$Cells_Standard[[input$Axis1_3D]],
          y = ~ X$Dim_Red$Cells_Standard[[input$Axis2_3D]],
          z = ~ X$Dim_Red$Cells_Standard[[input$Axis3_3D]],
          hoverinfo = "text",
          marker = list(
            opacity= input$Alpha_3D,
            size = input$Size_3D,
            line = list(color = 'rgba(0, 0, 0, .8)',
                        width = 2)
          )
        ) %>% add_markers() %>%
          layout(
            autosize = F,
            legend = list(x = 100, y = 0.5),
            scene = list(
              xaxis = list(title = input$Axis1_3D),
              yaxis = list(title = input$Axis2_3D),
              zaxis = list(title = input$Axis3_3D)
            )
          )
      )
    }
  )
  return(App)
}

Create_Shiny_Dim_Red <- function(X) {
  App <- shinyApp(
    ui = navbarPage(theme= shinytheme("cerulean"),
                    "Dimension Reduction",
                    tabPanel(
                      "Cell Space",
                      fluidRow(
                        column(
                          4, wellPanel(
                            selectInput(
                              "nAxis_Dim_Red",
                              label = "Representation",
                              choices = c("2 Axis", "3 Axis")
                            ),
                            selectInput(
                              "Type",
                              label = "Type",
                              choices = c("Principal", "Standard"),
                              selected = "Principal"
                            )
                          )),
                        column(
                          4, uiOutput("ui_Axis_Dim_Red")),
                        column(
                          4, uiOutput("ui_Option_Dim_Red"))
                      ),
                      mainPanel(width = 12,
                                plotlyOutput("CellSpace", width = "95%", height =
                                               "100%"))
                    ),
                    tabPanel(
                      title = "Gene Space",
                      fluidRow(
                        column(
                          5, wellPanel(
                            selectInput(
                              "Axis1_Gene",
                              label = "Select x Axis1",
                              choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                            "Sample") %>% select(contains("Axis")) %>%  colnames,
                              selected = "Axis1"
                            ),

                            selectInput(
                              "Axis2_Gene",
                              label = "Select y Axis",
                              choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                            "Sample") %>% select(contains("Axis")) %>%  colnames,
                              selected = "Axis2"
                            ),
                            selectInput(
                              "Type_Gene",
                              label = "Type",
                              choices = c("Principal", "Standard"),
                              selected = "Principal"
                            )
                          )
                        ),
                        column(
                          5, offset = 2,wellPanel(

                            sliderInput(
                              "Size_Gene",
                              label = "Point Size",
                              min = 0,
                              max = 5,
                              value = 2,
                              step = 1
                            ),

                            sliderInput(
                              "Alpha_Gene",
                              label = "Transparency",
                              min = 0,
                              max = 1,
                              value = 1,
                              step = 0.1
                            )
                          )
                        )
                      ),
                      mainPanel(width = 12,
                                plotlyOutput("GeneSpace"))
                    ),
                    tabPanel(
                      title = "Eigen Value",
                      selectInput(
                        "Mode",
                        label = "Select Mode",
                        choices = c("Cumul" = "Cumul", "EigenValue" =
                                      "EigenValue"),
                        selected = "EigenValue"
                      ),
                      sliderInput(
                        "amount_adjust",
                        label = "Amount",
                        min = 1,
                        max = X$Dim_Red$Cells_Principal %>% ncol ,
                        value = if(X$Dim_Red$Cells_Principal %>% ncol %>% is_greater_than(5)){5} else{1},
                        step = 1
                      ),
                      plotlyOutput("Eigen")
                    )
    ),
    server = function(input, output) {
      output$ui_Axis_Dim_Red <- renderUI(
        {
          switch(input$nAxis_Dim_Red,
                 "2 Axis" = wellPanel(
                   selectInput(
                     "Axis1",
                     label = "Select x Axis",
                     choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                   "Sample") %>% select(contains("Axis")) %>%  colnames,
                     selected = "Axis1"
                   ),

                   selectInput(
                     "Axis2",
                     label = "Select y Axis",
                     choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                   "Sample") %>% select(contains("Axis")) %>%  colnames,
                     selected = "Axis2"
                   )
                 ),
                 "3 Axis" =  wellPanel(
                   selectInput(
                     "Axis1_3D",
                     label = "Select x Axis",
                     choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                   "Sample") %>% select(contains("Axis")) %>%  colnames,
                     selected = "Axis1"
                   ),
                   selectInput(
                     "Axis2_3D",
                     label = "Select y Axis",
                     choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                   "Sample") %>% select(contains("Axis")) %>%  colnames,
                     selected = "Axis2"
                   ),
                   selectInput(
                     "Axis3_3D",
                     label = "Select z Axis",
                     choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
                                                                                   "Sample") %>% select(contains("Axis")) %>%  colnames,
                     selected = "Axis3"
                   )
                 )
          )
        }
      )

      output$ui_Option_Dim_Red <- renderUI({
        switch(input$nAxis_Dim_Red,
               "2 Axis" = wellPanel(
                 sliderInput(
                   "Size",
                   label = "Point Size",
                   min = 0,
                   max = 10,
                   value = 5,
                   step = 1
                 ),
                 sliderInput(
                   "Alpha",
                   label = "Transparency",
                   min = 0,
                   max = 1,
                   value = 1,
                   step = 0.1
                 )
               ),
               "3 Axis" =  wellPanel(
                 sliderInput(
                   "Size_3D",
                   label = "Point Size",
                   min = 0,
                   max = 10,
                   value = 5,
                   step = 1
                 ),
                 sliderInput(
                   "Alpha_3D",
                   label = "Transparency",
                   min = 0,
                   max = 1,
                   value = 1,
                   step = 0.1
                 )
               )
        )
      })

      output$CellSpace <- renderPlotly({
        switch(input$nAxis_Dim_Red,
               "2 Axis"={
                 if (input$Type == "Principal"){
                   d3 <- X$Dim_Red$Cells_Principal %>%  rownames_to_column(var = "Sample")
                   p<-plot_ly(data=d3, x=~d3[[input$Axis1]], y=~d3[[input$Axis2]], type = "scatter", mode = "markers", text=~Sample, hoverinfo="text",alpha= input$Alpha, marker = list(size = input$Size))%>% layout(xaxis = list(title=input$Axis1), yaxis = list(title=input$Axis2))
                   p
                 }
                 else{
                   d3 <- X$Dim_Red$Cells_Standard %>%  rownames_to_column(var = "Sample")
                   p<-plot_ly(data=d3, x=~d3[[input$Axis1]], y=~d3[[input$Axis2]], type = "scatter", mode = "markers", text=~Sample, hoverinfo="text", alpha= input$Alpha, marker = list(size = input$Size))%>% layout(xaxis = list(title=input$Axis1), yaxis = list(title=input$Axis2))
                   p
                 }},
               "3 Axis"=
                 ( if(input$Type == "Standard") {
                   plot_ly(
                     X$Dim_Red$Cells_Standard,
                     mode = 'markers',
                     text = ~ paste(
                       rownames(X$Dim_Red$Cells_Standard),
                       '</br>',
                       input$Axis1_3D,
                       ': ',
                       X$Dim_Red$Cells_Standard[[input$Axis1_3D]] %>%  signif(digits = 4),
                       '</br>',
                       input$Axis2_3D,
                       ': ',
                       X$Dim_Red$Cells_Standard[[input$Axis2_3D]] %>%  signif(digits = 4),
                       '</br>',
                       input$Axis3_3D,
                       ': ',
                       X$Dim_Red$Cells_Standard[[input$Axis3_3D]] %>%  signif(digits = 4)
                     )
                     ,
                     x = ~ X$Dim_Red$Cells_Standard[[input$Axis1_3D]],
                     y = ~ X$Dim_Red$Cells_Standard[[input$Axis2_3D]],
                     z = ~ X$Dim_Red$Cells_Standard[[input$Axis3_3D]],
                     hoverinfo = "text",
                     marker = list(
                       opacity=input$Alpha_3D,
                       size = input$Size_3D,
                       line = list(color = 'rgba(0, 0, 0, .8)',
                                   width = 2)
                     )
                   ) %>% add_markers() %>%
                     layout(
                       autosize = F,
                       scene = list(
                         xaxis = list(title = input$Axis1_3D),
                         yaxis = list(title = input$Axis2_3D),
                         zaxis = list(title = input$Axis3_3D)
                       )
                     )} else{plot_ly(
                       X$Dim_Red$Cells_Standard,
                       mode = 'markers',
                       text = ~ paste(
                         rownames(X$Dim_Red$Cells_Principal),
                         '</br>',
                         input$Axis1_3D,
                         ': ',
                         X$Dim_Red$Cells_Principal[[input$Axis1_3D]] %>%  signif(digits = 4),
                         '</br>',
                         input$Axis2_3D,
                         ': ',
                         X$Dim_Red$Cells_Principal[[input$Axis2_3D]] %>%  signif(digits = 4),
                         '</br>',
                         input$Axis3_3D,
                         ': ',
                         X$Dim_Red$Cells_Principal[[input$Axis3_3D]] %>%  signif(digits = 4)
                       )
                       ,
                       x = ~ X$Dim_Red$Cells_Principal[[input$Axis1_3D]],
                       y = ~ X$Dim_Red$Cells_Principal[[input$Axis2_3D]],
                       z = ~ X$Dim_Red$Cells_Principal[[input$Axis3_3D]],
                       width = 800,
                       height = 600,
                       hoverinfo = "text",
                       marker = list(
                         opacity=input$Alpha_3D,
                         size = input$Size_3D,
                         line = list(color = 'rgba(0, 0, 0, .8)',
                                     width = 2)
                       )
                     ) %>% add_markers() %>%
                         layout(
                           autosize = F,
                           scene = list(
                             xaxis = list(title = input$Axis1_3D),
                             yaxis = list(title = input$Axis2_3D),
                             zaxis = list(title = input$Axis3_3D)
                           )
                         )}
                 )
        )
      }
      )


      output$GeneSpace <- renderPlotly({
        if (input$Type_Gene == "Principal") {d3<-X$Dim_Red$Genes_Principal %>% rownames_to_column(var = "Genes")} else{d3<-X$Dim_Red$Genes_Standard %>% rownames_to_column(var = "Genes")}
        p<-plot_ly(data=d3, x=~d3[[input$Axis1_Gene]], y=~d3[[input$Axis2_Gene]], type = "scatter", mode = "markers", text=~Genes, hoverinfo="text", alpha= input$Alpha_Gene, marker = list(size = input$Size_Gene))%>% layout(xaxis = list(title=input$Axis1_Gene), yaxis = list(title=input$Axis2_Gene))
        p
      })
      output$Eigen <- renderPlotly({
        d4 <- X$Dim_Red$Cumul[X$Dim_Red$Cumul$Type == input$Mode, ]
        d4 <- d4[1:input$amount_adjust, ]
        (
          d4 %>% ggplot(aes(x = Axis, y = Value)) + geom_bar(stat = "identity") + scale_x_discrete(limits =
                                                                                                   d4$Axis)
        ) %>% ggplotly
      })
    }
  )
  return(App)
}


Create_Shiny_Functionnal_Analysis <- function(X) {
  App <- shinyApp(
    ui =
      navbarPage(theme= shinytheme("cerulean"),
                 "Functionnal Analysis",
                 tabPanel(title = "Enrichmentplot",
                          fluidPage(
                            titlePanel("Enrichmentplot"),
                            wellPanel(
                              fluidRow(column(
                                width = 6,
                                selectInput(
                                  "Geneset",
                                  "Choose a Geneset:",
                                  choices = (X$Functionnal_Analysis$GSEA_Results$pathway %>%  unique),
                                  selectize = T
                                )
                              ),
                              column(
                                6,
                                selectInput(
                                  "Choice_Func_Plot",
                                  "Choose a Cluster:",
                                  choices = X$Functionnal_Analysis$GSEA_Results$Cluster %>%  unique
                                ),
                                selectInput(
                                  "Mode_Func_Plot",
                                  "Choose",
                                  choices = c("Cluster", "Axis")
                                )
                              ))),
                            mainPanel(width=12,plotlyOutput(
                              "GSEA", width = "100%", height = "100%"
                            ))
                          )),
                 tabPanel(
                   title = "Enrichment Results",
                   fluidPage(titlePanel("Enrichment Table"),
                             wellPanel(fluidRow(column(
                               width = 6,
                               selectInput(
                                 "Mode_Func_DT",
                                 "Choose Data:",
                                 choices =
                                   (c("Cluster","Axis")),
                                 selectize = T
                               )
                             ),column(
                               width = 6,
                               selectInput(
                                 "Choice_Func_DT",
                                 "Choose a Geneset:",
                                 choices =
                                   (X$Functionnal_Analysis$GSEA_Results$Cluster %>%  unique),
                                 selectize = T
                               )
                             )))),
                   mainPanel(width=12,dataTableOutput(
                     "Table", width = "100%", height = "100%"
                   ))
                 ),
                 tabPanel(
                   title = "Boxplot Enrichment Cluster",
                   fluidPage(titlePanel("Genes Expression by Cluster"),
                             fluidRow(column(
                               width = 10,
                               selectInput(width = "100%",
                                 "Enrich_Heatmap_GeneSet",
                                 "Choose a Gene:",
                                 choices =
                                   (X$Functionnal_Analysis$GSEA_Results$pathway %>%  unique),
                                 selectize = T, multiple = TRUE, selected = (X$Functionnal_Analysis$GSEA_Results$pathway %>%  unique)[1]
                               )
                             ))),
                   mainPanel(width = 12, h3("Enrichment Score Heatmap"), plotlyOutput(
                     "Enrich_Heatmap", width = "90%", height = "90%"
                   ))
                 )
      )
    ,
    server = function(input, output, clientData, session) {
      observeEvent(input$Mode_Func_Plot,{
        updateSelectInput(session,"Choice_Func_Plot",
                          label= paste("Choose", input$Mode_Func_Plot),
                          choices={switch(input$Mode_Func_Plot,
                                          "Cluster"=X$Functionnal_Analysis$GSEA_Results$Cluster %>%  unique,
                                          "Axis"=X$Functionnal_Analysis$GSEA_Results_Axis$Axis %>%  unique
                          )})},ignoreNULL = TRUE)





      Data<-reactive(X$Functionnal_Analysis$Grouped %>% filter(Group==input$Choice_Func_Plot))
      RankinG <- reactive(Data() %>% select(Ranking) %>%  as.matrix %>%as.numeric %>% as.vector() %>% setNames(Data()$Genes %>% unique))
      output$GSEA<-renderPlotly(plotlyEnrichment(GMTfile[[input$Geneset]], RankinG(), gseaParam = 1))

      #Datatable of Enrichment results
      observeEvent(input$Mode_Func_DT,{
        updateSelectInput(session,"Choice_Func_DT",
                          label= paste("Choose", input$Mode_Func_DT),
                          choices={switch(input$Mode_Func_DT,
                                          "Cluster"=X$Functionnal_Analysis$GSEA_Results$Cluster %>%  unique,
                                          "Axis"=X$Functionnal_Analysis$GSEA_Results_Axis$Axis %>%  unique
                          )}
        )
      })
      Table_Enrich<-reactive({switch(input$Mode_Func_DT, "Cluster"=(X$Functionnal_Analysis$GSEA_Results %>% filter(Cluster==input$Choice_Func_DT) %>% select(-Cluster, -nMoreExtreme, -leadingEdge)), "Axis"=(X$Functionnal_Analysis$GSEA_Results_Axis %>% filter(Axis==input$Choice_Func_DT) %>% select(-Axis, -nMoreExtreme, -leadingEdge)))})
      output$Table <- renderDataTable(Table_Enrich()%>%  datatable(rownames=FALSE))

      Enrich_Boxplot<-X$Functionnal_Analysis$GSEA_Results
      output$Enrich_Heatmap <-
        (renderPlotly(Enrich_Boxplot %>% filter(pathway %in% input$Enrich_Heatmap_GeneSet) %>%   group_by(pathway) %>%  plot_ly(x=~pathway, y=~Cluster, z=~(ES %>% as.numeric)) %>% add_heatmap()
))




      }
  )
  return(App)
}
