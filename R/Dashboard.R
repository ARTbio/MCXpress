Create_Dashboard1 <- function(X) {
  dr_axis<-X$Dim_Red$Cells_Principal %>% select(contains("Axis")) %>%  colnames
  ui <- dashboardPage(
    dashboardHeader(title=h1("MCXpress")),
    dashboardSidebar(
      sidebarMenu(
        menuItem("MCA", tabName = "mca", icon = icon("arrows"))
        )
      ),
    dashboardBody(includeCSS("C:/Users/Akira/Documents/MCXpress/R/www/custom.css"),
#   ____________________________________________________________________________
#   MCA UI                                                                  ####
tabItems(
  tabItem(tabName = "mca", navbarPage(
    "MCA",
    navbarMenu("Cell Space",
      tabPanel(
        "2 Axis",
        wellPanel(fluidRow(
          column(
           6,selectInput(
             "DR_CS_2D_Type",
             label = "Type",
             choices = c("Principal", "Standard"),
             selected = "Principal"
             ),
           selectInput(
             "DR_CS_Axis_x",
             label = "Select x Axis",
             choices = dr_axis,
             selected = "Axis1"
             ),
           selectInput(
             "DR_CS_Axis_y",
             label = "Select y Axis",
             choices = dr_axis,
             selected = "Axis2"
             )
           ),
          column(
           6,
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
mainPanel(width = 12,
  plotlyOutput("CellSpace2D", width = "95%", height =
   "100%"))
), tabPanel(
"3 Axis",wellPanel(
 fluidRow(
   column(4,
     selectInput(
       "DR_CS_3D_Type",
       label = "Type",
       choices = c("Principal", "Standard"),
       selected = "Principal"
       ), selectInput(
       "DR_CS_Axis1_3D",
       label = "Select x Axis",
       choices = dr_axis,
       selected = "Axis1"
       )
       ),
   column(
     4,
     selectInput(
       "DR_CS_Axis2_3D",
       label = "Select y Axis",
       choices = dr_axis,
       selected = "Axis2"
       ),
     selectInput(
       "DR_CS_Axis3_3D",
       label = "Select z Axis",
       choices = dr_axis,
       selected = "Axis3"
       )
     ),
   column(
     4,
     sliderInput(
       "DR_CS_Size_3D",
       label = "Point Size",
       min = 0,
       max = 10,
       value = 5,
       step = 0.1
       ),
     sliderInput(
       "DR_CS_Alpha_3D",
       label = "Transparency",
       min = 0,
       max = 1,
       value = 1,
       step = 0.1
       )
     )
   )),mainPanel(width = 12,align="center",
plotlyOutput("CellSpace3D", width = "100%", height ="100%"
 ))
   ), tabPanel("Axis & Gene Correlation", wellPanel(fluidRow(
    column(
     4,selectInput(
       "DR_CS_AC_Type",
       label = "Type",
       choices = c("Principal", "Standard"),
       selected = "Principal"
       ),
     selectInput(
       "DR_CS_AC_Gene",
       label = "Gene",
       choices = X$ExpressionMatrix %>%  rownames
       )
     ),
    column(4,
      selectInput(
       "DR_CS_AC_Axis_x",
       label = "Select x Axis",
       choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
         "Sample") %>% select(contains("Axis")) %>%  colnames,
       selected = "Axis1"
       ),
      selectInput(
       "DR_CS_AC_Axis_y",
       label = "Select y Axis",
       choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
         "Sample") %>% select(contains("Axis")) %>%  colnames,
       selected = "Axis2"
       )
      ),
    column(
     4,
     sliderInput(
       "DR_CS_AC_Size",
       label = "Point Size",
       min = 0,
       max = 10,
       value = 5,
       step = 1
       ),
     sliderInput(
       "DR_CS_AC_Alpha",
       label = "Transparency",
       min = 0,
       max = 1,
       value = 1,
       step = 0.1
       )
     )
    )),mainPanel(width = 12,
plotlyOutput("CellSpaceGeneCor", width = "100%", height ="100%"
 )), mainPanel(width = 12,
    dataTableOutput("TableGeneCor", width = "100%", height ="100%"
     ))
 )
),
tabPanel(
  title = "Gene Space",
  fluidRow(
    column(
      5, wellPanel(
        selectInput(
          "Axis1_Gene",
          label = "Select x Axis1",
          choices = dr_axis,
          selected = "Axis1"
          ), selectInput(
          "Axis2_Gene",
          label = "Select y Axis",
          choices = dr_axis,
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
    choices = c("Explained_Variance" = "Explained_Variance", "Cumulative" =
      "Cumulative"),
    selected = "Explained_Variance"
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
)
)
)
)
)

server <- function(input, output,clientData, session) {
  DR_axis_name <- X$Dim_Red$Cells_Principal %>% select(contains("Axis")) %>%  colnames

  output$CellSpaceGeneCor<- renderPlotly({
   if (input$DR_CS_AC_Type == "Principal"){
     axis_cor <- X$Dim_Red$Cells_Principal %>%  rownames_to_column(var = "Sample") %>%  inner_join(X$ExpressionMatrix[input$DR_CS_AC_Gene,] %>% data.frame() %>% tibble::rownames_to_column() %>%  set_colnames(c("Sample", "Expression")), by="Sample")
     p<-plot_ly(data=axis_cor, x=~axis_cor[[input$DR_CS_AC_Axis_x]], y=~axis_cor[[input$DR_CS_AC_Axis_y]]) %>% add_markers(text=~Sample, hoverinfo="text",alpha= input$DR_CS_AC_Alpha, color=~Expression, marker = list(size = input$DR_CS_AC_Size)) %>% layout(xaxis = list(title=input$DR_CS_AC_Axis_x), yaxis = list(title=input$DR_CS_AC_Axis_y))
     p
   }
   else{
     axis_cor <- X$Dim_Red$Cells_Standard %>%  rownames_to_column(var = "Sample")%>%  inner_join(X$ExpressionMatrix[input$DR_CS_AC_Gene,] %>% data.frame() %>% tibble::rownames_to_column() %>%  set_colnames(c("Sample", "Expression")), by="Sample")
     p<-plot_ly(data=axis_cor, x=~axis_cor[[input$DR_CS_AC_Axis_x]], y=~axis_cor[[input$DR_CS_AC_Axis_y]]) %>% add_markers(text=~Sample, hoverinfo="text",alpha= input$DR_CS_AC_Alpha, color=~Expression, marker = list(size = input$DR_CS_AC_Size)) %>% layout(xaxis = list(title=input$DR_CS_Axis_x), yaxis = list(title=input$DR_CS_Axis_y))
     p
     }})
output$CellSpace2D <- renderPlotly({
 if (input$DR_CS_2D_Type == "Principal"){
   d3 <- X$Dim_Red$Cells_Principal %>%  rownames_to_column(var = "Sample")
   p<-plot_ly(data=d3, x=~d3[[input$DR_CS_Axis_x]], y=~d3[[input$DR_CS_Axis_y]], type = "scatter", mode = "markers", text=~Sample, hoverinfo="text",alpha= input$Alpha, marker = list(size = input$Size))%>% layout(xaxis = list(title=input$DR_CS_Axis_x), yaxis = list(title=input$DR_CS_Axis_y))
   p
 }
 else{
   d3 <- X$Dim_Red$Cells_Standard %>%  rownames_to_column(var = "Sample")
   p<-plot_ly(data=d3, x=~d3[[input$DR_CS_Axis_x]], y=~d3[[input$DR_CS_Axis_y]], type = "scatter", mode = "markers", text=~Sample, hoverinfo="text", alpha= input$Alpha, marker = list(size = input$Size))%>% layout(xaxis = list(title=input$DR_CS_Axis_x), yaxis = list(title=input$DR_CS_Axis_y))
   p
   }})

output$CellSpace3D<- renderPlotly(if(input$DR_CS_3D_Type == "Standard"){
 plot_ly(
   X$Dim_Red$Cells_Standard,
   mode = 'markers',
   text = ~ paste(
     rownames(X$Dim_Red$Cells_Standard),
     '</br>',
     input$DR_CS_Axis1_3D,
     ': ',
     X$Dim_Red$Cells_Standard[[input$DR_CS_Axis1_3D]] %>%  signif(digits = 4),
     '</br>',
     input$Axis2_3D,
     ': ',
     X$Dim_Red$Cells_Standard[[input$DR_CS_Axis2_3D]] %>%  signif(digits = 4),
     '</br>',
     input$Axis3_3D,
     ': ',
     X$Dim_Red$Cells_Standard[[input$DR_CS_Axis3_3D]] %>%  signif(digits = 4)
     )
   ,
   x = ~ X$Dim_Red$Cells_Standard[[input$DR_CS_Axis1_3D]],
   y = ~ X$Dim_Red$Cells_Standard[[input$DR_CS_Axis2_3D]],
   z = ~ X$Dim_Red$Cells_Standard[[input$DR_CS_Axis3_3D]],
   hoverinfo = "text",
   marker = list(
     opacity=input$DR_CS_Alpha_3D,
     size = input$DR_CS_Size_3D,
     line = list(color = 'rgba(0, 0, 0, .8)',
       width = 2)
     )
   ) %>% add_markers() %>%
 layout(
   autosize = F,
   scene = list(
     xaxis = list(title = input$DR_CS_Axis1_3D),
     yaxis = list(title = input$DR_CS_Axis2_3D),
     zaxis = list(title = input$DR_CS_Axis3_3D)
     )
   )} else{plot_ly(
     X$Dim_Red$Cells_Standard,
     mode = 'markers',
     text = ~ paste(
       rownames(X$Dim_Red$Cells_Principal),
       '</br>',
       input$DR_CS_Axis1_3D,
       ': ',
       X$Dim_Red$Cells_Principal[[input$DR_CS_Axis1_3D]] %>%  signif(digits = 4),
       '</br>',
       input$DR_CS_Axis2_3D,
       ': ',
       X$Dim_Red$Cells_Principal[[input$DR_CS_Axis2_3D]] %>%  signif(digits = 4),
       '</br>',
       input$DR_CS_Axis3_3D,
       ': ',
       X$Dim_Red$Cells_Principal[[input$DR_CS_Axis3_3D]] %>%  signif(digits = 4)
       )
     ,
     x = ~ X$Dim_Red$Cells_Principal[[input$DR_CS_Axis1_3D]],
     y = ~ X$Dim_Red$Cells_Principal[[input$DR_CS_Axis2_3D]],
     z = ~ X$Dim_Red$Cells_Principal[[input$DR_CS_Axis3_3D]],
     hoverinfo = "text",
     marker = list(
       opacity=input$DR_CS_Alpha_3D,
       size = input$DR_CS_Size_3D,
       line = list(color = 'rgba(0, 0, 0, .8)',
         width = 2)
       )
     ) %>% add_markers() %>%
   layout(margin=list(b=50, t=50, l=50, r=100, pad=0),
     autosize = F,
     scene = list(
       xaxis = list(title = input$DR_CS_Axis1_3D),
       yaxis = list(title = input$DR_CS_Axis2_3D),
       zaxis = list(title = input$DR_CS_Axis3_3D)
       ))})
output$TableGeneCor<- renderDataTable(X$Dim_Red$Axis_Gene_Cor %>% extract(1:6))
output$GeneSpace <- renderPlotly({
  if (input$Type_Gene == "Principal") {d3<-X$Dim_Red$Genes_Principal %>% rownames_to_column(var = "Genes")} else{d3<-X$Dim_Red$Genes_Standard %>% rownames_to_column(var = "Genes")}
  p<-plot_ly(data=d3, x=~d3[[input$Axis1_Gene]], y=~d3[[input$Axis2_Gene]], type = "scatter", mode = "markers", text=~Genes, hoverinfo="text", alpha= input$Alpha_Gene, marker = list(size = input$Size_Gene))%>% layout(xaxis = list(title=input$Axis1_Gene), yaxis = list(title=input$Axis2_Gene))
  p
  })
output$Eigen <- renderPlotly({
  Shiny_Eigen <- X$Dim_Red$Explained_Eigen_Variance
  Shiny_Eigen <- Shiny_Eigen[1:input$amount_adjust, ]
  plot_ly(Shiny_Eigen) %>% add_bars(x =~Axis, y=~Shiny_Eigen[[input$Mode]], hoverinfo="text", text=~paste0("Axis: ",Axis,'</br>Explained Variance: ', Explained_Variance %>% round(3),'</br>Cumulated Explained Variance: ',Cumulative %>% round(3))) %>% layout(yaxis = list(title=input$Mode), margin=list(b=100, t=25, l=50, r=50, pad=0))
  })}
return(shinyApp(ui, server))
}

Create_Dashboard2 <- function(X) {
  dr_axis<-X$Dim_Red$Cells_Principal %>% select(contains("Axis")) %>%  colnames
  ui <- dashboardPage(
    dashboardHeader(title=h1("MCXpress")),
    dashboardSidebar(
      sidebarMenu(
        menuItem("MCA", tabName = "mca", icon = icon("arrows")),
        menuItem("Clustering", tabName = "clus", icon = icon("object-group"))
        )
      ),
    dashboardBody(includeCSS("C:/Users/Akira/Documents/MCXpress/R/www/custom.css"),
#   ____________________________________________________________________________
#   MCA UI                                                                  ####
tabItems(
  tabItem(tabName = "mca", navbarPage(
    "MCA",
    navbarMenu("Cell Space",
      tabPanel(
        "2 Axis",
        wellPanel(fluidRow(
          column(
           6,selectInput(
             "DR_CS_2D_Type",
             label = "Type",
             choices = c("Principal", "Standard"),
             selected = "Principal"
             ),
           selectInput(
             "DR_CS_Axis_x",
             label = "Select x Axis",
             choices = dr_axis,
             selected = "Axis1"
             ),
           selectInput(
             "DR_CS_Axis_y",
             label = "Select y Axis",
             choices = dr_axis,
             selected = "Axis2"
             )
           ),
          column(
           6,
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
mainPanel(width = 12,
  plotlyOutput("CellSpace2D", width = "95%", height =
   "100%"))
), tabPanel(
"3 Axis",wellPanel(
 fluidRow(
   column(4,
     selectInput(
       "DR_CS_3D_Type",
       label = "Type",
       choices = c("Principal", "Standard"),
       selected = "Principal"
       ), selectInput(
       "DR_CS_Axis1_3D",
       label = "Select x Axis",
       choices = dr_axis,
       selected = "Axis1"
       )
       ),
   column(
     4,
     selectInput(
       "DR_CS_Axis2_3D",
       label = "Select y Axis",
       choices = dr_axis,
       selected = "Axis2"
       ),
     selectInput(
       "DR_CS_Axis3_3D",
       label = "Select z Axis",
       choices = dr_axis,
       selected = "Axis3"
       )
     ),
   column(
     4,
     sliderInput(
       "DR_CS_Size_3D",
       label = "Point Size",
       min = 0,
       max = 10,
       value = 5,
       step = 0.1
       ),
     sliderInput(
       "DR_CS_Alpha_3D",
       label = "Transparency",
       min = 0,
       max = 1,
       value = 1,
       step = 0.1
       )
     )
   )),mainPanel(width = 12,align="center",
plotlyOutput("CellSpace3D", width = "100%", height ="100%"
 ))
   ), tabPanel("Axis & Gene Correlation", wellPanel(fluidRow(
    column(
     4,selectInput(
       "DR_CS_AC_Type",
       label = "Type",
       choices = c("Principal", "Standard"),
       selected = "Principal"
       ),
     selectInput(
       "DR_CS_AC_Gene",
       label = "Gene",
       choices = X$ExpressionMatrix %>%  rownames
       )
     ),
    column(4,
      selectInput(
       "DR_CS_AC_Axis_x",
       label = "Select x Axis",
       choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
         "Sample") %>% select(contains("Axis")) %>%  colnames,
       selected = "Axis1"
       ),
      selectInput(
       "DR_CS_AC_Axis_y",
       label = "Select y Axis",
       choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
         "Sample") %>% select(contains("Axis")) %>%  colnames,
       selected = "Axis2"
       )
      ),
    column(
     4,
     sliderInput(
       "DR_CS_AC_Size",
       label = "Point Size",
       min = 0,
       max = 10,
       value = 5,
       step = 1
       ),
     sliderInput(
       "DR_CS_AC_Alpha",
       label = "Transparency",
       min = 0,
       max = 1,
       value = 1,
       step = 0.1
       )
     )
    )),mainPanel(width = 12,
plotlyOutput("CellSpaceGeneCor", width = "100%", height ="100%"
 )), mainPanel(width = 12,
    dataTableOutput("TableGeneCor", width = "100%", height ="100%"
     ))
 )
),
tabPanel(
  title = "Gene Space",
  fluidRow(
    column(
      5, wellPanel(
        selectInput(
          "Axis1_Gene",
          label = "Select x Axis1",
          choices = dr_axis,
          selected = "Axis1"
          ), selectInput(
          "Axis2_Gene",
          label = "Select y Axis",
          choices = dr_axis,
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
    choices = c("Explained_Variance" = "Explained_Variance", "Cumulative" =
      "Cumulative"),
    selected = "Explained_Variance"
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
)
),

#   ____________________________________________________________________________
#   Clustering UI                                                           ####

tabItem(tabName = "clus",
 navbarPage("Clustering", navbarMenu("Cell Space",
   tabPanel(
     "Clustering 2 Axis",  titlePanel("Cluster in the Cell Space"),
     wellPanel(fluidRow(
       column(
         5,
         selectInput(
           "Axis1_Clus",
           label = "Select x Axis",
           choices = dr_axis,
           selected = "Axis1"
           ),
         selectInput(
           "Axis2_Clus",
           label = "Select y Axis",
           choices = dr_axis,
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
           "Size_Clus",
           label = "Point Size",
           min = 0,
           max = 10,
           value = 5,
           step = 1
           ),
         sliderInput(
           "Alpha_Clus",
           label = "Transparency",
           min = 0,
           max = 1,
           value = 1,
           step = 0.1
           )
         )
       )),
     plotlyOutput("CellSpace_Clus", width = "100%", height = "100%")
     ), tabPanel(
     "Clustering 3 Axis",
     wellPanel(fluidRow(
       column(
         5,
         selectInput(
           "Axis1_3D_Clus",
           label = "Select x Axis",
           choices = dr_axis,
           selected = "Axis1"
           ),
         selectInput(
           "Axis2_3D_Clus",
           label = "Select y Axis",
           choices = dr_axis,
           selected = "Axis2"
           ),
         selectInput(
           "Axis3_3D_Clus",
           label = "Select z Axis",
           choices = dr_axis,
           selected = "Axis3"
           )
         ),
       column(
         5,
         offset = 1,
         sliderInput(
           "Size_3D_Clus",
           label = "Point Size",
           min = 0,
           max = 10,
           value = 5,
           step = 0.1
           ),
         sliderInput(
           "Alpha_3D_Clus",
           label = "Transparency",
           min = 0,
           max = 1,
           value = 1,
           step = 0.1
           )
         )
       )),mainPanel(width = 12, align="center",
     plotlyOutput("CellSpace3D_Clus", width = "100%", height ="100%"
       ))
       )),
tabPanel(
 title = "Gene Space",  titlePanel("Genespace with Cluster Centroids"),
 wellPanel(fluidRow(
   column(
     5,
     selectInput(
       "Axis1_Gene_Clus",
       label = "Select x Axis",
       choices = dr_axis,
       selected = "Axis1"
       ),
     selectInput(
       "Axis2_Gene_Clus",
       label = "Select y Axis",
       choices = dr_axis,
       selected = "Axis2"
       )
     ),
   column(
     5,
     offset = 1,
     sliderInput(
       "Size_Gene_Clus",
       label = "Point Size",
       min = 0,
       max = 2,
       value = 1,
       step = 0.1
       ),
     sliderInput(
       "Alpha_Gene_Clus",
       label = "Transparency",
       min = 0,
       max = 1,
       value = 1,
       step = 0.1
       )
     )
   )),
 mainPanel(width=12, plotlyOutput("GeneSpace_Clus"))
 ),
tabPanel(
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
)
)
)
)
)

server <- function(input, output,clientData, session) {
  DR_axis_name <- X$Dim_Red$Cells_Principal %>% select(contains("Axis")) %>%  colnames

  output$CellSpaceGeneCor<- renderPlotly({
   if (input$DR_CS_AC_Type == "Principal"){
     axis_cor <- X$Dim_Red$Cells_Principal %>%  rownames_to_column(var = "Sample") %>%  inner_join(X$ExpressionMatrix[input$DR_CS_AC_Gene,] %>% data.frame() %>% tibble::rownames_to_column() %>%  set_colnames(c("Sample", "Expression")), by="Sample")
     p<-plot_ly(data=axis_cor, x=~axis_cor[[input$DR_CS_AC_Axis_x]], y=~axis_cor[[input$DR_CS_AC_Axis_y]]) %>% add_markers(text=~Sample, hoverinfo="text",alpha= input$DR_CS_AC_Alpha, color=~Expression, marker = list(size = input$DR_CS_AC_Size)) %>% layout(xaxis = list(title=input$DR_CS_AC_Axis_x), yaxis = list(title=input$DR_CS_AC_Axis_y))
     p
   }
   else{
     axis_cor <- X$Dim_Red$Cells_Standard %>%  rownames_to_column(var = "Sample")%>%  inner_join(X$ExpressionMatrix[input$DR_CS_AC_Gene,] %>% data.frame() %>% tibble::rownames_to_column() %>%  set_colnames(c("Sample", "Expression")), by="Sample")
     p<-plot_ly(data=axis_cor, x=~axis_cor[[input$DR_CS_AC_Axis_x]], y=~axis_cor[[input$DR_CS_AC_Axis_y]]) %>% add_markers(text=~Sample, hoverinfo="text",alpha= input$DR_CS_AC_Alpha, color=~Expression, marker = list(size = input$DR_CS_AC_Size)) %>% layout(xaxis = list(title=input$DR_CS_Axis_x), yaxis = list(title=input$DR_CS_Axis_y))
     p
     }})
output$CellSpace2D <- renderPlotly({
 if (input$DR_CS_2D_Type == "Principal"){
   d3 <- X$Dim_Red$Cells_Principal %>%  rownames_to_column(var = "Sample")
   p<-plot_ly(data=d3, x=~d3[[input$DR_CS_Axis_x]], y=~d3[[input$DR_CS_Axis_y]], type = "scatter", mode = "markers", text=~Sample, hoverinfo="text",alpha= input$Alpha, marker = list(size = input$Size))%>% layout(xaxis = list(title=input$DR_CS_Axis_x), yaxis = list(title=input$DR_CS_Axis_y))
   p
 }
 else{
   d3 <- X$Dim_Red$Cells_Standard %>%  rownames_to_column(var = "Sample")
   p<-plot_ly(data=d3, x=~d3[[input$DR_CS_Axis_x]], y=~d3[[input$DR_CS_Axis_y]], type = "scatter", mode = "markers", text=~Sample, hoverinfo="text", alpha= input$Alpha, marker = list(size = input$Size))%>% layout(xaxis = list(title=input$DR_CS_Axis_x), yaxis = list(title=input$DR_CS_Axis_y))
   p
   }})

output$CellSpace3D<- renderPlotly(if(input$DR_CS_3D_Type == "Standard"){
 plot_ly(
   X$Dim_Red$Cells_Standard,
   mode = 'markers',
   text = ~ paste(
     rownames(X$Dim_Red$Cells_Standard),
     '</br>',
     input$DR_CS_Axis1_3D,
     ': ',
     X$Dim_Red$Cells_Standard[[input$DR_CS_Axis1_3D]] %>%  signif(digits = 4),
     '</br>',
     input$Axis2_3D,
     ': ',
     X$Dim_Red$Cells_Standard[[input$DR_CS_Axis2_3D]] %>%  signif(digits = 4),
     '</br>',
     input$Axis3_3D,
     ': ',
     X$Dim_Red$Cells_Standard[[input$DR_CS_Axis3_3D]] %>%  signif(digits = 4)
     )
   ,
   x = ~ X$Dim_Red$Cells_Standard[[input$DR_CS_Axis1_3D]],
   y = ~ X$Dim_Red$Cells_Standard[[input$DR_CS_Axis2_3D]],
   z = ~ X$Dim_Red$Cells_Standard[[input$DR_CS_Axis3_3D]],
   hoverinfo = "text",
   marker = list(
     opacity=input$DR_CS_Alpha_3D,
     size = input$DR_CS_Size_3D,
     line = list(color = 'rgba(0, 0, 0, .8)',
       width = 2)
     )
   ) %>% add_markers() %>%
 layout(
   autosize = F,
   scene = list(
     xaxis = list(title = input$DR_CS_Axis1_3D),
     yaxis = list(title = input$DR_CS_Axis2_3D),
     zaxis = list(title = input$DR_CS_Axis3_3D)
     )
   )} else{plot_ly(
     X$Dim_Red$Cells_Standard,
     mode = 'markers',
     text = ~ paste(
       rownames(X$Dim_Red$Cells_Principal),
       '</br>',
       input$DR_CS_Axis1_3D,
       ': ',
       X$Dim_Red$Cells_Principal[[input$DR_CS_Axis1_3D]] %>%  signif(digits = 4),
       '</br>',
       input$DR_CS_Axis2_3D,
       ': ',
       X$Dim_Red$Cells_Principal[[input$DR_CS_Axis2_3D]] %>%  signif(digits = 4),
       '</br>',
       input$DR_CS_Axis3_3D,
       ': ',
       X$Dim_Red$Cells_Principal[[input$DR_CS_Axis3_3D]] %>%  signif(digits = 4)
       )
     ,
     x = ~ X$Dim_Red$Cells_Principal[[input$DR_CS_Axis1_3D]],
     y = ~ X$Dim_Red$Cells_Principal[[input$DR_CS_Axis2_3D]],
     z = ~ X$Dim_Red$Cells_Principal[[input$DR_CS_Axis3_3D]],
     hoverinfo = "text",
     marker = list(
       opacity=input$DR_CS_Alpha_3D,
       size = input$DR_CS_Size_3D,
       line = list(color = 'rgba(0, 0, 0, .8)',
         width = 2)
       )
     ) %>% add_markers() %>%
   layout(margin=list(b=50, t=50, l=50, r=100, pad=0),
     autosize = F,
     scene = list(
       xaxis = list(title = input$DR_CS_Axis1_3D),
       yaxis = list(title = input$DR_CS_Axis2_3D),
       zaxis = list(title = input$DR_CS_Axis3_3D)
       ))})
output$TableGeneCor<- renderDataTable(X$Dim_Red$Axis_Gene_Cor %>% extract(1:6))
output$GeneSpace <- renderPlotly({
  if (input$Type_Gene == "Principal") {d3<-X$Dim_Red$Genes_Principal %>% rownames_to_column(var = "Genes")} else{d3<-X$Dim_Red$Genes_Standard %>% rownames_to_column(var = "Genes")}
  p<-plot_ly(data=d3, x=~d3[[input$Axis1_Gene]], y=~d3[[input$Axis2_Gene]], type = "scatter", mode = "markers", text=~Genes, hoverinfo="text", alpha= input$Alpha_Gene, marker = list(size = input$Size_Gene))%>% layout(xaxis = list(title=input$Axis1_Gene), yaxis = list(title=input$Axis2_Gene))
  p
  })
output$Eigen <- renderPlotly({
  Shiny_Eigen <- X$Dim_Red$Explained_Eigen_Variance
  Shiny_Eigen <- Shiny_Eigen[1:input$amount_adjust, ]
  plot_ly(Shiny_Eigen) %>% add_bars(x =~Axis, y=~Shiny_Eigen[[input$Mode]], hoverinfo="text", text=~paste0("Axis: ",Axis,'</br>Explained Variance: ', Explained_Variance %>% round(3),'</br>Cumulated Explained Variance: ',Cumulative %>% round(3))) %>% layout(yaxis = list(title=input$Mode), margin=list(b=100, t=25, l=50, r=50, pad=0))
  })


#   ____________________________________________________________________________
#   Clus server                                                             ####

options(warn=-1)
output$CellSpace_Clus <- renderPlotly({
  d3 <-
  X$Dim_Red$Cells_Principal %>%  rownames_to_column(var = "Sample") %>%  inner_join(X$cluster$Cluster_Quali, by = "Sample")
  d4 <-
  X$Dim_Red$Cells_Standard %>%  rownames_to_column(var = "Sample")  %>%  inner_join(X$cluster$Cluster_Quali, by = "Sample")

  if (input$Type == "Principal") {
    p<-plot_ly(data=d3, x=~d3[[input$Axis1_Clus]], y=~d3[[input$Axis2_Clus]], color =~Cluster, type = "scatter", mode = "markers", text=~Sample, hoverinfo="text",alpha= input$Alpha_Clus, marker = list(size = input$Size_Clus))%>% layout(xaxis = list(title=input$Axis1_Clus), yaxis = list(title=input$Axis2_Clus))
    p
  }
  else{
    p<-plot_ly(data=d4, x=~d4[[input$Axis1_Clus]], y=~d4[[input$Axis2_Clus]], color =~Cluster, type = "scatter", mode = "markers", text=~Sample, hoverinfo="text",alpha= input$Alpha_Clus, marker = list(size = input$Size_Clus))%>% layout(xaxis = list(title=input$Axis1_Clus), yaxis = list(title=input$Axis2_Clus))
    p
  }
  })

output$GeneSpace_Clus <- renderPlotly({
  Genes <- X$Dim_Red$Genes_Standard %>% rownames_to_column(var = "Genes") %>%  select_("Genes",input$Axis1_Gene_Clus, input$Axis2_Gene_Clus) %>% set_colnames(c("Genes", "AP1","AP2"))
  Centroids <- X$cluster$Coord_Centroids %>% select_("Cluster",input$Axis1_Gene_Clus, input$Axis2_Gene_Clus) %>%  set_colnames(c("Cluster", "AC1","AC2"))
  p<-plot_ly(data=Genes, x=~AP1, y=~AP2) %>%
  add_markers(name="Genes", text=~Genes, hoverinfo="text", marker=list(size=input$Size_Gene_Clus, color= "black", alpha= input$Alpha_Gene_Clus)) %>%
  add_markers(data=Centroids, x=~AC1, y=~AC2, color=~Cluster, text=~Cluster, hoverinfo="text", marker=list(size=10)) %>%
  layout(xaxis = list(title=input$Axis1_Gene_Clus), yaxis = list(title=input$Axis2_Gene_Clus))
  return(p)
  })

output$CellSpace3D_Clus <- renderPlotly(
  plot_ly(
    X$Dim_Red$Cells_Principal %>% rownames_to_column(var="Sample") %>%  inner_join(X$cluster$Cluster_Quali, by="Sample"),
    color = ~Cluster,
    mode = 'markers',
    text = ~paste(Cluster," ", Sample)
    ,
    x = ~ X$Dim_Red$Cells_Principal[[input$Axis1_3D_Clus]],
    y = ~ X$Dim_Red$Cells_Principal[[input$Axis2_3D_Clus]],
    z = ~ X$Dim_Red$Cells_Principal[[input$Axis3_3D_Clus]],
    hoverinfo = "text",
    marker = list(
      opacity= input$Alpha_3D_Clus,
      size = input$Size_3D_Clus,
      line = list(color = 'rgba(0, 0, 0, .8)',
        width = 2)
      )
    ) %>% add_markers() %>%
  layout(
    autosize = F,
    legend = list(x = 100, y = 0.5),
    scene = list(
      xaxis = list(title = input$Axis1_3D_Clus),
      yaxis = list(title = input$Axis2_3D_Clus),
      zaxis = list(title = input$Axis3_3D_Clus)
      )
    )
  )

output$DTBOX<-renderDataTable({
  DTboxplot<-data_frame()
  Data<-X$cluster$Gene_Cluster_Distance %>%  gather("Cluster", "Distance", -Genes)
  for(i in (Data$Cluster %>% unique))
  {
    Cluster<-i
    bin1<-Data %>% filter(Cluster==i) %>%  arrange_("Distance") %>%  separate(Genes, into=c("Genes", "bin"), sep="-bin")  %>%  filter(bin==1)  %>%  extract2("Genes")  %>% head(5)  %>% paste(collapse=" ")
    bin2<-Data %>% filter(Cluster==i) %>%  arrange_("Distance") %>%  separate(Genes, into=c("Genes", "bin"), sep="-bin")  %>%  filter(bin==2)  %>%  extract2("Genes")  %>% head(5)  %>% paste(collapse=" ")
    DTboxplot<-bind_rows(DTboxplot,data_frame(Cluster,bin1,bin2))
  }
  return(DTboxplot %>% datatable(rownames = FALSE))})

Boxplot<-X$ExpressionMatrix %>%  data.frame %>%  rownames_to_column(var="Genes") %>% set_colnames(c("Genes",X$ExpressionMatrix %>%  colnames)) %>%  gather("Sample","Expression",-Genes)%>% as_tibble %>%  arrange(Sample)  %>%  inner_join(X$cluster$Cluster_Quali, by="Sample") %>% as_tibble
output$Boxplot <-
(renderPlotly(Boxplot %>% filter(Genes %in% input$Genes_Boxplot) %>%  plot_ly(x=~Genes, y=~Expression) %>% add_trace(type="box", color=~Cluster, jitter=0.5, pointpos=0, boxpoints="all", boxmean="sd") %>%  layout(boxmode="group",margin=list(b=100, t=25, l=50, r=50, pad=0) )))
}
return(shinyApp(ui, server))
}


Create_Dashboard3 <- function(X) {
  dr_axis<-X$Dim_Red$Cells_Principal %>% select(contains("Axis")) %>%  colnames
  ui <- dashboardPage(
    dashboardHeader(title=h1("MCXpress"),downloadButton(icon=icon("file-pdf-o"),label = "Download", outputId=)),
    dashboardSidebar(
      sidebarMenu(
        menuItem("MCA", tabName = "mca", icon = icon("arrows")),
        menuItem("Clustering", tabName = "clus", icon = icon("object-group")),
        menuItem("GSEA", tabName = "gsea", icon = icon("gears"))
        )
      ),
    dashboardBody(includeCSS("C:/Users/Akira/Documents/MCXpress/R/www/custom.css"),
#   ____________________________________________________________________________
#   MCA UI                                                                  ####
tabItems(
  tabItem(tabName = "mca", navbarPage(
    "MCA",
    navbarMenu("Cell Space",
      tabPanel(
        "2 Axis",
        wellPanel(fluidRow(
          column(
           6,selectInput(
             "DR_CS_2D_Type",
             label = "Type",
             choices = c("Principal", "Standard"),
             selected = "Principal"
             ),
           selectInput(
             "DR_CS_Axis_x",
             label = "Select x Axis",
             choices = dr_axis,
             selected = "Axis1"
             ),
           selectInput(
             "DR_CS_Axis_y",
             label = "Select y Axis",
             choices = dr_axis,
             selected = "Axis2"
             )
           ),
          column(
           6,
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
mainPanel(width = 12,
  plotlyOutput("CellSpace2D", width = "95%", height =
   "100%"))
), tabPanel(
"3 Axis",wellPanel(
 fluidRow(
   column(4,
     selectInput(
       "DR_CS_3D_Type",
       label = "Type",
       choices = c("Principal", "Standard"),
       selected = "Principal"
       ), selectInput(
       "DR_CS_Axis1_3D",
       label = "Select x Axis",
       choices = dr_axis,
       selected = "Axis1"
       )
       ),
   column(
     4,
     selectInput(
       "DR_CS_Axis2_3D",
       label = "Select y Axis",
       choices = dr_axis,
       selected = "Axis2"
       ),
     selectInput(
       "DR_CS_Axis3_3D",
       label = "Select z Axis",
       choices = dr_axis,
       selected = "Axis3"
       )
     ),
   column(
     4,
     sliderInput(
       "DR_CS_Size_3D",
       label = "Point Size",
       min = 0,
       max = 10,
       value = 5,
       step = 0.1
       ),
     sliderInput(
       "DR_CS_Alpha_3D",
       label = "Transparency",
       min = 0,
       max = 1,
       value = 1,
       step = 0.1
       )
     )
   )),mainPanel(width = 12,align="center",
plotlyOutput("CellSpace3D", width = "100%", height ="100%"
 ))
   ), tabPanel("Axis & Gene Correlation", wellPanel(fluidRow(
    column(
     4,selectInput(
       "DR_CS_AC_Type",
       label = "Type",
       choices = c("Principal", "Standard"),
       selected = "Principal"
       ),
     selectInput(
       "DR_CS_AC_Gene",
       label = "Gene",
       choices = X$ExpressionMatrix %>%  rownames
       )
     ),
    column(4,
      selectInput(
       "DR_CS_AC_Axis_x",
       label = "Select x Axis",
       choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
         "Sample") %>% select(contains("Axis")) %>%  colnames,
       selected = "Axis1"
       ),
      selectInput(
       "DR_CS_AC_Axis_y",
       label = "Select y Axis",
       choices = X$Dim_Red$Cells_Principal %>%  rownames_to_column(var =
         "Sample") %>% select(contains("Axis")) %>%  colnames,
       selected = "Axis2"
       )
      ),
    column(
     4,
     sliderInput(
       "DR_CS_AC_Size",
       label = "Point Size",
       min = 0,
       max = 10,
       value = 5,
       step = 1
       ),
     sliderInput(
       "DR_CS_AC_Alpha",
       label = "Transparency",
       min = 0,
       max = 1,
       value = 1,
       step = 0.1
       )
     )
    )),mainPanel(width = 12,
plotlyOutput("CellSpaceGeneCor", width = "100%", height ="100%"
 )), mainPanel(width = 12,
    dataTableOutput("TableGeneCor", width = "100%", height ="100%"
     ))
 )
),
tabPanel(
  title = "Gene Space",
  fluidRow(
    column(
      5, wellPanel(
        selectInput(
          "Axis1_Gene",
          label = "Select x Axis1",
          choices = dr_axis,
          selected = "Axis1"
          ), selectInput(
          "Axis2_Gene",
          label = "Select y Axis",
          choices = dr_axis,
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
    choices = c("Explained_Variance" = "Explained_Variance", "Cumulative" =
      "Cumulative"),
    selected = "Explained_Variance"
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
)
),

#   ____________________________________________________________________________
#   Clustering UI                                                           ####

tabItem(tabName = "clus",
 navbarPage("Clustering", navbarMenu("Cell Space",
   tabPanel(
     "Clustering 2 Axis",  titlePanel("Cluster in the Cell Space"),
     wellPanel(fluidRow(
       column(
         5,
         selectInput(
           "Axis1_Clus",
           label = "Select x Axis",
           choices = dr_axis,
           selected = "Axis1"
           ),
         selectInput(
           "Axis2_Clus",
           label = "Select y Axis",
           choices = dr_axis,
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
           "Size_Clus",
           label = "Point Size",
           min = 0,
           max = 10,
           value = 5,
           step = 1
           ),
         sliderInput(
           "Alpha_Clus",
           label = "Transparency",
           min = 0,
           max = 1,
           value = 1,
           step = 0.1
           )
         )
       )),
     plotlyOutput("CellSpace_Clus", width = "100%", height = "100%")
     ), tabPanel(
     "Clustering 3 Axis",
     wellPanel(fluidRow(
       column(
         5,
         selectInput(
           "Axis1_3D_Clus",
           label = "Select x Axis",
           choices = dr_axis,
           selected = "Axis1"
           ),
         selectInput(
           "Axis2_3D_Clus",
           label = "Select y Axis",
           choices = dr_axis,
           selected = "Axis2"
           ),
         selectInput(
           "Axis3_3D_Clus",
           label = "Select z Axis",
           choices = dr_axis,
           selected = "Axis3"
           )
         ),
       column(
         5,
         offset = 1,
         sliderInput(
           "Size_3D_Clus",
           label = "Point Size",
           min = 0,
           max = 10,
           value = 5,
           step = 0.1
           ),
         sliderInput(
           "Alpha_3D_Clus",
           label = "Transparency",
           min = 0,
           max = 1,
           value = 1,
           step = 0.1
           )
         )
       )),mainPanel(width = 12, align="center",
     plotlyOutput("CellSpace3D_Clus", width = "100%", height ="100%"
       ))
       )),
tabPanel(
 title = "Gene Space",  titlePanel("Genespace with Cluster Centroids"),
 wellPanel(fluidRow(
   column(
     5,
     selectInput(
       "Axis1_Gene_Clus",
       label = "Select x Axis",
       choices = dr_axis,
       selected = "Axis1"
       ),
     selectInput(
       "Axis2_Gene_Clus",
       label = "Select y Axis",
       choices = dr_axis,
       selected = "Axis2"
       )
     ),
   column(
     5,
     offset = 1,
     sliderInput(
       "Size_Gene_Clus",
       label = "Point Size",
       min = 0,
       max = 2,
       value = 1,
       step = 0.1
       ),
     sliderInput(
       "Alpha_Gene_Clus",
       label = "Transparency",
       min = 0,
       max = 1,
       value = 1,
       step = 0.1
       )
     )
   )),
 mainPanel(width=12, plotlyOutput("GeneSpace_Clus"))
 ),
tabPanel(
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
)
),


#   ____________________________________________________________________________
#   GSEA UI                                                                 ####



tabItem(tabName = "gsea",
  navbarPage(
   "GSEA",
   tabPanel(title = "Enrichmentplot",
    fluidPage(
      titlePanel("Enrichmentplot"),
      wellPanel(
        fluidRow(column(
          width = 6,
          selectInput(
            "Geneset",
            "Choose a Geneset:",
            choices = (X$Functionnal_Analysis$Pathways),
            selectize = T
            )
          ),
        column(
          6,
          selectInput(
            "Choice_Func_Plot",
            "Choose a Cluster:",
            choices = X$Functionnal_Analysis$GSEA_Results %>% names
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
           (X$Functionnal_Analysis$GSEA_Results %>%  names),
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
           (X$Functionnal_Analysis$Pathways),
           selectize = T, multiple = TRUE, selected = (X$Functionnal_Analysis$Pathways)[1]
           )
         ))),
     mainPanel(width = 12, h3("Enrichment Score Heatmap"), plotlyOutput(
       "Enrich_Heatmap", width = "90%", height = "90%"
       ))
     )
   )
)
)


)

)

server <- function(input, output,clientData, session) {
  DR_axis_name <- X$Dim_Red$Cells_Principal %>% select(contains("Axis")) %>%  colnames

  output$CellSpaceGeneCor<- renderPlotly({
   if (input$DR_CS_AC_Type == "Principal"){
     axis_cor <- X$Dim_Red$Cells_Principal %>%  rownames_to_column(var = "Sample") %>%  inner_join(X$ExpressionMatrix[input$DR_CS_AC_Gene,] %>% data.frame() %>% tibble::rownames_to_column() %>%  set_colnames(c("Sample", "Expression")), by="Sample")
     p<-plot_ly(data=axis_cor, x=~axis_cor[[input$DR_CS_AC_Axis_x]], y=~axis_cor[[input$DR_CS_AC_Axis_y]]) %>% add_markers(text=~Sample, hoverinfo="text",alpha= input$DR_CS_AC_Alpha, color=~Expression, marker = list(size = input$DR_CS_AC_Size)) %>% layout(xaxis = list(title=input$DR_CS_AC_Axis_x), yaxis = list(title=input$DR_CS_AC_Axis_y))
     p
   }
   else{
     axis_cor <- X$Dim_Red$Cells_Standard %>%  rownames_to_column(var = "Sample")%>%  inner_join(X$ExpressionMatrix[input$DR_CS_AC_Gene,] %>% data.frame() %>% tibble::rownames_to_column() %>%  set_colnames(c("Sample", "Expression")), by="Sample")
     p<-plot_ly(data=axis_cor, x=~axis_cor[[input$DR_CS_AC_Axis_x]], y=~axis_cor[[input$DR_CS_AC_Axis_y]]) %>% add_markers(text=~Sample, hoverinfo="text",alpha= input$DR_CS_AC_Alpha, color=~Expression, marker = list(size = input$DR_CS_AC_Size)) %>% layout(xaxis = list(title=input$DR_CS_Axis_x), yaxis = list(title=input$DR_CS_Axis_y))
     p
     }})
output$CellSpace2D <- renderPlotly({
 if (input$DR_CS_2D_Type == "Principal"){
   d3 <- X$Dim_Red$Cells_Principal %>%  rownames_to_column(var = "Sample")
   p<-plot_ly(data=d3, x=~d3[[input$DR_CS_Axis_x]], y=~d3[[input$DR_CS_Axis_y]], type = "scatter", mode = "markers", text=~Sample, hoverinfo="text",alpha= input$Alpha, marker = list(size = input$Size))%>% layout(xaxis = list(title=input$DR_CS_Axis_x), yaxis = list(title=input$DR_CS_Axis_y))
   p
 }
 else{
   d3 <- X$Dim_Red$Cells_Standard %>%  rownames_to_column(var = "Sample")
   p<-plot_ly(data=d3, x=~d3[[input$DR_CS_Axis_x]], y=~d3[[input$DR_CS_Axis_y]], type = "scatter", mode = "markers", text=~Sample, hoverinfo="text", alpha= input$Alpha, marker = list(size = input$Size))%>% layout(xaxis = list(title=input$DR_CS_Axis_x), yaxis = list(title=input$DR_CS_Axis_y))
   p
   }})

output$CellSpace3D<- renderPlotly(if(input$DR_CS_3D_Type == "Standard"){
 plot_ly(
   X$Dim_Red$Cells_Standard,
   mode = 'markers',
   text = ~ paste(
     rownames(X$Dim_Red$Cells_Standard),
     '</br>',
     input$DR_CS_Axis1_3D,
     ': ',
     X$Dim_Red$Cells_Standard[[input$DR_CS_Axis1_3D]] %>%  signif(digits = 4),
     '</br>',
     input$Axis2_3D,
     ': ',
     X$Dim_Red$Cells_Standard[[input$DR_CS_Axis2_3D]] %>%  signif(digits = 4),
     '</br>',
     input$Axis3_3D,
     ': ',
     X$Dim_Red$Cells_Standard[[input$DR_CS_Axis3_3D]] %>%  signif(digits = 4)
     )
   ,
   x = ~ X$Dim_Red$Cells_Standard[[input$DR_CS_Axis1_3D]],
   y = ~ X$Dim_Red$Cells_Standard[[input$DR_CS_Axis2_3D]],
   z = ~ X$Dim_Red$Cells_Standard[[input$DR_CS_Axis3_3D]],
   hoverinfo = "text",
   marker = list(
     opacity=input$DR_CS_Alpha_3D,
     size = input$DR_CS_Size_3D,
     line = list(color = 'rgba(0, 0, 0, .8)',
       width = 2)
     )
   ) %>% add_markers() %>%
 layout(
   autosize = F,
   scene = list(
     xaxis = list(title = input$DR_CS_Axis1_3D),
     yaxis = list(title = input$DR_CS_Axis2_3D),
     zaxis = list(title = input$DR_CS_Axis3_3D)
     )
   )} else{plot_ly(
     X$Dim_Red$Cells_Standard,
     mode = 'markers',
     text = ~ paste(
       rownames(X$Dim_Red$Cells_Principal),
       '</br>',
       input$DR_CS_Axis1_3D,
       ': ',
       X$Dim_Red$Cells_Principal[[input$DR_CS_Axis1_3D]] %>%  signif(digits = 4),
       '</br>',
       input$DR_CS_Axis2_3D,
       ': ',
       X$Dim_Red$Cells_Principal[[input$DR_CS_Axis2_3D]] %>%  signif(digits = 4),
       '</br>',
       input$DR_CS_Axis3_3D,
       ': ',
       X$Dim_Red$Cells_Principal[[input$DR_CS_Axis3_3D]] %>%  signif(digits = 4)
       )
     ,
     x = ~ X$Dim_Red$Cells_Principal[[input$DR_CS_Axis1_3D]],
     y = ~ X$Dim_Red$Cells_Principal[[input$DR_CS_Axis2_3D]],
     z = ~ X$Dim_Red$Cells_Principal[[input$DR_CS_Axis3_3D]],
     hoverinfo = "text",
     marker = list(
       opacity=input$DR_CS_Alpha_3D,
       size = input$DR_CS_Size_3D,
       line = list(color = 'rgba(0, 0, 0, .8)',
         width = 2)
       )
     ) %>% add_markers() %>%
   layout(margin=list(b=50, t=50, l=50, r=100, pad=0),
     autosize = F,
     scene = list(
       xaxis = list(title = input$DR_CS_Axis1_3D),
       yaxis = list(title = input$DR_CS_Axis2_3D),
       zaxis = list(title = input$DR_CS_Axis3_3D)
       ))})
output$TableGeneCor<- renderDataTable(X$Dim_Red$Axis_Gene_Cor %>% extract(1:6))
output$GeneSpace <- renderPlotly({
  if (input$Type_Gene == "Principal") {d3<-X$Dim_Red$Genes_Principal %>% rownames_to_column(var = "Genes")} else{d3<-X$Dim_Red$Genes_Standard %>% rownames_to_column(var = "Genes")}
  p<-plot_ly(data=d3, x=~d3[[input$Axis1_Gene]], y=~d3[[input$Axis2_Gene]], type = "scatter", mode = "markers", text=~Genes, hoverinfo="text", alpha= input$Alpha_Gene, marker = list(size = input$Size_Gene))%>% layout(xaxis = list(title=input$Axis1_Gene), yaxis = list(title=input$Axis2_Gene))
  p
  })
output$Eigen <- renderPlotly({
  Shiny_Eigen <- X$Dim_Red$Explained_Eigen_Variance
  Shiny_Eigen <- Shiny_Eigen[1:input$amount_adjust, ]
  plot_ly(Shiny_Eigen) %>% add_bars(x =~Axis, y=~Shiny_Eigen[[input$Mode]], hoverinfo="text", text=~paste0("Axis: ",Axis,'</br>Explained Variance: ', Explained_Variance %>% round(3),'</br>Cumulated Explained Variance: ',Cumulative %>% round(3))) %>% layout(yaxis = list(title=input$Mode), margin=list(b=100, t=25, l=50, r=50, pad=0))
  })


#   ____________________________________________________________________________
#   Clus server                                                             ####

options(warn=-1)
output$CellSpace_Clus <- renderPlotly({
  d3 <-
  X$Dim_Red$Cells_Principal %>%  rownames_to_column(var = "Sample") %>%  inner_join(X$cluster$Cluster_Quali, by = "Sample")
  d4 <-
  X$Dim_Red$Cells_Standard %>%  rownames_to_column(var = "Sample")  %>%  inner_join(X$cluster$Cluster_Quali, by = "Sample")

  if (input$Type == "Principal") {
    p<-plot_ly(data=d3, x=~d3[[input$Axis1_Clus]], y=~d3[[input$Axis2_Clus]], color =~Cluster, type = "scatter", mode = "markers", text=~Sample, hoverinfo="text",alpha= input$Alpha_Clus, marker = list(size = input$Size_Clus))%>% layout(xaxis = list(title=input$Axis1_Clus), yaxis = list(title=input$Axis2_Clus))
    p
  }
  else{
    p<-plot_ly(data=d4, x=~d4[[input$Axis1_Clus]], y=~d4[[input$Axis2_Clus]], color =~Cluster, type = "scatter", mode = "markers", text=~Sample, hoverinfo="text",alpha= input$Alpha_Clus, marker = list(size = input$Size_Clus))%>% layout(xaxis = list(title=input$Axis1_Clus), yaxis = list(title=input$Axis2_Clus))
    p
  }
  })

output$GeneSpace_Clus <- renderPlotly({
  Genes <- X$Dim_Red$Genes_Standard %>% rownames_to_column(var = "Genes") %>%  select_("Genes",input$Axis1_Gene_Clus, input$Axis2_Gene_Clus) %>% set_colnames(c("Genes", "AP1","AP2"))
  Centroids <- X$cluster$Coord_Centroids %>% select_("Cluster",input$Axis1_Gene_Clus, input$Axis2_Gene_Clus) %>%  set_colnames(c("Cluster", "AC1","AC2"))
  p<-plot_ly(data=Genes, x=~AP1, y=~AP2) %>%
  add_markers(name="Genes", text=~Genes, hoverinfo="text", marker=list(size=input$Size_Gene_Clus, color= "black", alpha= input$Alpha_Gene_Clus)) %>%
  add_markers(data=Centroids, x=~AC1, y=~AC2, color=~Cluster, text=~Cluster, hoverinfo="text", marker=list(size=10)) %>%
  layout(xaxis = list(title=input$Axis1_Gene_Clus), yaxis = list(title=input$Axis2_Gene_Clus))
  return(p)
  })

output$CellSpace3D_Clus <- renderPlotly(
  plot_ly(
    X$Dim_Red$Cells_Principal %>% rownames_to_column(var="Sample") %>%  inner_join(X$cluster$Cluster_Quali, by="Sample"),
    color = ~Cluster,
    mode = 'markers',
    text = ~paste(Cluster," ", Sample)
    ,
    x = ~ X$Dim_Red$Cells_Principal[[input$Axis1_3D_Clus]],
    y = ~ X$Dim_Red$Cells_Principal[[input$Axis2_3D_Clus]],
    z = ~ X$Dim_Red$Cells_Principal[[input$Axis3_3D_Clus]],
    hoverinfo = "text",
    marker = list(
      opacity= input$Alpha_3D_Clus,
      size = input$Size_3D_Clus,
      line = list(color = 'rgba(0, 0, 0, .8)',
        width = 2)
      )
    ) %>% add_markers() %>%
  layout(
    autosize = F,
    legend = list(x = 100, y = 0.5),
    scene = list(
      xaxis = list(title = input$Axis1_3D_Clus),
      yaxis = list(title = input$Axis2_3D_Clus),
      zaxis = list(title = input$Axis3_3D_Clus)
      )
    )
  )

output$DTBOX<-renderDataTable({
  DTboxplot<-data_frame()
  Data<-X$cluster$Gene_Cluster_Distance %>%  gather("Cluster", "Distance", -Genes)
  for(i in (Data$Cluster %>% unique))
  {
    Cluster<-i
    bin1<-Data %>% filter(Cluster==i) %>%  arrange_("Distance") %>%  separate(Genes, into=c("Genes", "bin"), sep="-bin")  %>%  filter(bin==1)  %>%  extract2("Genes")  %>% head(5)  %>% paste(collapse=" ")
    bin2<-Data %>% filter(Cluster==i) %>%  arrange_("Distance") %>%  separate(Genes, into=c("Genes", "bin"), sep="-bin")  %>%  filter(bin==2)  %>%  extract2("Genes")  %>% head(5)  %>% paste(collapse=" ")
    DTboxplot<-bind_rows(DTboxplot,data_frame(Cluster,bin1,bin2))
  }
  return(DTboxplot %>% datatable(rownames = FALSE))})

Boxplot<-X$ExpressionMatrix %>%  data.frame %>%  rownames_to_column(var="Genes") %>% set_colnames(c("Genes",X$ExpressionMatrix %>%  colnames)) %>%  gather("Sample","Expression",-Genes)%>% as_tibble %>%  arrange(Sample)  %>%  inner_join(X$cluster$Cluster_Quali, by="Sample") %>% as_tibble
output$Boxplot <-
(renderPlotly(Boxplot %>% filter(Genes %in% input$Genes_Boxplot) %>%  plot_ly(x=~Genes, y=~Expression) %>% add_trace(type="box", color=~Cluster, jitter=0.5, pointpos=0, boxpoints="all", boxmean="sd") %>%  layout(boxmode="group",margin=list(b=100, t=25, l=50, r=50, pad=0) )))



#   ____________________________________________________________________________
#   GSEA server                                                             ####

observeEvent(input$Mode_Func_Plot,{
  updateSelectInput(session,"Choice_Func_Plot",
    label= paste("Choose", input$Mode_Func_Plot),
    choices={switch(input$Mode_Func_Plot,
      "Cluster"=X$Functionnal_Analysis$GSEA_Results %>%  names,
      "Axis"=X$Functionnal_Analysis$GSEA_Results_Axis %>%  names
      )})},ignoreNULL = TRUE)


Data<-reactive(X$Functionnal_Analysis$AllRanking[[input$Choice_Func_Plot]])
output$GSEA<-renderPlotly(plotlyEnrichment(X$Functionnal_Analysis$GMTfile[[input$Geneset]], Data(), gseaParam = X$Functionnal_Analysis$gseaParam))

      #Datatable of Enrichment results
      observeEvent(input$Mode_Func_DT,{
        updateSelectInput(session,"Choice_Func_DT",
          label= paste("Choose", input$Mode_Func_DT),
          choices={switch(input$Mode_Func_DT,
            "Cluster"=X$Functionnal_Analysis$GSEA_Results %>%  names,
            "Axis"=X$Functionnal_Analysis$GSEA_Results_Axis %>%  names
            )}
          )
        })
      Table_Enrich1<-reactive({switch(input$Mode_Func_DT, "Cluster"=(X$Functionnal_Analysis$GSEA_Results), "Axis"=(X$Functionnal_Analysis$GSEA_Results_Axis))})
      Table_Enrich2<-eventReactive(input$Choice_Func_DT,Table_Enrich1()[[input$Choice_Func_DT]] %>% select(-nMoreExtreme, -leadingEdge), ignoreNULL = TRUE)
      output$Table <- renderDataTable(Table_Enrich2()%>%  datatable(rownames=FALSE))

      Enrich_Boxplot<-X$Functionnal_Analysis$GSEA_Results %>%  map2(.y=X$Functionnal_Analysis$GSEA_Results %>% names, .f = function(x,y){mutate(.data=x, Cluster=y)}) %>%  bind_rows
      output$Enrich_Heatmap <-
      (renderPlotly(Enrich_Boxplot %>% filter(pathway %in% input$Enrich_Heatmap_GeneSet) %>%   group_by(pathway) %>%  plot_ly(x=~pathway, y=~Cluster, z=~(ES %>% as.numeric)) %>% add_heatmap(zmin=-1, zmax=1) %>%  layout(margin=list(
        b=100,
        t=25, l=150, r=150, pad=0
        ))
      ))

      output$downloadData <- downloadHandler(
  filename = function() {
    paste('data-', Sys.Date(), '.csv', sep='')
  },
  content = function(con) {
    write.csv(data, con)
  }
)
    }
    return(shinyApp(ui, server))
  }
