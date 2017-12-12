TabMCA <- function(X, dr_axis) {
  tabItem(
    tabName = "mca",
    navbarPage(
      "MCA",
      navbarMenu(
        "Cell Space",
        tabPanel(
          "2 Axis",
          wellPanel(fluidRow(
            column(
              6,
              selectInput(
                "MCA_CS_2D_Type",
                label = "Type",
                choices = c("Principal", "Standard"),
                selected = "Principal"
              ),
              selectInput(
                "MCA_CS_Axis_x",
                label = "Select x Axis",
                choices = dr_axis,
                selected = "Axis1"
              ),
              selectInput(
                "MCA_CS_Axis_y",
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
          mainPanel(
            width = 12,
            plotlyOutput(
              "CellSpace2D",
              width = "95%", height = "100%"
            )
          )
        ),
        tabPanel(
          "3 Axis",
          wellPanel(fluidRow(
            column(
              4,
              selectInput(
                "MCA_CS_3D_Type",
                label = "Type",
                choices = c("Principal", "Standard"),
                selected = "Principal"
              ),
              selectInput(
                "MCA_CS_Axis1_3D",
                label = "Select x Axis",
                choices = dr_axis,
                selected = "Axis1"
              )
            ),
            column(
              4,
              selectInput(
                "MCA_CS_Axis2_3D",
                label = "Select y Axis",
                choices = dr_axis,
                selected = "Axis2"
              ),
              selectInput(
                "MCA_CS_Axis3_3D",
                label = "Select z Axis",
                choices = dr_axis,
                selected = "Axis3"
              )
            ),
            column(
              4,
              sliderInput(
                "MCA_CS_Size_3D",
                label = "Point Size",
                min = 0,
                max = 10,
                value = 5,
                step = 0.1
              ),
              sliderInput(
                "MCA_CS_Alpha_3D",
                label = "Transparency",
                min = 0,
                max = 1,
                value = 1,
                step = 0.1
              )
            )
          )),
          mainPanel(
            width = 12,
            align = "center",
            plotlyOutput("CellSpace3D", width = "100%", height = "100%")
          )
        ),
        tabPanel(
          "Axis & Gene Correlation",
          wellPanel(fluidRow(
            column(
              4,
              selectInput(
                "MCA_CS_AC_Type",
                label = "Type",
                choices = c("Principal", "Standard"),
                selected = "Principal"
              ),
              selectInput(
                "MCA_CS_AC_Gene",
                label = "Gene",
                choices = X$ExpressionMatrix %>% rownames()
              )
            ),
            column(
              4,
              selectInput(
                "MCA_CS_AC_Axis_x",
                label = "Select x Axis",
                choices = X$MCA$cells_principal %>%
                  rownames_to_column(var = "Sample") %>%
                  select(contains("Axis")) %>%
                  colnames(),
                selected = "Axis1"
              ),
              selectInput(
                "MCA_CS_AC_Axis_y",
                label = "Select y Axis",
                choices = X$MCA$cells_principal %>%
                  rownames_to_column(var = "Sample") %>%
                  select(contains("Axis")) %>%
                  colnames(),
                selected = "Axis2"
              )
            ),
            column(
              4,
              sliderInput(
                "MCA_CS_AC_Size",
                label = "Point Size",
                min = 0,
                max = 10,
                value = 5,
                step = 1
              ),
              sliderInput(
                "MCA_CS_AC_Alpha",
                label = "Transparency",
                min = 0,
                max = 1,
                value = 1,
                step = 0.1
              )
            )
          )),
          mainPanel(
            width = 12,
            plotlyOutput("CellSpaceGeneCor", width = "100%", height = "100%")
          ),
          mainPanel(
            width = 12,
            DT::dataTableOutput("TableGeneCor", width = "100%", height = "100%")
          )
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
              ),
              selectInput(
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
            5, offset = 2, wellPanel(
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
        mainPanel(width = 12, plotlyOutput("GeneSpace"))
      ),
      tabPanel(
        title = "Eigen Value",
        selectInput(
          "Mode",
          label = "Select Mode",
          choices = c(
            Explained_Variance = "Explained_Variance",
            Cumulative = "Cumulative"
          ),
          selected = "Explained_Variance"
        ),
        sliderInput(
          "amount_adjust",
          label = "Amount",
          min = 1,
          max = X$MCA$cells_principal %>% ncol(),
          value = if (X$MCA$cells_principal %>% ncol() %>% is_greater_than(5)) {
            5
          } else {
            1
          },
          step = 1
        ),
        plotlyOutput("Eigen")
      )
    )
  )
}

TabClus <- function(X, dr_axis) {
  tabItem(
    tabName = "clus",
    navbarPage(
      "Clustering",
      navbarMenu(
        "Cell Space",
        tabPanel(
          "Clustering 2 Axis",
          titlePanel("Cluster in the Cell Space"),
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
        ),
        tabPanel(
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
          )),
          mainPanel(
            width = 12,
            align = "center",
            plotlyOutput("CellSpace3D_Clus", width = "100%", height = "100%")
          )
        )
      ),
      tabPanel(
        title = "Gene Space",
        titlePanel("Genespace with Cluster Centroids"),
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
        mainPanel(width = 12, plotlyOutput("GeneSpace_Clus"))
      ),
      navbarMenu(
        title = "Expression",
        tabPanel(
          title = "Boxplot",
          fluidPage(titlePanel("Genes Expression by Cluster"), fluidRow(column(
            width = 6,
            selectInput(
              "Genes_Boxplot",
              "Choose a Gene:",
              choices = (X$ExpressionMatrix %>%
                rownames() %>%
                sort()),
              selectize = TRUE,
              multiple = TRUE,
              selected = (X$ExpressionMatrix %>%
                rownames())[1]
            )
          ))),
          mainPanel(
            width = 12,
            h3("Boxplot"),
            plotlyOutput(
              "Boxplot",
              width = "100%", height = "100%"
            ),
            h3("Genes Specific to Cluster"),
            DT::dataTableOutput("DTBOX")
          )
        ),
        tabPanel(
          "Expression Heatmap",
          wellPanel(fluidRow(column(
            5,
            sliderInput(
              "Num_Top_Genes",
              label = "Number of Genes",
              min = 1,
              max = 10,
              value = 3,
              step = 1
            )
          ))),
          mainPanel(
            width = 12,
            align = "center",
            plotlyOutput(
              "Heatmap_Expression_Clus",
              width = "100%",
              height = "100%"
            )
          )
        )
      )
    )
  )
}


TabGSEA <- function(X, dr_axis) {
  tabItem(
    tabName = "gsea",
    navbarPage(
      "GSEA",
      tabPanel(
        title = "Enrichmentplot",
        fluidPage(
          titlePanel("Enrichmentplot"),
          wellPanel(fluidRow(
            column(
              width = 6,
              selectInput(
                inputId = "Geneset",
                label = "Choose a Geneset:",
                choices = (X$GSEA$Pathways),
                selectize = TRUE
              )
            ),
            column(
              6,
              selectInput(
                inputId = "Choice_Func_Plot",
                label = "Choose a Cluster:",
                choices = X$GSEA$GSEA_Results %>% names()
              ),
              selectInput(
                "Mode_Func_Plot",
                "Choose", choices = c("Cluster", "Axis")
              )
            )
          )),
          mainPanel(width = 12, plotlyOutput(
            "GSEA",
            width = "100%",
            height = "100%"
          ))
        )
      ),
      tabPanel(
        title = "Enrichment Results",
        fluidPage(titlePanel("Enrichment Table"), wellPanel(fluidRow(
          column(
            width = 6,
            selectInput(
              inputId = "Mode_Func_DT",
              label = "Choose Data:",
              choices = (c("Cluster", "Axis")),
              selectize = TRUE
            )
          ),
          column(
            width = 6,
            selectInput(
              inputId = "Choice_Func_DT",
              label = "Choose a Geneset:",
              choices = (X$GSEA$GSEA_Results %>% names()),
              selectize = TRUE
            )
          )
        ))),
        mainPanel(
          width = 12,
          DT::dataTableOutput("Table", width = "100%", height = "100%")
        )
      ),
      tabPanel(
        title = "GSEA Heatmap",
        fluidPage(
          titlePanel("Genes Enrichment by Cluster"),
          wellPanel(fluidRow(
            column(
              width = 5,
              sliderInput(
                inputId = "GSEA_N_Path",
                label = "Number of Pathway to represent",
                min = 1,
                max = 50,
                value = 5,
                step = 1
              ),
              selectInput(
                inputId = "GSEA_Heatmap_Metrics",
                label = "Metrics:",
                choices = c("NES", "ES"),
                selectize = TRUE
              ),
              selectInput(
                inputId = "GSEA_Heatmap_Row_Side_Color",
                label = "Row Side Color",
                choices = colvector <- c("Yellow_Red", "Grey", "RG", "Spectral", "Cool_Warm"),
                selectize = TRUE
              ),
              selectInput(
                inputId = "GSEA_Heatmap_Color",
                label = "Color",
                choices = colvector <- c("Yellow_Red", "Grey", "RG", "Spectral", "Cool_Warm"),
                selectize = TRUE
              ),
              selectInput(
                inputId = "GSEA_Heatmap_Level",
                label = "Level",
                choices = {
                  if (X$SC_GSEA %>% is.null() %>% not()) {
                    c("Cluster", "Cell-Cluster", "Cell")
                  } else {
                    "Cluster"
                  }
                },
                selectize = TRUE
              )
            ),
            column(
              width = 7,
              fluidRow(
                column(
                  width = 4,
                  numericInput(
                    inputId = "GSEA_Heatmap_Pval",
                    label = "Adjusted Pvalue:",
                    min = 0,
                    max = 1,
                    step = 0.05,
                    value = 0.05
                  )
                ),
                column(
                  width = 4,
                  numericInput(
                    inputId = "GSEA_Heatmap_ES",
                    label = "ES Threshold",
                    min = -1,
                    max = 1,
                    step = 0.1,
                    value = 0
                  )
                ),
                column(
                  width = 4,
                  numericInput(
                    inputId = "GSEA_Heatmap_NES",
                    label = "NES Threshold",
                    min = -Inf,
                    max = Inf,
                    value = 0,
                    step = 0.5
                  )
                )
              ),
              fluidRow(
                column(
                  width = 6,
                  numericInput(
                    "GSEA_Heatmap_Cex_Col",
                    "Col Character Size:",
                    value = 1,
                    step = 0.1,
                    min = 0,
                    max = 10
                  )
                ),
                column(
                  width = 6,
                  numericInput(
                    "GSEA_Heatmap_Cex_Row",
                    "Row Character Size:",
                    value = 1,
                    step = 0.1,
                    min = 0,
                    max = 10
                  )
                )
              )
            )
          ))
        ),
        mainPanel(
          width = 12,
          h3("GSEA Enrichment Heatmap"),
          plotlyOutput(
            "GSEA_Heatmap", width = "90%",
            height = "90%"
          )
        )
      )
    )
  )
}
