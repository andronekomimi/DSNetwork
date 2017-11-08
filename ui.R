require(shinydashboard)

dashboardPage(
  header = dashboardHeader(title = "DSNetwork"),
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("Data", tabName = "load_data", icon = icon("th")),
      menuItem("Network", tabName = "network", icon = icon("snowflake-o"), selected = T)
    )
  ),
  body = dashboardBody(
    tabItems(
      tabItem(tabName = "load_data",
              h3("Data loading"),
              tags$ul(
                tags$li("Load new dataframe"), 
                tags$li("Load custom object")
              ),
              h3("Data processing"),
              tags$ul(
                tags$li("Get snp position from id (hg19)"), 
                tags$li("Run ANNOVAR annotation"), 
                tags$li("Fetch LD information"),
                tags$li("Combine for network visualization")
              ),
              h3("Save data for external or further usage (Recommended)")
      ),
      tabItem(tabName = "network",
              fluidRow(
                column(width = 9,
                       box(width = NULL, title = "Network", collapsible = TRUE, 
                           visNetworkOutput("network_hello",
                                            height = "800", width = "auto")
                       ),
                       box(width = NULL, status = "info", title = "Legends", height = "auto",
                           collapsible = T, 
                           fluidRow(
                             column(width = 12, plotOutput(outputId = "color_key", height = "200"))
                           )
                       )
                ),
                column(width = 3,
                       box(width = NULL, status = "warning",
                           checkboxGroupInput("annotations", "Predictors",
                                              choices = c(
                                                CADD13_PHRED = 1,
                                                Eigen = 2,
                                                FATHMM_noncoding = 3,
                                                `gerp++gt2` = 4,
                                                GWAVA_region_score = 5,
                                                LINSIGHT = 6
                                              ),
                                              selected = 1:6
                           ),
                           p(
                             class = "text-muted",
                             paste("This option enables to select the set of prediction and",
                                   "scoring algorithms used to create the network."
                             )
                           ),
                           actionButton("update_annotations", "Update")
                       ),
                       box(width = NULL, status = "warning",
                           selectInput("focus", "Focus on",
                                       choices = c(
                                         "None" = -1
                                       ),
                                       selected = -1
                           ),
                           p(class = "text-muted",
                             br(),
                             "This option enables to focus the network on a particular variant"
                           )
                       ),
                       box(width = NULL, status = "warning",
                           sliderInput("ld_range", "LD range",
                                       min = 0, max = 1, value = c(0, 1)),
                           p(class = "text-muted",
                             br(),
                             "This option enables to select the interval of LD values represented between variants"
                           ),
                           actionButton("update_ld", "Update")
                       )
                )
              )
      )
    )
  ), 
  title = "DSNetwork",
  skin = "black"
)