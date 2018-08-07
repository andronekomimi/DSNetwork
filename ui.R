library(shinydashboard)
library(visNetwork)
library(shinyBS)
library(d3heatmap)

source('ui_modules.R')

dashboardPage(
  header = dashboardHeader(title = "DSNetwork"),
  sidebar = dashboardSidebar(
    collapsed = TRUE,
    sidebar_content()
  ),
  body = dashboardBody(
    tags$head(tags$style(
      HTML(".shiny-notification {
             position:fixed;
           top: calc(50%);;
           left: calc(50%);;
            width: calc(25%);;
           }
           "
      )
    )),
    tabItems(
      tabItem(tabName = "main",
              box(width = "100%", height = "100%",
                  fluidRow(
                    column(width = 4, 
                           input_data_module(),
                           input_network_module(),
                           ld_mapping_module()
                    ),
                    column(width = 8, 
                           conditionalPanel(condition="input.fetch_annotations",
                                            output_plot_row(),
                                            network_results_modules()
                           )
                           
                    )
                  ),
                  hr(),
                  fluidRow(
                    column(width = 12,
                           conditionalPanel(condition="input.buildNetwork",
                                            network_modifiers_row()
                           )
                    )
                  )
              )
      ),
      tabItem(tabName = "readme",
              includeMarkdown("README.Rmd")
      )
    )
  ), 
  title = "DSNetwork",
  skin = "black"
)