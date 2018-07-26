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
              fluidRow(
                box(width = 12, 
                    tabsetPanel(id = "input_tabset",
                                input_data_module()
                                #output_network_results_modules()
                                #input_ld_module(),
                                #input_network_module()
                    )
                )
              ),
              fluidRow(
                box(title = "Data", status = "primary", width = 12, 
                    fluidRow(
                      #input_data_module(), 
                      input_ld_module()
                      #input_network_module()
                    )
                )
              ),
              fluidRow(
                box(width = 12,
                    tabsetPanel(id = "results_tabset",
                                #output_raw_results_module(),
                                output_ld_results_module(),
                                output_network_results_modules(),
                                output_predictors_results_modules()
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