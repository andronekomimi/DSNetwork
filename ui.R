library(shinydashboard)
library(visNetwork)
library(shinyBS)
library(d3heatmap)

source('ui_modules.R')

dashboardPage(
  header = dashboardHeader(title = "DSNetwork"),
  sidebar = dashboardSidebar(
    sidebar_content()
  ),
  body = dashboardBody(
    tabItems(
      tabItem(tabName = "main",
              fluidRow(
                box(title = "Data", status = "primary", width = 12, 
                    fluidRow(
                      input_data_module(), 
                      input_ld_module(),
                      input_network_module()
                    )
                )
              ),
              fluidRow(
                box(width = 12,
                    tabsetPanel(id = "results_tabset",
                                output_raw_results_module(),
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