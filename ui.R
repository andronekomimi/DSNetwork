require(shinydashboard)
require(visNetwork)
require(shinyBS)
#require(d3heatmap)
require(shinyjs)
require(shinyalert)
#require(markdown)
require(ggplot2)

source('ui_modules.R')

jscode <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
};
shinyjs.myFunction = function(idx) {
    alert('Select '+idx);
};
shinyjs.reset = function() {history.go(0)}
"

dashboardPage(
  header = dashboardHeader(title = "DSNetwork"),
  sidebar = dashboardSidebar(
    collapsed = TRUE,
    sidebar_content()
  ),
  body = dashboardBody(
    shinyalert::useShinyalert(),
    useShinyjs(),
    extendShinyjs(text = jscode, functions = c("collapse","myFunction","reset")),
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
              presentation_module(),
              fluidRow(box(id = "input_box", title = "Request panel", status = "primary",
                           solidHeader = FALSE, background = NULL, width = 12,
                           collapsible = TRUE, collapsed = FALSE,
                           fluidRow(
                             column(width = 3, 
                                    conditionalPanel(condition = "input.fetch_annotations == 0",
                                                     input_data_module()
                                    ),
                                    conditionalPanel(condition = "input.fetch_annotations == 1",
                                                     bsButton(inputId = "reload", 
                                                              label = "Reset for new query", 
                                                              icon = icon("redo"), disabled = FALSE),
                                                     br(),
                                                     br(),
                                                     bsAlert("alert_conv"),
                                                     bsAlert("alert_res"))
                             ),
                             column(width = 9, 
                                    textOutput(outputId = 'can_run'),
                                    conditionalPanel(condition="output.can_run == 'CDTS plot'",
                                                     output_plot_row()
                                    )       
                             )
                           )
                           
              )), # end request panel
              fluidRow(box(id = "selection_box", title = "Selection panel", status = "info",
                           solidHeader = FALSE, background = NULL, width = 12,
                           collapsible = TRUE, collapsed = FALSE,
                           conditionalPanel(condition="output.can_run == 'CDTS plot'",
                                            fluidRow(
                                              column(width = 3, 
                                                     selection_module()
                                              ),
                                              column(width = 9, 
                                                     raw_results_row()
                                              )
                                            )
                           )
                           
              )), # end selection panel
              fluidRow(box(title = "Network panel", status = "success",
                           solidHeader = FALSE, background = NULL, width = 12,
                           collapsible = TRUE, collapsed = FALSE,
                           conditionalPanel(condition="input.buildNetwork",
                                            fluidRow(
                                              column(width = 3, 
                                                     focus_module(),
                                                     ld_mapping_module(),
                                                     nodes_modifiers_box()
                                              ),
                                              column(width = 8, 
                                                     network_results_modules()
                                              ),
                                              column(width = 1,
                                                     scales_module()
                                              )
                                              
                                            )
                           )
                           
              )) # end network
              
      ),
      # tabItem(tabName = "readme",
      #         shiny::includeHTML(path = "README.Rhtml")
      # ),
      tabItem(tabName = "description",
              score_desc_module()
      ),
      tabItem(tabName = "about",
              shiny::includeHTML("ABOUT.Rhtml")
      )
    ),#end dashboardPage
    tags$footer(htmlOutput(outputId = "version_footer"), align = "center")
  ), 
  title = "DSNetwork",
  skin = "black"
)
