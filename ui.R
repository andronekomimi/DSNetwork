require(shinydashboard)
require(visNetwork)
require(shinyBS)
#require(d3heatmap)
require(shinyjs)
require(shinyalert)
#require(markdown)

source('ui_modules.R')

jscode <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
};
shinyjs.myFunction = function(idx) {
    alert('Select '+idx);
}
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
    extendShinyjs(text = jscode, functions = c("collapse","myFunction")),
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
                             column(width = 4, 
                                    input_data_module()
                             ),
                             column(width = 8, 
                                    conditionalPanel(condition="input.fetch_annotations",
                                                     output_plot_row()
                                    )       
                             )
                           )
                           
              )), # end request panel
              fluidRow(box(id = "selection_box", title = "Selection panel", status = "info",
                           solidHeader = FALSE, background = NULL, width = 12,
                           collapsible = TRUE, collapsed = FALSE,
                           conditionalPanel(condition="input.fetch_annotations",
                                            fluidRow(
                                              column(width = 4, 
                                                     selection_module()
                                              ),
                                              column(width = 8, 
                                                     raw_results_row()
                                              )
                                            )
                           )
                           
              )), # end selection panel
              fluidRow(box(title = "Network panel", status = "success",
                           solidHeader = FALSE, background = NULL, width = 12,
                           collapsible = TRUE, collapsed = FALSE,
                           conditionalPanel(condition="input.fetch_annotations",
                                            fluidRow(
                                              column(width = 4, 
                                                     box(width = "100%",
                                                         ld_mapping_module(),
                                                         nodes_modifiers_box()
                                                     ) 
                                              ),
                                              column(width = 8, 
                                                     network_results_modules()
                                              )
                                            )
                           )
                           
              )) # end network
              
      ),
      tabItem(tabName = "readme",
              shiny::includeHTML(path = "README.Rhtml")
      ),
      tabItem(tabName = "description",
              score_desc_module()
      ),
      tabItem(tabName = "about",
              shiny::includeHTML("ABOUT.Rhtml")
      )
    )
  ), 
  title = "DSNetwork",
  skin = "black"
)