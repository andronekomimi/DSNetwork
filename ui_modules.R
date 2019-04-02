require(shinydashboard)
require(visNetwork)
require(shinyBS)
#require(d3heatmap)
require(plotly)
#require(ggrepel)
require(shinyjqui)
require(shinyjs)

source('helper.R', local = TRUE)

populations <- c("African Caribbean in Barbados (ACB)" = 'ACB',
                 "African Ancestry in Southwest US (ASW)" = 'ASW',
                 "Bengali in Bangladesh (BEB)" = 'BEB',
                 "Chinese Dai in Xishuangbanna, China (CDX)" = 'CDX',
                 "Utah residents with Northern and Western Europen ancestry (CEU)" = 'CEU',
                 "Han Chinese in Bejing, China (CHB)" = 'CHB',
                 "Southern Han Chinese, China (CHS)" = 'CHS',
                 "Colombian in Medellin, Colombia (CLM)" = 'CLM',
                 "Esan in Nigeria (ESN)" = 'ESN',
                 "Finnish in Finland (FIN)" = 'FIN',
                 "British in England and Scotland (GBR)" = 'GBR',
                 "Gujarati Indian in Houston, TX (GIH)" = 'GIH',
                 "Iberian populations in Spain (IBS)" = 'IBS',
                 "Indian Telugu in the UK (ITU)" = 'ITU',
                 "Japanese in Tokyo, Japan (JPT)" = 'JPT',
                 "Kinh in Ho Chi Minh City, Vietnam (KHV)" = 'KHV',
                 "Luhya in Webuye, Kenya (LWK)" = 'LWK',
                 "Mandinka in The Gambia (MAG)" = 'MAG',
                 "Mende in Sierra Leone (MSL)" = 'MSL',
                 "Mexican Ancestry in Los Angeles, California (MXL)" = 'MXL',
                 "Peruvian in Lima, Peru (PEL)" = 'PEL',
                 "Punjabi in Lahore, Pakistan (PJL)" = 'PJL',
                 "Puerto Rican in Puerto Rico (PUR)" = 'PUR',
                 "Sri Lankan Tamil in the UK (STU)" = 'STU',
                 "Toscani in Italy (TSI)" = 'TSI',
                 "Yoruba in Ibadan, Nigeria (YRI)" = 'YRI')

pie_types <- list(
  `Intra-predictor ranking` = c("Relative rank" = 'pie_scores',
                                "Relative rank (group by color)" = 'pie_scores_group',
                                "Absolute rank" = 'pie_scores_abs',
                                "Absolute rank (group by color)" = 'pie_scores_abs_group'),
  `Global ranking` = c("Mean relative rank (NA last)" = 'pie_rank_na_last',
                       "Mean relative rank (NA mean)" = 'pie_rank_na_mean',
                       "Mean relative rank (NA median)" = 'pie_rank_na_median')
)

sidebar_content <- function(){
  sidebarMenu(
    menuItem("Get your network", tabName = "main", icon = icon("hand-o-up")),
    menuItem("Scores description", tabName = "description", icon = icon("info-circle")),
    menuItem("Read Me", tabName = "readme", icon = icon("mortar-board")),
    menuItem("About", tabName = "about", icon = icon("question"))
  )
}

input_data_module <- function(){
  box(width = "100%",
      textAreaInput("query", "Enter variant ids", "", rows = 5, 
                    placeholder = "Please enter one variant id per line (rs123455 or 1:1324:A:C)"),
      div(style = "text-align:-webkit-right", 
          actionLink(inputId = "load_demo1", label = "Load 1p36 data, "),
          actionLink(inputId = "load_demo2", label = "load 1p34 data, "),
          actionLink(inputId = "load_demo3", label = "load 7q22 data, "),
          actionLink(inputId = "load_demo4", label = "load 11p15 data.")),
      br(),
      fileInput("query_file", "or load text file (one variant id per line)", 
                multiple = FALSE, 
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      checkboxInput(inputId = "fetch_snpnexus", label = "Fetch annotations from SNPnexus (significatively increases fetching duration)"),
      conditionalPanel(condition = "input.fetch_snpnexus == 1",
                       sliderInput(inputId = "waiting", 
                                   label = "How long are you willing to wait ? (default: 5 min)", 
                                   min = 1, max = 10, value = 5, step = 1)),
      bsButton(inputId = "fetch_annotations", 
               label = "Fetch Annotations", 
               icon = icon("search"), disabled = FALSE),
      conditionalPanel(condition = "input.fetch_annotations",
                       br(),
                       bsAlert("alert_conv"),
                       bsAlert("alert_res"),
                       downloadButton('downloadRawTable', 'Download results (TSV)'))
  )
}

output_plot_row <- function(){
  list(
    div(style = "color:gray", 
        HTML(paste0("This plot represents the requested variants along the map of sequence constraint "),
             "- <b>C</b>ontect-<b>D</b>ependent <b>T</b>olerence <b>S</b>core (CDTS) - ",
             "determined throught alignment of thousands of human genomes.")),
    shinyjqui::jqui_resizable(plotly::plotlyOutput(outputId = "my_plot",
                                                   height = "400px", width = "auto"))
  )
}

raw_results_row <- function(){
  list(
    DT::dataTableOutput("raw_data"),
    div(style = "text-align:-webkit-right", 
        bsButton(inputId = 'buildNetwork', label = 'Build Network',
                 disabled = TRUE, icon = icon("gear")))
  )
}

selection_module <- function(){
  box(width = "100%",
      selectInput(inputId = 'network_type',
                  label = 'Choose the network to build', selected = "regul",
                  choices = c("non synonymous variants" = "non_syn",
                              "synonymous and non-coding variants" = "regul"),
                  width = "100%")
  )
}

nodes_modifiers_box <- function(){
  list(
    selectInput(inputId = "snv_nodes_type",
                label = "Prediction visualization",
                multiple = FALSE,
                choices = pie_types,
                selected = 'pie_scores'
    ),
    conditionalPanel(condition="input.snv_nodes_type!='metascores'",
                     fluidRow(
                       column(width = 12,
                              selectInput("selected_scores", 
                                          'Predictors selection',
                                          choices = c(),
                                          selectize = TRUE, multiple = TRUE))
                     ),
                     bsButton(inputId = 'update_metascore', label = 'Update',  
                              icon = icon("gear"), disabled = TRUE),
                     p(
                       class = "text-muted", br(),
                       paste("This option enables to select the set of prediction and",
                             "scoring algorithms used to compute the metascore (color of the database-shaped nodes)"
                       )
                     )
    )
  )
}

# network_modifiers_row <- function(){
#   box(width = NULL, status = "warning", 
#       height = "500px", collapsible = TRUE,
#       textOutput(outputId = "snv_score_details_id"),
#       DT::dataTableOutput(outputId = "snv_score_details")
#   )
# }

network_results_modules <- function(){
  list(
    shinyjqui::jqui_resizable(visNetworkOutput("my_network", height = "600px",
                                               width = "100%"))
  )
}

ld_mapping_module <- function(){
  list(
    selectInput(inputId = 'population',
                label = 'Choose the population to use',
                choices = populations),
    fluidRow(column(width = 6, 
                    bsButton(inputId = 'runLD', label = 'Add LD Information',  
                             icon = icon("gear"), disabled = TRUE)),
             column(width = 6, 
                    bsButton(inputId = 'removeLD', label = 'Remove LD Information',  
                             icon = icon("gear"), disabled = TRUE))
    ),
    br(),
    bsAlert(anchorId = "alert_ld"),
    sliderInput("ld_range", "LD range",
                min = 0, max = 1, value = c(0, 1), step = 0.1), 
    bsButton(inputId = 'update_ld', label = 'Update',  
             icon = icon("gear"), disabled = TRUE),
    p(class = "text-muted",
      br(),
      "This option enables to select the interval of LD values represented between variants"
    )
  )
}

presentation_module <- function(){
  fluidRow(box(id = "intro_box", title = "DSNetwork", solidHeader = FALSE, background = NULL, width = 12,
               collapsible = TRUE, collapsed = TRUE,
               column(width = 12, 
                      helpText("DSNetwork provides a user-friendly interface integrating predictors for both coding and non-coding variants in an easy-to-interpret visualization to assist prioritization process. The usage of DSNetwork greatly facilitate the selection process by aggregating the results from nearly sixty prediction approaches and easily highlights the best candidate variants for further functional analysis."),
                      helpText("Each node corresponds to an annotated variant and edges between the nodes can be used to materialize LD levels between two variants.
Nodes filling are dedicated to prediction scores display. The default display is a pie chart where each slice represents the variant score for a particular predictor.
                               The slice color gradient from green to red within each predictor reflects the ranking of this variant compared to the other variants in the selected list, a red slice indicating a variant more likely to be damaging with regard to a particular predictor.
                               LD levels are mapped upon request and based on a user-chosen 1000 genomes population after building the network. LD (squared correlation R2) is represented by an absolute color gradient from yellow to red. Red indicating a high disequilibrium. We did not use the usual white-to-red gradient to preserve the contrast with the white background.
                               ")
               )
  ))
}

score_desc_module <- function(){
  box(id = "description", title = "Predictors", solidHeader = FALSE, background = NULL, width = 12,
      collapsible = FALSE, collapsed = FALSE,
      DT::dataTableOutput("scores_description_table"),
      bsModal("annotModal", "Annotations", "", tags$p("Annotations"), size = "small")
  )
}

scales_module <- function(){
  list(
    plotOutput(outputId = "scale",
               #height = "600px",
               width = "100px"),
    plotOutput(outputId = "ld_scale",
               #height = "600px",
               width = "100px")
  )
}

#### OBSOLETE ####
# predictors_selection <- function(){
#   fluidRow(
#     column(width = 12,
#            selectInput("predictors", 'Predictors',width = "50%",
#                        choices = c(), selectize = TRUE, multiple = TRUE),
#            bsButton(inputId = 'get_predictors_info', label = "Get info",
#                     icon = icon("search"))
#     )
#   )
# }
# output_predictors_results_modules <- function(){
#   tabPanel(title = h5("Scores information"),
#            value = "scores_stats",
#            conditionalPanel(condition="input.buildNetwork",
#                             predictors_selection(),
#                             output_predictors_descript()
#            )
#   )
# }
# 
# output_predictors_descript <- function(){
#   fluidRow(
#     column(width = 6,
#            box(width = NULL, status = "success",
#                collapsible = TRUE,
#                plotOutput(outputId = "scores_stats",height = "500px")
#            ),
#            box(width = NULL, status = "success",
#                d3heatmapOutput(outputId = "scores_corr", height = "600px")
#            )
#     ),
#     column(width = 6,
#            box(width = NULL, status = "success",
#                DT::dataTableOutput(outputId = "scores_missing_data")
#            )
#     )
#   )
# }
# 
# predictors_modifiers_box <- function(){
#   box(width = NULL, status = "warning",
#       height = "500px", collapsible = TRUE,
#       fluidRow(
#         column(width = 12,
#                selectInput("selected_scores", 'Predictors',
#                            choices = c(),
#                            selectize = TRUE, multiple = TRUE))
#       ),
#       actionButton("update_metascore", "Update"),
#       p(
#         class = "text-muted", br(),
#         paste("This option enables to select the set of prediction and",
#               "scoring algorithms used to compute the metascore (color of the database-shaped nodes)"
#         )
#       )
#   )
# }
# 
# edges_modifiers_box <- function(){
#   box(width = NULL, status = "warning",  
#       height = "500px", collapsible = TRUE,
#       selectInput("snv_edges_type", "SNV Edges",
#                   choices = c(
#                     "Distance" = 0,
#                     "Linkage" = 1
#                   ),
#                   selected = 0
#       ),
#       conditionalPanel(condition="input.snv_edges_type=='0'",
#                        sliderInput("dist_range", "Distance (kb)",
#                                    min = 0, max = 1000, step = 1, value = 1000),
#                        actionButton("update_dist", "Update"),
#                        p(class = "text-muted",
#                          br(),
#                          "This option enables to select the interval of 
#                          distance represented between variants"
#                        )
#                        ),
#       conditionalPanel(condition="input.snv_edges_type=='1'",
#                        selectInput(inputId = 'population',
#                                    label = 'Choose the population to use',
#                                    choices = populations),
#                        bsButton(inputId = 'runLD', label = 'Add LD Information',  
#                                          icon = icon("gear"), disabled = FALSE),
#                        br(),
#                        br(),
#                        verbatimTextOutput("runld_res"),
#                        tags$style(type="text/css", "#runld_res {white-space: pre-wrap;}"),
#                        br(),
#                        sliderInput("ld_range", "LD range",
#                                    min = 0, max = 1, value = c(0, 1)),
#                        actionButton("update_ld", "Update"),
#                        p(class = "text-muted",
#                          br(),
#                          "This option enables to select the interval of LD values represented between variants"
#                        )
#       )
#       )
# }
# 
# output_raw_results_module <- function(){
#   column(width = 4, 
#          br(),
#          conditionalPanel(condition="input.fetch_annotations",
#                           #DT::dataTableOutput("raw_data"),
#                           fluidRow(column(width = 12, C3PieChartOutput(outputId = "consequences"))),
#                           fluidRow(column(width = 12, C3PieChartOutput(outputId = "annotations"))),
#                           fluidRow(column(width = 12, C3PieChartOutput(outputId = "genenames")))
#                           
#          )
#   )
# } 
# 
# output_raw_results_module_2 <- function(){
#   column(width = 4, 
#          br(),
#          conditionalPanel(condition="input.fetch_annotations",
#                           fluidRow(column(width = 12, C3PieChartOutput(outputId = "cdts_scores", height = "400px"))),
#                           fluidRow(column(width = 12, C3PieChartOutput(outputId = "cdts_percentile")))
#                           
#          )
#   )
# }
# 
# output_ld_results_module <- function(){
#   tabPanel(title = h5("LD Plots"),
#            value = "ld_results",
#            conditionalPanel(condition="input.runLD",
#                             selectInput("ld_regions", "LD Regions",
#                                         width = "20%",
#                                         choices = c()
#                             ),
#                             imageOutput("ld_plot", height = 600)
#            )
#   )
# }
#
# output_plot_row.old <- function(){
#   fluidRow(
#     column(width = 12,
#            br(),
#            shinyjqui::jqui_resizable(plotlyOutput(outputId = "my_plot",
#                                        height = "400px", width = "auto")),
#            box(title = "Variant selection",
#                width = 12,
#                collapsible = TRUE, 
#                collapsed = FALSE, 
#                htmlOutput("selection")))
#   )
# }