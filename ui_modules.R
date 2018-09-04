library(shinydashboard)
library(visNetwork)
library(shinyBS)
library(d3heatmap)
require(plotly)
require(ggrepel)
library(shinyjqui)
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

preload_loci <- c("FGFR2" = "locus_0",
                  "Locus 2 [chr1:113948389-14948389]" = "locus_2",
                  "Locus 6 [chr1:201687176-202687176]" = "locus_6",
                  "Locus 38 [chr7:143574929-144574929]" = "locus_38",
                  "Locus 60 [chr13:32468810-33472626]" = "locus_60",
                  "Locus 70 [chr17:77281387-78281725]" = "locus_70",
                  "Locus 78 [chr22:40376234-41527870]" = "locus_78",
                  "Locus 80 [chr3:86537543-8753754]" = "locus_80"
)

sidebar_content <- function(){
  sidebarMenu(
    menuItem("How to get a network", tabName = "main", icon = icon("hand-o-up")),
    helpText("Help text"),
    menuItem("ReadMe", tabName = "readme", icon=icon("mortar-board")),
    menuItem("About", tabName = "about", icon = icon("question"))
  )
}

input_data_module <- function(){
  box(width = "100%",
      textAreaInput("query", "Enter variant ids", "", rows = 5, 
                    placeholder = "Please enter one variant id per line (rs123455 or 1:1324:A:C)"),
      fileInput("query_file", "or load text file (one variant id per line)", 
                multiple = FALSE, 
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      checkboxInput(inputId = "fetch_snpnexus", label = "Fetch annotations from SNPnexus (significatively increase fetching duration)"),
      conditionalPanel(condition = "input.fetch_snpnexus == 1",
                       sliderInput(inputId = "waiting", 
                                   label = "How long are you willing to wait ?", 
                                   min = 1, max = 10, value = 2, step = 1),
                       verbatimTextOutput("snpnexus_res"),
                       tags$style(type="text/css", "#snpnexus_res {white-space: pre-wrap;}")),
      shinyBS::bsButton(inputId = "fetch_annotations", 
                        label = "Fetch Annotations", 
                        icon = icon("search"), disabled = FALSE),
      verbatimTextOutput("transform_res"),
      tags$style(type="text/css", "#transform_res {white-space: pre-wrap;}"),
      br(),
      verbatimTextOutput("query_res"),
      tags$style(type="text/css", "#query_res {white-space: pre-wrap;}"),
      selectInput(inputId = 'preload',
                  label = '0) Load preset query',
                  choices = preload_loci),
      actionButton("preload_loci", "Load preset query", icon = icon("caret-right"))
  )
}

output_plot_row <- function(){
  list(
    helpText("Qu'est-ce qu'on regarde bordel ??!"),
    jqui_resizable(plotlyOutput(outputId = "my_plot",
                                height = "400px", width = "auto"))
  )
}

raw_results_row <- function(){
  list(
    DT::dataTableOutput("raw_data"),
    fluidRow(column(width = 6, downloadButton('downloadRawTable', 'Download all (TSV)')),
             column(width = 6, shinyBS::bsButton(inputId = 'buildNetwork', label = 'Build Network',
                                                 disabled = FALSE, icon = icon("gear"))))
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
                label = "5) Predictors selection",
                multiple = FALSE,
                choices = c(
                  "Scores Pie" = 'pie_scores',
                  list(
                    `Relative metascores` = c("Rank (NA last)" = 'pie_rank_na_last',
                                              "Rank (NA mean)" = 'pie_rank_na_mean'),
                    `Absolute metascores` = c("BayesDel" = 'pie_bayesdel', 
                                              "LINSIGHT" = 'pie_linsight', 
                                              "IW-Scoring Known (K10)" = 'pie_iwscoring_known',
                                              "IW-Scoring Novel (K6)" = 'pie_iwscoring_novel',
                                              "All absolute metascores" = 'pie_all')
                  )
                ),
                selected = 'pie_scores'
    ),
    conditionalPanel(condition="input.snv_nodes_type=='pie_scores' | input.snv_nodes_type=='pie_rank_na_last' | input.snv_nodes_type=='pie_rank_na_mean'",
                     fluidRow(
                       column(width = 12,
                              selectInput("selected_scores", 'Predictors',
                                          choices = c(),
                                          selectize = TRUE, multiple = TRUE))
                     ),
                     actionButton("update_metascore", "Update"),
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
    jqui_resizable(visNetworkOutput("my_network", height = "600px")),
    bsModal(id = "modalExample", title = "Details", trigger = "current_node_id", size = "small",
            DT::dataTableOutput(outputId = "snv_score_details"))
  )
}

ld_mapping_module <- function(){
  list(
    selectInput(inputId = 'population',
                label = 'Choose the population to use',
                choices = populations),
    shinyBS::bsButton(inputId = 'runLD', label = 'Add LD Information',  
                      icon = icon("gear"), disabled = FALSE),
    br(),
    br(),
    verbatimTextOutput("runld_res"),
    tags$style(type="text/css", "#runld_res {white-space: pre-wrap;}")
  )
}

presentation_module <- function(){
  fluidRow(box(id = "intro_box", title = "DSNetwork", solidHeader = FALSE, background = NULL, width = 12,
               collapsible = TRUE, collapsed = TRUE,
               column(width = 12, 
                      helpText("Text de presentation")
               )
  ))
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
#                        shinyBS::bsButton(inputId = 'runLD', label = 'Add LD Information',  
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
#            jqui_resizable(plotlyOutput(outputId = "my_plot",
#                                        height = "400px", width = "auto")),
#            box(title = "Variant selection",
#                width = 12,
#                collapsible = TRUE, 
#                collapsed = FALSE, 
#                htmlOutput("selection")))
#   )
# }