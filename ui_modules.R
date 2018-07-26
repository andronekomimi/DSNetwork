library(shinydashboard)
library(visNetwork)
library(shinyBS)
library(d3heatmap)

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

preload_loci <- c("Locus 2 [chr1:113948389-14948389]" = "locus_2",
                  "Locus 6 [chr1:201687176-202687176]" = "locus_6",
                  "Locus 38 [chr7:143574929-144574929]" = "locus_38",
                  "Locus 60 [chr13:32468810-33472626]" = "locus_60",
                  "Locus 70 [chr17:77281387-78281725]" = "locus_70",
                  "Locus 78 [chr22:40376234-41527870]" = "locus_78",
                  "Locus 80 [chr3:86537543-8753754]" = "locus_80")

sidebar_content <- function(){
  sidebarMenu(
    menuItem("How to get a network", tabName = "main", icon = icon("hand-o-up")),
    helpText("Help text"),
    menuItem("ReadMe", tabName = "readme", icon=icon("mortar-board")),
    menuItem("About", tabName = "about", icon = icon("question"))
  )
}



input_data_module <- function(){
  tabPanel(title = h5("INPUT"),
           value = "input",
           fluidRow(
             column(width = 4,
                    br(),
                    textAreaInput("query", "1) Enter variant ids", "", rows = 5, 
                                  placeholder = "Please enter one variant id per line (rs123455 or 1:1324:A:C)"),
                    fileInput("query_file", "or load text file (one variant id per line)", 
                              multiple = FALSE, 
                              accept = c(
                                "text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv")
                    ),
                    selectInput(inputId = 'preload',
                                label = '0) Load preset query',
                                choices = preload_loci),
                    actionButton("preload_loci", "Load preset query", icon = icon("caret-right")),
                    shinyBS::bsButton(inputId = "fetch_annotations", 
                                      label = "Fetch Annotations", 
                                      icon = icon("search"), disabled = TRUE),
                    verbatimTextOutput("transform_res"),
                    tags$style(type="text/css", "#transform_res {white-space: pre-wrap;}"),
                    br(),
                    verbatimTextOutput("query_res"),
                    tags$style(type="text/css", "#query_res {white-space: pre-wrap;}"),
                    br()
             ), 
             output_raw_results_module()
           )
  )
}

input_ld_module <- function(){
  column(width = 5,
         selectInput(inputId = 'population',
                     label = '2) LD data : choose the population to use',
                     choices = populations),
         shinyBS::bsButton(inputId = 'runLD', label = 'Add LD Information',  
                           icon = icon("gear"), disabled = TRUE),
         br(),
         br(),
         verbatimTextOutput("runld_res"),
         tags$style(type="text/css", "#runld_res {white-space: pre-wrap;}"),
         br()
  )
}

input_network_module <- function(){
  column(width = 2, 
         tags$b('3) Build your network'),br(),
         shinyBS::bsButton(inputId = 'buildNetwork', label = 'Build Network',
                           disabled = TRUE, icon = icon("gear")),
         br()
  )
}

output_raw_results_module <- function(){
  # tabPanel(title = h5("Raw results"),
  #          value = "raw_results",
  #          br(),
  column(width = 8, 
           br(),
           conditionalPanel(condition="input.fetch_annotations",
                            #DT::dataTableOutput("raw_data"),
                            fluidRow(column(width = 6, C3PieChartOutput(outputId = "PieBranch1")),  column(width = 6, C3PieChartOutput(outputId = "PieBranch2"))),
                            fluidRow(column(width = 6, C3PieChartOutput(outputId = "PieBranch3")),  column(width = 6, C3PieChartOutput(outputId = "PieBranch4"))),
                            downloadButton('downloadRawTable', 'Download raw results (csv)')
                            
           )
  )
} 

output_ld_results_module <- function(){
  tabPanel(title = h5("LD Plots"),
           value = "ld_results",
           conditionalPanel(condition="input.runLD",
                            selectInput("ld_regions", "LD Regions",
                                        width = "20%",
                                        choices = c()
                            ),
                            imageOutput("ld_plot", height = 600)
           )
  )
}

output_network_row <- function(){
  fluidRow(
    column(width = 9,
           fluidRow(
             column(width = 6,
                    box(width = NULL, status = "success",
                        collapsible = TRUE, height = "600px",
                        visNetworkOutput("my_network",
                                         height = "550px", width = "auto")
                    )),
             column(width = 6,
                    box(width = NULL, status = "success",
                        collapsible = TRUE, height = "600px",
                        visNetworkOutput("my_network_2",
                                         height = "550px", width = "auto")
                    )))
    ),
    column(width = 3,
           box(width = NULL, status = "success", height = "600px",
               textOutput("snv_score_details_id"),
               DT::dataTableOutput(outputId = "snv_score_details")
           )
    )
  )
} 

nodes_modifiers_box <- function(){
  box(width = NULL, status = "warning", 
      height = "500px", collapsible = TRUE,
      selectInput("snv_nodes_type", "SNV Nodes color",
                  choices = c(
                    "Scores Pie" = 'pie_scores',
                    "Metascore" = 'pie_metascores', 
                    "Rank (NA last)" = 'pie_rank_na_last',
                    "Rank (NA mean)" = 'pie_rank_na_mean'
                  ),
                  selected = 'pie_scores'
      ),
      #),
      #box(width = NULL, status = "warning",
      selectInput("focus", "Focus on",
                  choices = c(
                    "None" = -1
                  ),
                  selected = -1
      ),
      p(class = "text-muted",
        br(),
        "This option enables to focus the network on a particular variant"
      ),
      #),
      #box(width = NULL, status = "warning",
      selectInput("highlight", "Highlight variants",
                  choices = c(
                    "None"
                  ),
                  selected = "None"
      ),
      p(class = "text-muted",
        br(),
        paste("This option enables to highlight particular variants according to their type.",
              " Replacing the lightblue circles by red stars.")
      )
  )
}

edges_modifiers_box <- function(){
  box(width = NULL, status = "warning",  
      height = "500px", collapsible = TRUE,
      selectInput("snv_edges_type", "SNV Edges",
                  choices = c(
                    "Distance" = 0,
                    "Linkage" = 1
                  ),
                  selected = 1
      ),
      conditionalPanel(condition="input.snv_edges_type=='0'",
                       sliderInput("dist_range", "Distance (kb)",
                                   min = 0, max = 1000, step = 1, value = 1000),
                       actionButton("update_dist", "Update"),
                       p(class = "text-muted",
                         br(),
                         "This option enables to select the interval of 
                      distance represented between variants"
                       )
      ),
      conditionalPanel(condition="input.snv_edges_type=='1'",
                       sliderInput("ld_range", "LD range",
                                   min = 0, max = 1, value = c(0.8, 1)),
                       actionButton("update_ld", "Update"),
                       p(class = "text-muted",
                         br(),
                         "This option enables to select the interval of LD values represented between variants"
                       )
      )
  )
}

predictors_modifiers_box <- function(){
  box(width = NULL, status = "warning",
      height = "500px", collapsible = TRUE,
      fluidRow(
        column(width = 6,
               selectInput("adj_scores", 'Ranked scores', 
                           choices = c(), 
                           selectize = TRUE, multiple = TRUE)),
        column(width = 6,
               selectInput("raw_scores", "Raw scores",
                           choices = c(), 
                           selectize = TRUE, multiple = TRUE
               ))
      ),
      fluidRow(
        column(width = 12,
               selectInput("non_na_scores", 'Non-NA scores', 
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
}


network_modifiers_row <- function(){
  fluidRow(
    column(width = 3,
           nodes_modifiers_box()
    ),
    column(width = 3, 
           edges_modifiers_box()   
    ),
    column(width = 6,
           predictors_modifiers_box()
    )
  )
}


output_network_results_modules <- function(){
  tabPanel(title = h5("Network"),
           value = "network_results",
           conditionalPanel(condition="input.buildNetwork",
                            output_network_row(),
                            network_modifiers_row()
           )
  )
}


predictors_selection <- function(){
  fluidRow(
    column(width = 12,
           selectInput("predictors", 'Predictors',width = "50%",
                       choices = c(), selectize = TRUE, multiple = TRUE),
           bsButton(inputId = 'get_predictors_info', label = "Get info",
                    icon = icon("search")),
           br(),br()
    )
  )
}

output_predictors_descript <- function(){
  fluidRow(
    column(width = 6,
           box(width = NULL, status = "success",
               collapsible = TRUE,
               plotOutput(outputId = "scores_stats",height = "500px")
           ),
           box(width = NULL, status = "success",
               d3heatmapOutput(outputId = "scores_corr", height = "600px")
           )
    ),
    column(width = 6,
           box(width = NULL, status = "success",
               DT::dataTableOutput(outputId = "scores_missing_data")
           )
    )
  )
}

output_predictors_results_modules <- function(){
  tabPanel(title = h5("Scores information"),
           value = "scores_stats",
           conditionalPanel(condition="input.buildNetwork",
                            predictors_selection(),
                            output_predictors_descript()
           )
  )
}



