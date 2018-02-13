library(shinydashboard)
library(visNetwork)
library(shinyBS)
library(d3heatmap)

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
                  "Locus 80 [chr3:86537543-8753754]" = "locus_80"
                  )

dashboardPage(
  header = dashboardHeader(title = "DSNetwork"),
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("How to get a network", tabName = "main", icon = icon("hand-o-up")),
      helpText("Help text"),
      menuItem("ReadMe", tabName = "readme", icon=icon("mortar-board")),
      menuItem("About", tabName = "about", icon = icon("question"))
    )
  ),
  body = dashboardBody(
    tabItems(
      tabItem(tabName = "main",
              fluidRow(
                box(title = "Data", status = "primary", width = 12, 
                    fluidRow(
                      column(width = 5,
                             textAreaInput("query", "1) Query", "", rows = 5, 
                                           placeholder = "Please enter one variant id per line (rs123455 or 1:1324:A:C)"),
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
                      ),
                      column(width = 2, 
                             shinyBS::bsButton(inputId = 'buildNetwork', label = 'Build Network',
                                               disabled = TRUE, icon = icon("gear")),
                             br()
                      )
                    )
                )
              ),
              fluidRow(
                tabBox( width = 12,
                        tabsetPanel(id = "results_tabset",
                        tabPanel(title = h5("Raw results"),
                                 value = "raw_results",
                                 conditionalPanel(condition="input.fetch_annotations",
                                                  DT::dataTableOutput("raw_data"),
                                                  downloadButton('downloadRawTable', 'Download')
                                 )
                        ),
                        # tabPanel(h5("Population Frequencies"),
                        #          DT::dataTableOutput("populations"),
                        #          downloadButton('downloadFreqTable', 'Download')
                        # ),
                        tabPanel(title = h5("Network"),
                                 value = "network_results",
                                 conditionalPanel(condition="input.buildNetwork",
                                                  fluidRow(
                                                    column(width = 9,
                                                           box(width = NULL, title = "Network", collapsible = TRUE, 
                                                               visNetworkOutput("my_network",
                                                                                height = "500", width = "auto")
                                                           ), 
                                                           box(width = NULL, status = "info", title = "Scores Stats", height = "auto",
                                                               collapsible = T,
                                                               fluidRow(
                                                                 column(width = 6, DT::dataTableOutput(outputId = "scores_stats")),
                                                                 column(width = 6, d3heatmapOutput(outputId = "scores_corr"))
                                                               )
                                                           )
                                                           # ,
                                                           # box(width = NULL, status = "info", title = "Legends", height = "auto",
                                                           #     collapsible = T, 
                                                           #     fluidRow(
                                                           #       column(width = 12, plotOutput(outputId = "color_key", height = "200"))
                                                           #     )
                                                           # )
                                                    ),
                                                    column(width = 3,
                                                           box(width = NULL, status = "warning",
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
                                                           )
                                                           ,
                                                           box(width = NULL, status = "warning",
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
                                                           # ,
                                                           # box(width = NULL, status = "warning",
                                                           #     checkboxGroupInput("adj_scores", "Adjusted scores",
                                                           #                        choices = c()
                                                           #     ),
                                                           #     checkboxGroupInput("raw_scores", "Raw scores",
                                                           #                        choices = c()
                                                           #     ),
                                                           #     actionButton("update_metascore", "Update"),
                                                           #     p(
                                                           #       class = "text-muted", br(),
                                                           #       paste("This option enables to select the set of prediction and",
                                                           #             "scoring algorithms used to compute the metascore (color of the database-shaped nodes)"
                                                           #       )
                                                           #     )
                                                           # )
                                                    )
                                                  )
                                 )
                        ),
                        tabPanel(title = h5("Scores details"),
                                 value = "scores_details",
                                 conditionalPanel(condition="input.buildNetwork",
                                                  fluidRow(
                                                    column(width = 9,
                                                           box(width = NULL, title = "Network", collapsible = TRUE, 
                                                               visNetworkOutput("my_subnetwork",
                                                                                height = "500", width = "auto")
                                                           ), 
                                                           box(width = NULL, status = "info", title = "Variant Stats", height = "auto",
                                                               collapsible = T,
                                                               fluidRow(
                                                                 column(width = 12, plotOutput(outputId = "critical_dist"))
                                                               )
                                                           )
                                                    ),
                                                    column(width = 3,
                                                           box(width = NULL, status = "warning",
                                                               sliderInput("neg_cor_range", "Negative correlation range",
                                                                           min = -1, max = 0, value = c(-1, -0.8)),
                                                               sliderInput("pos_cor_range", "Positive correlation range",
                                                                           min = 0, max = 1, value = c(0.8, 1)),
                                                               actionButton("update_cor", "Update"),
                                                               p(class = "text-muted",
                                                                 br(),
                                                                 "This option enables to select the interval of correlation values represented between annotations"
                                                               )
                                                           ),
                                                           box(width = NULL, status = "warning",
                                                               selectInput("cor_focus", "Focus on",
                                                                           choices = c(
                                                                             "None" = -1
                                                                           ),
                                                                           selected = -1
                                                               ),
                                                               p(class = "text-muted",
                                                                 br(),
                                                                 "This option enables to focus the network on a particular annotation"
                                                               )
                                                           ),
                                                           box(width = NULL, status = "warning",
                                                               checkboxGroupInput("adj_scores", "Adjusted scores",
                                                                                  choices = c()
                                                               ),
                                                               checkboxGroupInput("raw_scores", "Raw scores",
                                                                                  choices = c()
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
                                                  )
                                                  
                                 )
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