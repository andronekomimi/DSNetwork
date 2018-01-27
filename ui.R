library(shinydashboard)
library(visNetwork)

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

dashboardPage(
  header = dashboardHeader(title = "DSNetwork"),
  sidebar = dashboardSidebar(
    sidebarMenu(
      textAreaInput("query", "Query", "", rows = 5, 
                    placeholder = "Please enter one variant id per line (rs123455 or 1:1324:A:C)"),
      actionButton("transform_query", "Search for annotations", icon = icon("search")),
          selectInput(inputId = 'population',
                      label = 'Choose the population to use',
                      choices = populations),
          shinyBS::bsButton(inputId = 'runLD', label = 'Add LD Information',  
                   icon = icon("gear")),
      menuItem("Data", tabName = "load_data", icon = icon("th")),
      menuItem("Network", tabName = "network", icon = icon("snowflake-o"))
    )
  ),
  body = dashboardBody(
    tabItems(
      tabItem(tabName = "load_data",
              fluidRow(
                box(title = "Infos", status = "warning", width = 4, 
                    solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                    br(),
                    verbatimTextOutput("transform_res"),
                    tags$style(type="text/css", "#transform_res {white-space: pre-wrap;}"),
                    br(),
                    verbatimTextOutput("runld_res"),
                    tags$style(type="text/css", "#runld_res {white-space: pre-wrap;}"),
                    br(),
                    verbatimTextOutput("query_res"),
                    tags$style(type="text/css", "#query_res {white-space: pre-wrap;}")
                ),
                box(title = "Populations data", status = "primary", width = 8, 
                    solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                    br(),
                    DT::dataTableOutput("populations")
                )
              ),
              fluidRow(
                
              ),
              fluidRow(
                box(title = "Network", status = "success", width = 12, 
                    solidHeader = TRUE, collapsible = TRUE, 
                    shinyBS::bsButton(inputId = 'buildNetwork', label = 'Build Network',  
                                      icon = icon("gear")),
                    br(),
                    h1('N  E  T  W  O  R  K'),
                    visNetworkOutput("my_network",
                                     height = "800", width = "auto")
                )
              ),
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
                           checkboxGroupInput("metascore", "Meta-score",
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
                           actionButton("update_metascore", "Update"),
                           p(
                             class = "text-muted",
                             paste("This option enables to select the set of prediction and",
                                   "scoring algorithms used to compute the metascore (color of the database-shaped nodes)"
                             )
                           )
                       ),
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
                           actionButton("update_annotations", "Update"),
                           p(
                             class = "text-muted",
                             paste("This option enables to select the set of prediction and",
                                   "scoring algorithms represented in the network "
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
                       ),
                       box(width = NULL, status = "warning",
                           sliderInput("ld_range", "LD range",
                                       min = 0, max = 1, value = c(0, 1)),
                           actionButton("update_ld", "Update"),
                           p(class = "text-muted",
                             br(),
                             "This option enables to select the interval of LD values represented between variants"
                           )
                       )
                )
              )
      )
    )
  ), 
  title = "DSNetwork",
  skin = "black"
)