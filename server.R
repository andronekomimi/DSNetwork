require(shiny)
#require(d3heatmap)
require(visNetwork)
require(myvariant)
require(shinyBS)
require(plotly) # do not update to 4.8.0, keep 4.7.1
require(shinyjs)
#require(tableHTML)
#require(magrittr)
require(shinyalert)

options(shiny.trace = FALSE)

server <- function(input, output, session) {
  
  tmpDir <- tempdir() 
  print(tmpDir)
  #dir.create(path = paste0(tmpDir, "ld_figures"), showWarnings = F)
  #dir.create(path = paste0(tmpDir, "objects"), showWarnings = F)
  
  source('helper.R', local = TRUE)
  
  if(is.local()){
    #### CONF ####
    appDir <- "/home/nekomimi/Workspace/DSNetwork/"
    dataDir <- "/home/nekomimi/Workspace/DSNetwork/data/"
    path_to_victor <- paste0(appDir, "softs/VICTOR/")
    python_path <- "/usr/bin/python"
    tabix_path <- '/usr/bin/tabix'
    path_to_vcf_converter <- paste0(appDir, "scripts/vcf_to_ped_converter.pl")
    load(paste0(appDir, "demo/preload.rda"))
  } else {
    #### CONF ####
    appDir <- "/srv/shiny-server/dsnetwork/"
    dataDir <- "/mnt/apps_data/dsnetwork/"
    path_to_victor <- "/mnt/apps_softs/dsnetwork/VICTOR/"
    python_path <- "/bin/python" ###
    tabix_path <- '/mnt/apps_softs/dsnetwork/TABIX/tabix' ###
    path_to_vcf_converter <- paste0(appDir, "scripts/vcf_to_ped_converter.pl")
    load(paste0(appDir, "demo/preload.rda"))
  }
  
  app.conf <- list(TABIX = tabix_path,
                   VCF =  paste0(dataDir, '1000Genomes/')) 
  path_to_images <- paste0(appDir, 'www/scores_figures/')
  
  values <- reactiveValues()
  values$annotations <- as.matrix(data.frame(waiting = ""))
  values$all_edges <- NULL
  values$all_nodes <- NULL
  values$current_edges <- NULL
  values$current_nodes <- NULL
  values$scores_data <- NULL
  values$selected_node <- ''
  values$ld_regions <- NULL
  values$notfound_id <- NULL
  values$can_run <- FALSE
  values$all_scores_data <- NULL
  
  # A notification ID
  id <- NULL
  
  cat("Session start\n")
  
  observeEvent(c(input$query_file, input$query), ({
    
    if(!is.null(input$query_file) && input$query_file$size > 0){
      updateButton(session = session, inputId = "fetch_annotations",
                   disabled = FALSE)
    } else {
      if(nchar(input$query) > 0){
        updateButton(session = session, inputId = "fetch_annotations",
                     disabled = FALSE)
      } else {
        updateButton(session = session, inputId = "fetch_annotations",
                     disabled = TRUE)
      }
    }
  }))
  
  #### LOAD PRESET QUERY ####
  observeEvent(input$preload_loci, {
    updateTextAreaInput(session = session, inputId = "query", value = preload[[input$preload]])
  })
  
  #### LOAD DEMO QUERY ####
  observeEvent(input$load_demo1, {
    updateTextAreaInput(session = session, inputId = "query", value = preload[["locus_1p36"]])
  })
  
  observeEvent(input$load_demo2, {
    updateTextAreaInput(session = session, inputId = "query", value = preload[["locus_1p34"]])
  })
  
  observeEvent(input$load_demo3, {
    updateTextAreaInput(session = session, inputId = "query", value = preload[["locus_7q22"]])
  })
  
  observeEvent(input$load_demo4, {
    updateTextAreaInput(session = session, inputId = "query", value = preload[["locus_11p15"]])
  })
  
  #### INPUT DATA MANAGEMENT ####
  transform_query <- function(string){
    if(nchar(input$query) > 0){
      string <- unlist(strsplit(x = string, split = "\n"))
      string <- gsub(x = string, pattern = " +$", replacement = "")
      
      modstring <- sapply(X = string,  FUN = function(x) {
        #rsids
        if(grepl(pattern = "^rs\\d+$", x = x)){
          return(x)
        }  else { 
          x <- unlist(strsplit(x = x, split = ":"))
          if(length(x) != 4){
            return('FAIL')
          } else {
            if(suppressWarnings(!is.na(as.numeric(x[2])))){
              myvariant::formatSingleHgvs(x[1], 
                                          as.numeric(x[2]), 
                                          toupper(x[3]), 
                                          toupper(x[4]))
            } else {
              return('FAIL')
            }
          }
        }
      })
      return(modstring)
    }
  }
  
  transform_query_file <- function(filepath){
    if(file.exists(filepath)){
      string <- read.csv(file = filepath, header = F, stringsAsFactors = F)$V1
      string <- tolower(gsub(x = string, pattern = " +$", replacement = ""))
      
      modstring <- sapply(X = string,  FUN = function(x) {
        # rsids
        if(grepl(pattern = "^rs\\d+$", x = x)){
          return(x)
        }  else { 
          x <- unlist(strsplit(x = x, split = ":"))
          if(length(x) != 4){
            return('FAIL')
          } else {
            if(suppressWarnings(!is.na(as.numeric(x[2])))){
              myvariant::formatSingleHgvs(x[1], 
                                          as.numeric(x[2]), 
                                          toupper(x[3]), 
                                          toupper(x[4]))
            } else {
              return('FAIL')
            }
          }
        }
      })
      return(modstring)
    }
  }
  
  #### FETCH ANNOTATIONS ####
  
  # observeEvent(input$fetch_annotations, {
  #   js$collapse("intro_box")
  #   js$collapse("selection_box")
  #   js$collapse("network_box")
  # })
  
  observeEvent(input$fetch_annotations, {
    # fermeture de tous les graphiques existants
    graphics.off()
    
    updateButton(session = session, inputId = "fetch_annotations", 
                 disabled = TRUE)
    
    withProgress(message = 'Fetching annotations', value = 0, {
      n <- 10
      modstring <- c()
      if(!is.null(input$query_file))
        modstring <- c(modstring, transform_query_file(input$query_file$datapath))
      
      if(nchar(input$query) > 0)
        modstring <- c(modstring, transform_query(tolower(input$query)))
      
      valid_transfo <- unique(modstring[!grepl(x = modstring, pattern = 'FAIL')])
      
      save(modstring, file = paste0(tmpDir, '/modstring.rda'))
      
      fail_transfo <- names(modstring[grep(x = modstring, pattern = 'FAIL')])
      if(length(modstring) > 0){
        if(length(fail_transfo) > 0){
          res0 <- paste0('Id recognition fails for: ', paste(fail_transfo, collapse = ","), '.')
          print(res0)
          createAlert(session = session, "alert_conv", "alert1", 
                      title = "Id recognition",
                      content = res0, append = TRUE)
        }
      }
      
      if(length(valid_transfo) > 0){
        #### run myvariant ####
        incProgress(1/n, detail = "Interrogating MyVariant.info...")
        res <- as.data.frame(myvariant::getVariants(hgvsids = valid_transfo,
                                                    verbose = F, return.as = "DataFrame"))
        
        
        
        #### add metascores columns ####
        res$linsight <- res$cdts_score <- res$cdts_percentile <- NA
        res$bayesdel <- res$iwscoring_known <- res$iwscoring_novel <- NA
        
        values$res <- res
        remove(res)
        if(!is.null(values$res) && nrow(values$res) > 0){
          values$can_run <- TRUE
          #### renommage des doublons ####
          duplicated_entries <- unique(values$res$query[duplicated(values$res$query)])
          row_2_delete <- c() 
          
          if(length(duplicated_entries) > 0){
            for(duplicated_entry in duplicated_entries){
              values$res[values$res$query == duplicated_entry,]$query <- paste0(values$res[values$res$query == duplicated_entry,]$query,"_",
                                                                                values$res[values$res$query == duplicated_entry,]$dbsnp.ref,".",
                                                                                values$res[values$res$query == duplicated_entry,]$dbsnp.alt)
            }
          }
          
          #### order ####
          values$res <- values$res[order(values$res$dbsnp.chrom, values$res$dbsnp.hg19.start),]
          
          
          #### global range ####
          global_ranges <- apply(X = values$res, MARGIN = 1, 
                                 FUN = function(x) hgvsToGRange(hgvs_id = x[names(x) == "X_id"]$X_id, 
                                                                query_id = x[names(x) == "query"]$query))
          names(global_ranges) <- NULL
          global_ranges <- do.call("c", global_ranges) 
          global_ranges <- sort(global_ranges)
          
          save(global_ranges, file = paste0(tmpDir, "/global_ranges.rda"))
          values$global_ranges <- global_ranges
          
          ### remove nonfound variant lines 
          if("notfound" %in% colnames(values$res)){
            values$notfound_id <- paste(values$res[!is.na(values$res$notfound),"query"], collapse = ", ")
            values$res <- values$res[-which(values$res$notfound == TRUE), ]
          }
          
          
          requested_chromosomes <- GenomeInfoDb::seqlevelsInUse(global_ranges)
          
          #### fetch linsight ####
          incProgress(1/n, detail = "Extracting LINSIGHT scores...")
          for(requested_chr in requested_chromosomes){
            if(file.exists(paste0(dataDir,'LINSIGHT/LINSIGHT_',requested_chr,'.rda'))){
              load(paste0(dataDir,'LINSIGHT/LINSIGHT_',requested_chr,'.rda'))
              hits <- findOverlaps(query = values$global_ranges, subject = gr)
              for(i in seq_along(hits)){ 
                hit <- hits[i]
                query <- values$global_ranges[queryHits(hit)]$query
                value <- gr[subjectHits(hit)]$score
                values$res[values$res$query == query,"linsight"] <- value
              }
            }
          }
          
          #### fetch cdts ####
          incProgress(1/n, detail = "Extracting CDTS data...")
          for(requested_chr in requested_chromosomes){
            if(file.exists(paste0(dataDir,'CDTS/CDTS_hg19/CDTS_',requested_chr,'.rda'))){
              load(paste0(dataDir,'CDTS/CDTS_hg19/CDTS_',requested_chr,'.rda'))
              hits <- findOverlaps(query = values$global_ranges, subject = CDTS)
              for(i in seq_along(hits)){ 
                hit <- hits[i]
                query <- global_ranges[queryHits(hit)]$query
                value1 <- CDTS[subjectHits(hit)]$CDTS
                value2 <- CDTS[subjectHits(hit)]$percentile
                values$res[values$res$query == query,"cdts_score"] <- value1
                values$res[values$res$query == query,"cdts_percentile"] <- value2
              }
            }
          }
          
          #### store CDTS scores on the complete region
          cdts_hits <- findOverlaps(query = range(values$global_ranges), subject = CDTS, maxgap = 100)
          if(length(cdts_hits) > 0){
            cdts_hits <- as.data.frame(CDTS[subjectHits(cdts_hits)])
            values$cdts_region <- cdts_hits
          } else {
            values$cdts_region <- NULL
          }
          
          #### fetch bayesdel
          incProgress(1/n, detail = "Extracting BayesDel scores...")
          filename <- tempfile(tmpdir = tmpDir, fileext = ".vcf")
          print(filename)
          createVCF(session_values = values, filename = filename)
          value <- extractBayesDel(path_to_victor, filename)
          
          if(!is.null(value)){
            values$res$bayesdel <- value
          }
          
          #### fetch SNPnexus from HTML request ####
          if(input$fetch_snpnexus){
            t1 <- Sys.time()
            incProgress(4/n, detail = "Interrogating SNPNexus platform...")
            path_to_snpnexus <- paste0(appDir, "scripts/")
            
            # create temp snps file
            filename <- tempfile(tmpdir = tmpDir, fileext = ".txt")
            print(filename)
            save(values, file = paste0(tmpDir, "/res.rda"))
            createSNPnexusInput(session_values = values, filename = filename)
            value <- runSNPnexus(python_path = python_path,
                                 path_to_snpnexus = path_to_snpnexus, 
                                 filename = filename, 
                                 waiting_time = input$waiting)
            print(Sys.time() - t1)
            save(value, file = paste0(tmpDir, "/value.rda"))
            
            if(!is.null(value)){
              pre_res <- values$res
              save(pre_res, file = paste0(tmpDir, "/pre_res.rda"))
              #`id,cadd_phred,deepsea,eigen,eigen_pc,fathmm,fitcons,funseq2,gwava_region,gwava_tss,gwava_unmatched,remm,iwscorek11,pvalk11,iwscorek10,pvalk10`
              x1 <- value[[1]] 
              #`id,cadd_phred,deepsea,eigen,eigen_pc,fathmm,fitcons,funseq2,remm,iwscoren8,pvaln8,iwscoren6,pvaln6`
              x2 <- value[[2]]
              snpnexus_res <- merge(x = x1, y = x2, by.x = "id", by.y = "id", all = T)
              
              #convert id to submit to myvariant info to get a common id between snpnexus and res
              x3 <- sapply(X = snpnexus_res$id, FUN = function(x) snpnexusIDconversion(x))
              # deletion mal convertie
              for(j in which(x3[2,] == "deletion")){
                #chr9:g.110893722_110893725del => #chr9:g.110893721_110893724del to fit with myvariant.info format
                chrom_part <- gsub(x = x3[1,j], pattern = "(chr.*:g.).*", replacement = "\\1")
                start_part <- as.numeric(gsub(x = x3[1,j], pattern = "chr.*:g.(\\d+)_.*", replacement = "\\1")) - 1
                stop_part <- as.numeric(gsub(x = x3[1,j], pattern = ".*_(\\d+)del", replacement = "\\1")) - 1
                if(start_part == stop_part){
                  x3[1,j] <- paste0(chrom_part, start_part,"del")
                } else {
                  x3[1,j] <- paste0(chrom_part, start_part, "_", stop_part,"del")
                }
              }
              converted_id <- as.character(x3[1,])
              
              snpnexus_res$converted_id <- converted_id
              # merge results based on common id
              pre_res <- merge(y = snpnexus_res, x = pre_res, by.y = "converted_id", by.x = "X_id", all = T)
              # metascores
              pre_res$iwscoring_novel <- pre_res$iwscoren6
              pre_res$iwscoring_known <- pre_res$iwscorek10
              # non-coding scores
              colnames(pre_res)[colnames(pre_res) == "deepsea.x"] <- "deepseq_sig_log2" #deppsea : Lower score indicates higher likelihood of functional significance of the variant. But log2 is in the "good" direction
              colnames(pre_res)[colnames(pre_res) == "eigen.x"] <- "eigen_nc"
              colnames(pre_res)[colnames(pre_res) == "eigen_pc.x"] <- "eigen_pc_nc" #higher score more likelihood of the variant to be functional. 
              colnames(pre_res)[colnames(pre_res) == "fathmm.x"] <- "fathmm_nc" #Scores above 0.5 are predicted to be deleterious
              colnames(pre_res)[colnames(pre_res) == "fitcons.x"] <- "fitcons_nc" #with higher scores indicating more potential
              colnames(pre_res)[colnames(pre_res) == "funseq2.x"] <- "funseq" #higher scores indicating more likely to be functional.
              colnames(pre_res)[colnames(pre_res) == "remm.x"] <- "remm" #higher scores indicating more likely to be deleterious.
              #gwava : higher scores indicating variants predicted as more likely to be functional.
              
              # suppress unsued column
              pre_res <- pre_res[,!grepl(x = colnames(pre_res), pattern = "\\.x$|\\.y$|^iwscore|converted_id|pval")]
              
              # re-insert res table
              values$res <- pre_res
            } else {
              print("SNPnexus taking too long to answer, sorry...")
              createAlert(session = session, anchorId = "alert_res",
                          alertId = "alert4", title = "Data retrieval",
                          content = "SNPnexus taking too long to answer, sorry..." ,
                          append = TRUE, style = "warning")
            }
            print(Sys.time() - t1)
          } else {
            incProgress(4/n, detail = "Skipping SNPNexus...")
          }
          
          incProgress(2/n, detail = "Aggregating data...")
          
          common_fields <-c("query","X_id")
          
          included_scores <- read.csv(file = paste0(appDir, '/scores_description.tsv'), header = T, sep = "\t", stringsAsFactors = F)
          
          #### convert scores (toNumeric and mean if needed)
          values$all_res <- values$res[,(colnames(values$res) %in% c(common_fields, included_scores$id))]
          
          if(nrow(values$all_res) > 0){
            
            #### suppress empty scores
            values$all_res <- values$all_res[,sapply(X = values$all_res, FUN = function(x) sum(is.na(x)) != length(x))]
            
            adjusted_scores <- sapply(X = c(common_fields, included_scores$id), 
                                      FUN = function(x) extract_score_and_convert(values$all_res, score_name = x))
            
            switch (class(adjusted_scores),
                    matrix = {
                      rownames(adjusted_scores) <- NULL
                      adjusted_scores <- as.list(data.frame(adjusted_scores))
                    },
                    numeric = {
                      names(adjusted_scores) <- gsub(x = names(adjusted_scores), pattern = "(.*[a-z]).\\d+.?\\d+$", replacement = "\\1")
                      adjusted_scores <- lapply(split(adjusted_scores, names(adjusted_scores)), unname)
                    }
            )
            
            adjusted_scores <- adjusted_scores[sapply(X = adjusted_scores, FUN = function(x) sum(is.na(x)) != length(x))] 
            
            for(s in names(adjusted_scores)){
              values$all_res[[s]] <- adjusted_scores[[s]]
            }
            
          } else {
            adjusted_scores <- NULL
          }
          
          save(adjusted_scores, file = paste0(tmpDir, "/scores.rda"))
          values$adjusted_scores <- names(adjusted_scores)
          
          #### split res in 2 -> non-syn vs other
          non_syn_res <- sapply(values$res$cadd.consequence, function(x) "NON_SYNONYMOUS" %in% x)
          values$all_nonsyn_res <- values$all_res[non_syn_res, ]
          values$all_regul_res <- values$all_res[!non_syn_res, ]
          
          #### add the overall ranking  
          temp_df <- NULL
          if(nrow(values$all_nonsyn_res) > 0){
            overall_mean_rank <- basic_ranking(score_dataframe = values$all_nonsyn_res[,values$adjusted_scores])
            values$all_nonsyn_res$overall_mean_rank_na_last <- overall_mean_rank$na_last
            values$all_nonsyn_res$overall_mean_rank_na_mean <- overall_mean_rank$na_mean
            
            temp_df <- values$all_nonsyn_res[,c("X_id",
                                                "overall_mean_rank_na_last",
                                                "overall_mean_rank_na_mean")]
          }
          
          if(nrow(values$all_regul_res) > 0){
            overall_mean_rank <- basic_ranking(score_dataframe = values$all_regul_res[,values$adjusted_scores])
            values$all_regul_res$overall_mean_rank_na_last <- overall_mean_rank$na_last
            values$all_regul_res$overall_mean_rank_na_mean <- overall_mean_rank$na_mean
            
            if(is.null(temp_df)){
              temp_df <- values$all_regul_res[,c("X_id",
                                                 "overall_mean_rank_na_last",
                                                 "overall_mean_rank_na_mean")]
            } else {
              temp_df <- rbind(temp_df, values$all_regul_res[,c("X_id",
                                                                "overall_mean_rank_na_last",
                                                                "overall_mean_rank_na_mean")])
            }
          }
          
          values$res <- merge(x = values$res, y = temp_df, by = "X_id", all = T)
          
          res <- values$res
          save(res, file = paste0(tmpDir, '/res.rda'))
          remove(temp_df)
          
          if(is.null(values$notfound_id)){
            print("Annotations found for all variants")
            createAlert(session = session, anchorId = "alert_res",
                        alertId = "alert2", title = "Data retrieval",
                        content = "Annotations found for all variants" , append = TRUE, style = "success")
          } else {
            print( paste0("No Annotations for the following variants: ",values$notfound_id,"."))
            createAlert(session = session, anchorId = "alert_res",
                        alertId = "alert2", title = "Data retrieval",
                        content = paste0("No Annotations for the following variants: ",
                                         values$notfound_id,".") , append = TRUE, style = "warning")
          }
        } else {
          values$can_run <- FALSE
          print("ERROR : No Annotations found!")
          createAlert(session = session, anchorId = "alert_res",
                      alertId = "alert3", title = "Data retrieval",
                      content = "ERROR : No Annotations found!",
                      append = TRUE, style = "danger")
        }
      } else {
        values$can_run <- FALSE
        print("ERROR : No valid entry!")
        createAlert(session = session, anchorId = "alert_conv",
                    alertId = "alert1", title = "Id recognition",
                    content = "ERROR : No valid entry!",
                    append = TRUE, style = "danger")
      }
      
      if(values$can_run){
        local({
          #### DISPLAY RESULTS TABLE FOR VARIANTS SELECTION ####
          annotations_fields <- c("query","X_id", "cadd.consequence", 
                                  "cdts_score","overall_mean_rank_na_last",
                                  "overall_mean_rank_na_mean")
          
          annotations_infos <- values$res[,annotations_fields]
          annotations_infos$cadd.consequence <- sapply(annotations_infos$cadd.consequence, 
                                                       FUN = function(x) paste(x, collapse = ","))
          # annotations_infos$more <- shinyInput(actionButton, nrow(annotations_infos), 
          #                                      'button_', label = "Fire", 
          #                                      onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
          values$annotations <- annotations_infos
          save(annotations_infos, file = paste0(tmpDir, "/annotations.rda"))
          
          output$raw_data <- DT::renderDataTable({
            
            if(input$network_type == "regul"){
              dt <- values$annotations[!non_syn_res,]
              max_var <- min(sum(!non_syn_res), MAX_VAR)
            } else {
              dt <- values$annotations[non_syn_res,]
              max_var <- min(sum(non_syn_res), MAX_VAR)
            }
            
            n <- order(dt$overall_mean_rank_na_last)[1:max_var]
            
            # round scores
            for (i in 1:ncol(dt)) { if(is.numeric(dt[,i])) {dt[,i] <- round(x = dt[,i], 3)} }
            
            DT::datatable(dt, 
                          escape = FALSE,
                          rownames = FALSE,
                          extensions = c('Scroller','Buttons'),
                          selection = list(mode = 'multiple', selected = n),
                          options = list(dom = 'Bfrtip',
                                         buttons = list(
                                           list(
                                             extend = 'copy',
                                             text = 'Copy',
                                             title = paste0('DSNetwork_',format(Sys.time(), "%m%d%y%H%M%S"))
                                           ), 
                                           list(
                                             extend = 'print',
                                             text = 'Print',
                                             title = paste0('DSNetwork_',format(Sys.time(), "%m%d%y%H%M%S"))
                                           ),
                                           list(
                                             extend = 'csv',
                                             filename = paste0('DSNetwork_',format(Sys.time(), "%m%d%y%H%M%S")),
                                             text = 'Download CSV'
                                           ),
                                           list(
                                             extend = 'excel',
                                             filename = paste0('DSNetwork_',format(Sys.time(), "%m%d%y%H%M%S")),
                                             text = 'Download XLSX'
                                           ),
                                           list(
                                             extend = 'pdf',
                                             title = paste0('DSNetwork_',format(Sys.time(), "%m%d%y%H%M%S")),
                                             text = 'Download PDF',
                                             orientation = 'landscape'
                                           )
                                         ),
                                         scrollX = TRUE,
                                         ordering = TRUE,
                                         deferRender = TRUE,
                                         scrollY = 400,
                                         scroller = TRUE)
            )
          })
          
          output$downloadRawTable <- downloadHandler(
            filename = "annotations_table.tsv",
            content = function(file) {
              readr::write_delim(x = values$res, path = file, delim = "\t")
            }
          )
        })
        
        available_network_types <- c()
        if(sum(!non_syn_res) > 0) available_network_types <- c(available_network_types, 
                                                               "synonymous and non-coding variants" = "regul")
        if(sum(non_syn_res) > 0) available_network_types <- c(available_network_types, 
                                                              "non synonymous variants" = "non_syn")
        
        updateSelectInput(session = session, inputId = "network_type", 
                          choices = available_network_types)
        
        #### BUILDING DATA FOR FIRST DEFAULT PLOT ####
        my_data <- data.frame(values$global_ranges, stringsAsFactors = F)
        my_data$cdts_score <- values$res$cdts_score
        my_data$non_synonymous <- sapply(values$res$cadd.consequence, 
                                         FUN = function(x) "NON_SYNONYMOUS" %in% x)
        my_data$consequences <- sapply(values$res$cadd.consequence, 
                                       FUN = function(x) paste(x, collapse = ","))
        my_data$no_cdts <- is.na(my_data$cdts_score)
        my_data$selected <- TRUE
        
        if(input$network_type == "regul"){
          if(sum(my_data$non_synonymous) > 0){
            my_data[my_data$non_synonymous,]$selected <- FALSE
          }
        } else {
          if(sum(!my_data$non_synonymous) > 0){
            my_data[!my_data$non_synonymous,]$selected <- FALSE
          }
        }
        
        if(nrow(my_data) > MAX_VAR){
          my_data[(MAX_VAR+1):nrow(my_data),]$selected <- FALSE
        }
        
        if(sum(is.na(my_data$cdts_score)) > 0){
          my_data[is.na(my_data$cdts_score),]$cdts_score <- 0
        }
        
        # reste inchangé
        my_data$shape <- "x"
        my_data$shape[my_data$non_synonymous] <- "o"
        my_data$size <- 10
        
        # sera MAJ
        my_data$color <- "blue"
        my_data$color[my_data$selected] <- "red"
        
        values$my_data <- my_data
        save(my_data, file = paste0(tmpDir, "/my_data.rda"))
      }
    })
    
    if(!is.null(values$res) && nrow(values$res) > 0){
      updateButton(session = session, inputId = "buildNetwork", 
                   disabled = FALSE, style = "success")
    } else {
      updateButton(session = session, inputId = "buildNetwork",
                   disabled = TRUE, style = "primary")
    }
  })
  
  output$query_res <- renderText({ 
    if(!is.null(values$res) && nrow(values$res) > 0){
      ifelse(test = is.null(values$notfound_id), 
             no = paste0("No Annotations for the following variants: ", 
                         values$notfound_id,"."), 
             yes = "Annotations found for all variants")  
    }
  })
  
  #### FETCH LD INFOS ####
  observeEvent(input$runLD,{
    my_res <- values$my_res
    if(is.null(my_res))
      return(NULL)
    
    local({
      output$ld_scale <- renderPlot({
        draw_ld_palette()
      })
    })
    
    updateSliderInput(session = session, inputId = "ld_range", value = c(0,1))
    
    id <<- showNotification(paste("Computing linkage disequilibrium..."), duration = 0, type = "message")
    
    
    if(nrow(my_res) > 1){
      my_ranges <- apply(X = my_res, MARGIN = 1, 
                         FUN = function(x) hgvsToGRange(hgvs_id = x[names(x) == "X_id"], 
                                                        query_id = x[names(x) == "query"]))
      names(x = my_ranges) <- NULL
      my_ranges <- do.call("c", my_ranges) 
      my_ranges <- sort(my_ranges)
      
    } else {
      my_ranges <- hgvsToGRange(hgvs_id = my_res$X_id, 
                                query_id = my_res$query)
    }
    
    # build ilets
    hits <- findOverlaps(query = my_ranges, subject = my_ranges, maxgap = 1000000) #1Mbp
    ilets <- c()
    
    for(i in queryHits(hits)){
      ilets <- c(ilets, paste(subjectHits(hits[queryHits(hits) == i]), collapse = "_"))
    }
    
    ilets <- unique(ilets)
    
    ld_results <- sapply(ilets, function(selection){
      current_range  <- setStudyRange(my_ranges, selection)
      ld <- matrix(current_range$query) # notfound set to all
      
      if(length(current_range) > 1){
        region <- setStudyRegion(current_range)
        tryCatch({
          ld <- computeLDHeatmap(region = region,
                                 requested_variants = current_range$query, 
                                 results_dir = tmpDir, 
                                 vcf_dir = app.conf$VCF, 
                                 tabix_path = app.conf$TABIX, 
                                 population = input$population, 
                                 tempDir = tmpDir, i = selection)
        }, error = function(e) {
          print(e)
        })
      }
      return(ld)
    })
    
    save(ld_results, file = paste0(tmpDir, "/ld_results.rda"))
    values$ld <- ld_results
    
    notfound <- do.call("c", values$ld[1,])
    
    #### map LD edges ####
    snv_edges <- build_snv_edges(values, "1", NULL, network_type = input$network_type)
    values$all_edges <- snv_edges
    values$current_edges <- snv_edges
    
    visNetworkProxy("my_network") %>%
      visUpdateEdges(edges = snv_edges)
    
    if(length(notfound) > 0){
      print(paste0('No LD data for the following variants: ', paste(notfound, collapse = ', ')))
      createAlert(session = session, anchorId = "alert_ld",
                  alertId = "alert5", title = "Data retrieval",
                  content = paste0('No LD data for the following variants: ',
                                   paste(notfound, collapse = ', ')) ,
                  append = TRUE, style = "warning")
    } else {
      print("LD computation succeed for all variants")
      createAlert(session = session, anchorId = "alert_ld",
                  alertId = "alert5", title = "Data retrieval",
                  content = "LD computation succeed for all variants",
                  append = TRUE, style = "success")
    }
    
    if(!is.null(ld_results)){
      updateButton(session = session, inputId = "update_ld", 
                   disabled = FALSE)
    } else {
      updateButton(session = session, inputId = "update_ld", 
                   disabled = TRUE)
    } 
    
    if (!is.null(id))
      removeNotification(id)
    
  })
  
  #### remove LD infos ####
  observeEvent(input$removeLD, {
    
    local({
      output$ld_scale <- renderPlot({
        NULL
      })
    })
    
    updateButton(session = session, inputId = "update_ld", 
                 disabled = TRUE)
    
    max_ld <- as.numeric(input$ld_range[2])
    min_ld <- as.numeric(input$ld_range[1])
    
    # update edges
    ld_edges <- values$all_edges #set all color on
    ld_edges$color <- "rgba(0, 0, 0, 0)" # set every edges at 
    ld_edges$xvalue <- 0
    ld_edges$title <- NA
    
    values$current_edges <- ld_edges
    # update network
    visNetworkProxy("my_network") %>%
      visUpdateEdges(edges = values$current_edges)
    
  })
  
  #### update LD infos ####
  observeEvent(input$update_ld, {

    max_ld <- as.numeric(input$ld_range[2])
    min_ld <- as.numeric(input$ld_range[1])
    
    # update edges
    ld_edges <- values$all_edges #set all color on
    
    if(sum((as.numeric(ld_edges$xvalue) > max_ld) | (as.numeric(ld_edges$xvalue) < min_ld)) > 0){
      ld_edges[(as.numeric(ld_edges$xvalue) > max_ld) | (as.numeric(ld_edges$xvalue) < min_ld),]$color <- "rgba(0, 0, 0, 0)"
    }
    values$current_edges <- ld_edges
    # update network
    visNetworkProxy("my_network") %>%
      visUpdateEdges(edges = values$current_edges)
    
  })
  
  
  #### BUILDING NETWORK ####
  observeEvent(input$buildNetwork, {
    js$collapse("input_box")
    js$collapse("selection_box")
    updateSelectInput(session = session, inputId = "snv_nodes_type", selected = 'pie_scores')
    updateButton(session = session, inputId = "buildNetwork", disabled = TRUE)
  })
  
  observeEvent(input$raw_data_rows_selected,{
    
    #setdiff(input$tbl_rows_selected, input$tbl_row_last_clicked)
    # updateButton(session = session, inputId = "buildNetwork", 
    #                       disabled = FALSE, style = "warning", label = "Update Network")
    
    if(length(input$raw_data_rows_selected) > MAX_VAR){
      shinyalert::shinyalert(title = "Selection limit", html = TRUE, type = "warning",
                             text = as.character(tags$div(style = "text-align:-webkit-center",
                                                          paste0("Variants selection from ",
                                                                 (length(input$raw_data_rows_selected) - MAX_VAR),
                                                                 " the Network visualisation is limit (",MAX_VAR, ")")
                             )))
      updateButton(session = session, inputId = "buildNetwork", 
                   disabled = TRUE, style = "danger")
    } else {
      updateButton(session = session, inputId = "buildNetwork", 
                   disabled = FALSE, style = "success")
    }
  })
  
  buildNetwork <- eventReactive(input$buildNetwork,{
    print(paste0("buildNetwork:", input$buildNetwork))
    
    if(is.null(values$my_data))
      return(NULL)
    
    my_data <- values$my_data
    
    if(input$network_type == "regul"){
      absolute_idx <- which(!my_data$non_synonymous)
      my_res <- values$all_regul_res
    } else {
      absolute_idx <- which(my_data$non_synonymous)
      my_res <- values$all_nonsyn_res
    }
    
    s1 = input$raw_data_rows_selected
    
    if(!(is.null(s1))){
      my_res <- my_res[sort(s1),]
    }
    
    # start by destroying everything
    if(!is.null(values$current_nodes)){
      print("cleaning nodes")
      visNetworkProxy("my_network") %>%
        visRemoveNodes(id = values$current_nodes$id)
    }
    
    if(!is.null(values$current_edges)){
      print("cleaning edges")
      visNetworkProxy("my_network") %>%
        visRemoveEdges(id = values$current_edges$id)
    }
    
    id <<- showNotification(paste("Building your wonderful network (~30sec)..."), duration = 0, type = "message")
    
    values$my_res <- my_res
    save(my_res, file = paste0(tmpDir, "/my_res.rda"))
    snv_edges <- build_snv_edges(values, "1", NULL, network_type = input$network_type) #create edges without real ld info
    save(snv_edges, file = paste0(tmpDir, "/snv_edges.rda"))
    snv_nodes <- build_snv_nodes(session_values = values, network_type = input$network_type, net = input$buildNetwork)
    save(snv_nodes, file = paste0(tmpDir, "/snv_nodes.rda"))
    
    all_scores_data <- build_score_pies(session_values = values,
                                        selected_scores = input$selected_scores,
                                        net = input$buildNetwork, inc = NULL)
    
    values$all_scores_data <- all_scores_data
    save(all_scores_data, file = paste0(tmpDir, "/all_scores_data.rda"))
    
    scores_data <- switch (input$snv_nodes_type,
                           # - pie_scores : relative score ordered by predictors
                           pie_scores = all_scores_data$relative$original,
                           # - pie_scores_group : relative score ordered by color
                           pie_scores_group = all_scores_data$relative$group_by,
                           # - pie_scores_abs : absolute score ordered by predictors
                           pie_scores_abs = all_scores_data$absolute$original,
                           # - pie_scores_abs_group : absolute score ordered by color
                           pie_scores_abs_group = all_scores_data$absolute$group_by,
                           # - pie_rank_na_last : global ranking NA last
                           pie_rank_na_last = all_scores_data$global$na_last,
                           # - pie_rank_na_mean : global ranking NA mean
                           pie_rank_na_mean = all_scores_data$global$na_mean
    )
    
    values$scores_data <- scores_data
    values$current_edges <- snv_edges
    values$current_nodes <- snv_nodes
    
    meta_values_mapped <- "gray"
    
    print(paste0("network_type : ",input$network_type, " ; current_nodes : ", nrow(values$current_nodes), " ; current_edges : ", nrow(values$current_edges)))
    
    updateButton(session = session, inputId = "buildNetwork", disabled = TRUE)
    
    if(!is.null(values$current_edges) && nrow(values$current_edges) > 0){
      updateButton(session = session, inputId = "runLD", 
                   disabled = FALSE)
      updateButton(session = session, inputId = "removeLD", 
                   disabled = FALSE)
    } else {
      updateButton(session = session, inputId = "runLD",
                   disabled = TRUE)
      updateButton(session = session, inputId = "removeLD", 
                   disabled = TRUE)
    } 
    
    if(!is.null(values$current_nodes) && nrow(values$current_nodes) > 0){
      updateButton(session = session, inputId = "update_metascore", 
                   disabled = FALSE)
    } else {
      updateButton(session = session, inputId = "update_metascore",
                   disabled = TRUE)
    }
    
    if (!is.null(id))
      removeNotification(id)
    
    return(list(nodes = values$current_nodes, edges = values$current_edges))
  })
  
  output$my_network <- renderVisNetwork({
    t0 <- Sys.time()
    vn_components <- buildNetwork()
    
    if(is.null(vn_components))
      return(NULL)
    
    local({
      output$scale <- renderPlot({
        draw_rank_palette(nbr_variants = nrow(vn_components$nodes), 
                          nodes_type = input$snv_nodes_type)
      })
    })
    
    save(vn_components, file = paste0(tmpDir, "/vn_components.rda"))
    print(Sys.time() - t0)
    visNetwork(vn_components$nodes, 
               vn_components$edges) %>%
      visEvents(doubleClick = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}") %>%
      visInteraction(tooltipDelay = 0, hideEdgesOnDrag = TRUE, navigationButtons = FALSE) %>%
      visOptions(highlightNearest = TRUE, clickToUse = FALSE, manipulation = FALSE, nodesIdSelection = FALSE) %>%
      visNodes(shapeProperties = list(useBorderWithImage = TRUE)) %>%
      visLayout(randomSeed = 2018) %>% 
      visPhysics("repulsion", repulsion = list(nodeDistance = 1000))
  })
  
  #### UPDATE PLOT ####
  buildPlot <- eventReactive(c(input$fetch_annotations, input$raw_data_rows_selected, input$network_type), {
    if(is.null(values$my_data))
      return(NULL)
    
    my_data <- values$my_data
    
    if(input$network_type == "regul"){
      absolute_idx <- which(!my_data$non_synonymous)
    } else {
      absolute_idx <- which(my_data$non_synonymous)
    }
    
    s1 = input$raw_data_rows_selected
    s2 = input$raw_data_rows_current
    
    if(!(is.null(s1))){
      my_data$selected <- FALSE
      my_data$color <- "blue"
      my_data[sort(absolute_idx[s1]),]$selected <- TRUE
      my_data$color[my_data$selected] <- "red"
    }
    
    return(my_data)
  })
  
  buildPlot_d <- buildPlot %>% debounce(5000)
  
  output$my_plot <- plotly::renderPlotly({
    pdf(NULL) # to avoid the production of Rplots.pdf
    
    #my_data <- buildPlot_d()
    my_data <- buildPlot()
    if(is.null(my_data) || nrow(my_data) == 0)
      return(NULL)
    
    cdts_region_line <- values$cdts_region
    cdts_region_line$color = "blue"
    cdts_region_line$size = 1
    cdts_region_line$shape = "line-ew"
    
    text_snp <- ifelse(test = my_data$no_cdts, yes = paste0(my_data$query," (no CDTS data)"), no = my_data$query)
    text_snp <- paste(text_snp, paste0("\nConsequences: ", my_data$consequences))
    xa <- list(title = levels(my_data$seqnames),
               zeroline = FALSE,
               showline = TRUE,
               showticklabels = TRUE,
               showgrid = TRUE)
    pdf(NULL) # to avoid the production of Rplots.pdf
    plotly::plot_ly(data = my_data,
                    type = "scatter",
                    mode = "markers",
                    x = ~start,
                    y = ~cdts_score,
                    symbol = ~I(shape), color = ~I(color), size = ~I(size),
                    marker = list(
                      opacity = 0.5
                    ),
                    showlegend = F, source = "subset") %>%
      plotly::add_trace(data = cdts_region_line,
                        type = 'scatter',
                        mode = 'lines',
                        x = ~start,
                        y = ~CDTS,
                        hoverinfo = "none",
                        line = list(color = 'gray'), 
                        opacity = 0.3,
                        showlegend = F) %>%
      plotly::add_trace(data = my_data,
                        x = ~start,
                        y = ~cdts_score,
                        text =  text_snp,
                        type = "scatter", 
                        mode = "markers",
                        #split = ~non_synonymous,
                        hoverinfo = 'text',
                        # legendgroup = 'In network',
                        # name = 'In network',
                        showlegend = T) %>%
      plotly::layout(xaxis = xa, dragmode =  "select")
  })
  
  #### UPDATE NODES CONTENT ####
  observe({
    
    if(!is.null(values$current_nodes) && nrow(values$current_nodes) > 0){
      svn_nodes <- values$current_nodes[values$current_nodes$group == "Variants",]
      
      if(input$update_metascore > 0)
        svn_nodes$image <- paste0("scores_figures/",input$snv_nodes_type, "_", svn_nodes$id,input$buildNetwork,input$update_metascore,".png")
      else 
        svn_nodes$image <- paste0("scores_figures/",input$snv_nodes_type, "_", svn_nodes$id,input$buildNetwork,".png")
      
      values$current_nodes$image <- as.character(values$current_nodes$image)
      values$current_nodes[values$current_nodes$group == "Variants",] <- svn_nodes #update current_nodes
      
      visNetworkProxy("my_network") %>%
        visUpdateNodes(nodes = svn_nodes)
      
    }
    
    if(!is.null(values$all_scores_data) && length(values$all_scores_data) > 0){
      
      all_scores_data <- values$all_scores_data
      
      scores_data <- switch (input$snv_nodes_type,
                             # - pie_scores : relative score ordered by predictors
                             pie_scores = all_scores_data$relative$original,
                             # - pie_scores_group : relative score ordered by color
                             pie_scores_group = all_scores_data$relative$group_by,
                             # - pie_scores_abs : absolute score ordered by predictors
                             pie_scores_abs = all_scores_data$absolute$original,
                             # - pie_scores_abs_group : absolute score ordered by color
                             pie_scores_abs_group = all_scores_data$absolute$group_by,
                             # - pie_rank_na_last : global ranking NA last
                             pie_rank_na_last = all_scores_data$global$na_last,
                             # - pie_rank_na_mean : global ranking NA mean
                             pie_rank_na_mean = all_scores_data$global$na_mean
      )
      
      values$scores_data <- scores_data
    }
  })
  
  #### UPDATE SCORES UI ####
  observe({
    if(is.null(values$my_data) || is.null(input$raw_data_rows_selected))
      return(NULL)
    
    # suppress cdts infos
    unwanted_cols <- c("cdts_score", "cdts_percentile")
    
    # suppress colnames with only NA values regarding the selected nodes
    my_data <- values$my_data
    
    if(input$network_type == "regul"){
      absolute_idx <- which(!my_data$non_synonymous)
      my_res <- values$all_regul_res
    } else {
      absolute_idx <- which(my_data$non_synonymous)
      my_res <- values$all_nonsyn_res
    }
    
    my_res <- my_res[sort(input$raw_data_rows_selected),]
    x <- values$res[values$res$query %in% my_res$query, colnames(values$res) %in% values$adjusted_scores]
    null_predictors <- names(which(apply(x, 2, function(i) sum(is.na(i))) == nrow(x)))
    
    predictors <- sort(values$adjusted_scores[!values$adjusted_scores %in% c(unwanted_cols, null_predictors)])
    
    if(length(values$adjusted_scores) > 0){
      updateSelectizeInput(session, "selected_scores",
                           choice = predictors,
                           selected = predictors)
    }
  })
  
  
  #### METASCORES ####
  observeEvent(input$update_metascore, {
    id <<- showNotification(paste("Updating your wonderful network..."), duration = 0, type = "message")
    t0 <- Sys.time()
    all_scores_data <- build_score_pies(session_values = values,
                                        selected_scores = input$selected_scores,
                                        net = input$buildNetwork, 
                                        inc = input$update_metascore)
    values$all_scores_data <- all_scores_data
    
    save(all_scores_data, file = paste0(tmpDir, "/all_scores_data.rda"))
    
    scores_data <- switch (input$snv_nodes_type,
                           # - pie_scores : relative score ordered by predictors
                           pie_scores = all_scores_data$relative$original,
                           # - pie_scores_group : relative score ordered by color
                           pie_scores_group = all_scores_data$relative$group_by,
                           # - pie_scores_abs : absolute score ordered by predictors
                           pie_scores_abs = all_scores_data$absolute$original,
                           # - pie_scores_abs_group : absolute score ordered by color
                           pie_scores_abs_group = all_scores_data$absolute$group_by,
                           # - pie_rank_na_last : global ranking NA last
                           pie_rank_na_last = all_scores_data$global$na_last,
                           # - pie_rank_na_mean : global ranking NA mean
                           pie_rank_na_mean = all_scores_data$global$na_mean
    )
    
    #save(scores_data, file = paste0(tmpDir, "/scores_data.rda"))
    values$scores_data <- scores_data
    svn_nodes <- values$current_nodes[values$current_nodes$group == "Variants",]
    svn_nodes$image <- paste0("scores_figures/",input$snv_nodes_type, "_", 
                              svn_nodes$id,input$buildNetwork, input$update_metascore,".png")
    
    visNetworkProxy("my_network") %>% 
      visUpdateNodes(nodes = svn_nodes)
    
    if (!is.null(id))
      removeNotification(id)
    
    print(Sys.time() - t0)
  })
  
  
  observeEvent(input$update_metascore, {
    
    # supprimer les nodes qui n'appartiennent pas aux scores selectionnés
    selected_scores <- c(input$adj_scores, input$raw_scores)
    score_nodes <- values$nodes[values$nodes$shape == "triangle",]
    score_edges <- values$edges[values$nodes$type == "score_edges",]
    nodes_id_2_remove <- score_nodes[!grepl(x = score_nodes$id, pattern = paste(selected_scores, collapse = "|")),]$id
    
    # nodes a garder
    nodes_2_keep <- values$nodes[!values$nodes$id %in% nodes_id_2_remove,]    
    
    visNetworkProxy("my_network") %>%
      visUpdateNodes(nodes = nodes_2_keep) %>%
      visUpdateEdges(edges = score_edges)
    
    visNetworkProxy("my_network") %>%
      visRemoveNodes(id = nodes_id_2_remove)
  })
  
  #### SCORES POPUP ####
  observeEvent(input$current_node_id, {
    values$selected_node <- input$current_node_id
    
    snv_score_details <- values$scores_data[values$scores_data$group == paste0("Scores_",values$selected_node), 
                                            c("id","color","label")]
    if(nrow(snv_score_details) > 0){
      color_code <- as.character(snv_score_details$color)
      text_color_code <- sapply(X = color_code, FUN = contrasting_text_color)
      snv_score_details$id <- gsub(x = snv_score_details$id, 
                                   pattern = paste0("_", values$selected_node), 
                                   replacement = "")
      snv_score_details <- snv_score_details[,-2] #suppression de la colonne 'color' 
      rownames(snv_score_details) <- NULL
      
      if(grepl(x = input$snv_nodes_type, pattern = "pie_rank")){
        colnames(snv_score_details) <- c("predictors", "rank")
        snv_score_details$rank <- round(x = as.numeric(as.character(snv_score_details$rank)),
                                        digits = 4)
      } else {
        colnames(snv_score_details) <- c("predictors", "value")
        snv_score_details$value <- round(x = as.numeric(as.character(snv_score_details$value)),
                                         digits = 4)
      }
      
      tab <- tableHTML::tableHTML(snv_score_details, theme = 'scientific', 
                                  rownames = FALSE)
      for (i in 1:nrow(snv_score_details)) { 
        tab <- tab %>% tableHTML::add_css_row(css = list(c('background-color','font-weight','color'), 
                                                         c(color_code[i],'bold', text_color_code[i])), rows = i+1) # first column is header
      }
      tab <- tab %>% tableHTML::add_css_column(css = list(c('background-color','text-align','font-weight','color'),
                                                          c('white','left','normal','gray')), columns = 1)
    } else {
      tab <- "No data for the selected predictors"
    }
    shinyalert::shinyalert(title = input$current_node_id, html = TRUE, 
                           closeOnEsc = TRUE,closeOnClickOutside = TRUE,
                           text = as.character(tags$div(style = "text-align:-webkit-center",
                                                        tab)))
  })
  
  
  observeEvent(input$select_button, {
    selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    print(paste('click on',selectedRow))
  })
  
  #### SCORES DESCRIPTION TABLE ####
  output$scores_description_table <- DT::renderDataTable({
    
    X <- readr::read_tsv(file = paste0(appDir, '/scores_description.tsv'))
    X <- split(X, X$group_name)
    
    new_df <- NULL
    for (n in names(X)){
      x <- X[[n]]
      new_row <- data.frame(Score = n, 
                            Id = paste(x$id, collapse = ";</br>"),
                            Range = paste(paste0(paste0("[",x$tolerate,"_",x$deleterious),paste0("]")), collapse = ";</br>"),
                            Ref = paste(unique(x$reference), 
                                        collapse = " "),
                            Description = paste(unique(x$description),
                                                collapse = " "),
                            Annotations = paste(unique(x$annotations),
                                                collapse = " ")
      )
      
      if(!is.null(new_df)){
        new_df <- rbind(new_df, new_row)
      } else {
        new_df <- new_row
      }
    }
    
    DT::datatable(cbind(' ' = '&oplus;', new_df), 
                  escape = FALSE,
                  rownames = FALSE,
                  extensions = c('Scroller','Buttons'),
                  selection = 'none',
                  options = list(
                    columnDefs = list(
                      list(visible = FALSE, targets = c(5,6)),
                      list(orderable = FALSE, className = 'details-control', targets = 0)
                    ),
                    dom = 'Bfrtip',
                    buttons = list(
                      list(
                        extend = 'print',
                        text = 'Print',
                        title = 'Predictors_description'
                      ),
                      list(
                        extend = 'csv',
                        filename = 'Predictors_description',
                        text = 'Download CSV'
                      ),
                      list(
                        extend = 'excel',
                        filename = 'Predictors_description',
                        text = 'Download XLSX'
                      ),
                      list(
                        extend = 'pdf',
                        title = 'Predictors_description',
                        text = 'Download PDF',
                        orientation = 'landscape'
                      )
                    ),
                    scrollX = TRUE,
                    scrollY = 500,
                    ordering = TRUE,
                    deferRender = TRUE,
                    scroller = TRUE),
                  callback = JS("
                                table.column(1).nodes().to$().css({cursor: 'pointer'});
                                var format = function(d) {
                                  return '<div style=\"background-color:#eee; padding: .5em;\"> <b>Description</b><br/>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;' +
                                          d[5] + '</br><b>Annotations</b><br/>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;' + d[6] + '</div>';
                                };
                                table.on('click', 'td.details-control', function() {
                                  var td = $(this), row = table.row(td.closest('tr'));
                                  if (row.child.isShown()) {
                                    row.child.hide();
                                    td.html('&oplus;');
                                  } else {
                                    row.child(format(row.data())).show();
                                    td.html('&CircleMinus;');
                                  }
                                });"
                  )
    )
  })
  
  
  #### ON STOP ####
  onStop(function() { 
    old_figures <- dir(path = "www/scores_figures/", full.names = T)
    file.remove(old_figures)
    temp_files <- dir(path = tmpDir, full.names = T, recursive = T)
    file.remove(temp_files)
    graphics.off()
    cat("Session stopped\n")
  })
}


#### OBSOLETE ####
# observeEvent(input$runLD, {
#   updateTabsetPanel(session, "results_tabset",
#                     selected = "ld_results"
#   )
# })
# output$ld_plot <- renderImage({  
#   filename <- "processing.gif"
#   if (!is.null(id))
#     removeNotification(id)
#   id <<- NULL
#   if(!is.null(values$ld)){
#     filename <- paste0(tmpDir,"/ld_figures/LD_plot_",input$ld_regions,".png")
#   }
#   
#   list(src = filename,
#        alt = "This is alternate text")
#   
# },deleteFile = FALSE)
# 
# #### BUILDING NETWORK ####
# observeEvent(input$buildNetwork, {
#   updateTabsetPanel(session, "results_tabset",
#                     selected = "network_results"
#   )
# })
#### LD ####
# observeEvent(input$update_ld, {
#   
#   # extraire les bon edges
#   max_ld <- input$ld_range[2]
#   min_ld <- input$ld_range[1]
#   
#   # supprimer les edges de type dist
#   dist_edges_id_2_remove <- values$current_edges[values$current_edges$type == "snv_dist_edges",]$id
#   
#   # supprimer les edges de type dist
#   empty_ld_edges_id_2_remove <- values$current_edges[(values$current_edges$type == "snv_ld_edges") & 
#                                                        (values$current_edges$title == "NA"),]$id
#   
#   
#   
#   
#   # supprimer les edges de type ld qui ne sont pas dans le range
#   ld_edges_id_2_remove <- values$current_edges[(values$current_edges$type == "snv_ld_edges") & 
#                                                  ((as.numeric(values$current_edges$xvalue) < min_ld | as.numeric(values$current_edges$xvalue) > max_ld)),]$id
#   
#   # ajouter les edges de type ld qui sont dans le range mais qui n'étaient pas dans le current_edges
#   ld_edges <- values$all_edges[values$all_edges$type == "snv_ld_edges",]
#   ld_edges_id_2_add <- ld_edges[(as.numeric(ld_edges$xvalue) >= min_ld) & (as.numeric(ld_edges$xvalue) <= max_ld),]
#   
#   edges_id_2_remove <- c(dist_edges_id_2_remove, ld_edges_id_2_remove)
#   
#   edges_2_keep <- values$current_edges[! values$current_edges$id %in% edges_id_2_remove,]
#   
#   # add newly computed LD
#   snv_ld_edges <- build_snv_edges(values, "1", NULL, network_type = input$network_type)
#   snv_ld_edges <- snv_ld_edges[((as.numeric(snv_ld_edges$xvalue) < min_ld | as.numeric(snv_ld_edges$xvalue) > max_ld)),]
#   edges_2_keep <- rbind(edges_2_keep, ld_edges_id_2_add, snv_ld_edges)
#   
#   values$current_edges <- edges_2_keep
#   
#   visNetworkProxy("my_network") %>%
#     visUpdateEdges(edges = edges_2_keep)
#   
#   visNetworkProxy("my_network") %>%
#     visRemoveEdges(id = edges_id_2_remove)
# })

#### DISTANCE ####
# observeEvent(input$update_dist, {
#   # supprimer les edges de type ld
#   ld_edges_id_2_remove <- values$current_edges[values$current_edges$type == "snv_ld_edges",]$id
#   
#   # supprimer les edges entre les nodes trop eloignés
#   max_dist <- 1000 * input$dist_range
#   dist_edges_id_2_remove <- values$current_edges[(values$current_edges$type == "snv_dist_edges") & (as.numeric(values$current_edges$xvalue) > max_dist),]$id
#   
#   #dist_edges_id_2_remove <- dist_edges_id_2_remove[as.numeric(dist_edges_id_2_remove$xvalue) > max_dist,]$id
#   
#   # ajouter les edges de type dist qui sont dans le range mais qui n'étaient pas dans le current_edges
#   dist_edges <- values$all_edges[values$all_edges$type == "snv_dist_edges",]
#   dist_edges_2_add <- dist_edges[(as.numeric(dist_edges$xvalue) <= max_dist),]
#   
#   edges_id_2_remove <- c(ld_edges_id_2_remove, dist_edges_id_2_remove)
#   
#   # edges a garder
#   edges_2_keep <- values$current_edges[! values$current_edges$id %in% edges_id_2_remove,]
#   edges_2_keep <- rbind(edges_2_keep, dist_edges_2_add)
#   
#   values$current_edges <- edges_2_keep
#   
#   visNetworkProxy("my_network") %>%
#     visUpdateEdges(edges = edges_2_keep)
#   
#   visNetworkProxy("my_network") %>%
#     visRemoveEdges(id = edges_id_2_remove)
# })

#### FOCUS ####
# observe({
#   updateSelectInput(session, "focus",
#                     choices = c("None", as.character(values$current_nodes$id)),
#                     selected = "None"
#   )
# })
# 
# observe({
#   if(input$focus != "None"){
#     visNetworkProxy("my_network") %>%
#       visFocus(id = input$focus, scale = 1)
#   } else {
#     visNetworkProxy("my_network") %>%
#       visFit()
#   }
# })

#### HIGHLIGHT ####
# observe({
#   if("cadd.consdetail" %in% colnames(values$annotations)){
#     choices <- unique(unlist(strsplit(x = as.character(values$annotations$cadd.consdetail), split = ","))) 
#     updateSelectInput(session, "highlight",
#                       choices = c("None", choices),
#                       selected = "None")
#   }
# })

# observe({
#   svn_nodes <- values$current_nodes[values$current_nodes$group == "Variants",]
#   
#   # remise à zero
#   svn_nodes$shape <- "circularImage"
#   
#   if(input$highlight != "None"){
#     # nodes to highlight
#     svn_id_2_highlight <- as.character(values$annotations$query)[grepl(x = as.character(values$annotations$cadd.consdetail),
#                                                                        pattern = input$highlight)]
#     svn_nodes[svn_nodes$id %in% svn_id_2_highlight,]$shape <- "image"
#     values$current_nodes$shape <- as.character(values$current_nodes$shape)
#     values$current_nodes[values$current_nodes$group == "Variants",] <- svn_nodes #update current_nodes
#   }
#   
#   visNetworkProxy("my_network") %>%
#     visUpdateNodes(nodes = svn_nodes)
# })
#### SCORES STATS ####

# compute_scores_missing_data <- eventReactive(input$get_predictors_info, {
#   if(length(input$predictors) > 0){
#     # non_null_raw_scores <- values$raw_scores[!sapply(values$raw_scores, is.null)]
#     # non_null_adj_scores <- values$adjusted_scores[!sapply(values$adjusted_scores, is.null)]
#     non_null_raw_scores <- values$raw_scores
#     non_null_adj_scores <- values$adjusted_scores
#     all_scores <- c(non_null_adj_scores, non_null_raw_scores)
#     all_scores <- all_scores[input$predictors]
#     all_scores <- data.frame(all_scores)
#     
#     scores_stats <- do.call("rbind",
#                             lapply(X = all_scores,
#                                    FUN = function(x) {
#                                      data.frame(Missing_values = sum(is.na(x)),
#                                                 Min = signif(x = min(x, na.rm = T), digits = 3),
#                                                 #First_Qt = quantile(x, probs = 0.25, names = F, na.rm = T),
#                                                 Mean = signif(x = mean(x, na.rm = TRUE), digits = 3),
#                                                 #Median = median(x, na.rm = TRUE),
#                                                 #Third_Qt =  quantile(x, probs = 0.75, names = F, na.rm = T),
#                                                 Max = signif(x = max(x, na.rm = TRUE), digits = 3))
#                                    }
#                             )
#     )
#     
#     pretty_rownames <- row.names(scores_stats)
#     pretty_rownames <- gsub(x = pretty_rownames, pattern = "score|rawscore|rankscore", replacement = "")
#     pretty_rownames <- gsub(x = pretty_rownames, replacement = "", pattern = "\\.$")
#     pretty_rownames <- gsub(x = pretty_rownames, replacement = " ", pattern = "\\.")
#     pretty_rownames <- gsub(x = pretty_rownames, pattern = "_(.*)_", replacement = " (\\1)")
#     
#     row.names(scores_stats) <- pretty_rownames
#     
#     return(scores_stats)
#   }
# })
# 
# compute_scores_stats <- eventReactive(input$get_predictors_info,{
#   
#   # non_null_raw_scores <- values$raw_scores[!sapply(values$raw_scores, is.null)]
#   # non_null_adj_scores <- values$adjusted_scores[!sapply(values$adjusted_scores, is.null)]
#   non_null_raw_scores <- values$raw_scores
#   non_null_adj_scores <- values$adjusted_scores
#   all_scores <- c(non_null_adj_scores, non_null_raw_scores)
#   all_scores <- all_scores[input$predictors]
#   all_scores <- data.frame(all_scores)
#   
#   df.m <- reshape2::melt(all_scores)
#   p <- ggplot(data = df.m, aes(x=variable, y=value)) + geom_boxplot()
#   p <- p + facet_wrap( ~ variable, scales="free", strip.position = "left")
#   p <- p + theme_bw() + 
#     theme(axis.title.y = element_blank(), 
#           axis.title.x = element_blank(), 
#           axis.text.x = element_blank())
#   
#   return(p)
# })
# 
# output$scores_stats <- renderPlot({
#   #suppressWarnings(compute_scores_stats())
# })
# 
# output$scores_missing_data <- DT::renderDataTable({
#   # details <- as.matrix(compute_scores_missing_data())
#   # # if(is.null(detail)){
#   # #   detail <- as.matrix(data.frame(waiting = ""))
#   # # }
#   # DT::datatable(details,
#   #               rownames = TRUE,
#   #               options = list(scrollX = TRUE,
#   #                              lengthChange = FALSE,
#   #                              searching = FALSE))
# })
# 
# 
# #### SCORES CORRELATION MATRICE ####
# compute_scores_corr <- eventReactive(input$get_predictors_info, {
#   sub_matrice <- diag(nrow = 10, ncol = 10)
#   
#   if(length(input$predictors) > 0){
#     sub_matrice <- scores_correlation_matrice[rownames(scores_correlation_matrice) %in% input$predictors,
#                                               colnames(scores_correlation_matrice) %in% input$predictors]
#   }
#   
#   return(sub_matrice)
# })
# 
# output$scores_corr <- renderD3heatmap({
#   # m <- compute_scores_corr()
#   # if(!is.null(m)){
#   #   d3heatmap(x = m, Colv = "Rowv", 
#   #             Rowv = NULL, colors = cor_color_breaks, 
#   #             cexCol = .7, na.rm = F, cexRow = .7)
#   # }
#   
# })
#cadd coonsequences
# output$consequences <- renderC3PieChart({
#   if(!is.null(values$res) && nrow(values$res) > 0){
#     d1 <- data.frame(non_synonymous = as.numeric(sum(sapply(values$res$cadd.consequence, function(y) "NON_SYNONYMOUS" %in% y))), 
#                      other = as.numeric((nrow(values$res) - sum(sapply(values$res$cadd.consequence, function(j) "NON_SYNONYMOUS" %in% j)))))
#     d1 %>% c3() %>% c3_pie()
#   }
# })
# 
# #cadd annotations
# output$annotations <- renderC3PieChart({
#   if(!is.null(values$res) && nrow(values$res) > 0){
#     annotype_counts <- as.numeric(table(unlist(values$res$cadd.annotype)))
#     annotyp_names <- names(table(unlist(values$res$cadd.annotype)))
#     d1 <- data.frame(matrix(annotype_counts, nrow = 1))
#     colnames(d1) <-  annotyp_names
#     d1 %>% c3() %>% c3_pie()
#   }
# })
# 
# #genes
# output$genenames <- renderC3PieChart({
#   if(!is.null(values$res) && nrow(values$res) > 0){
#     if(is.null(values$res$cadd.gene)){
#       genenames <- values$res$cadd.gene.genename
#       genenames[is.na(genenames)] <- "none"
#       genenames[is.null(genenames)] <- "none"
#     } else {
#       genenames <- sapply(X = values$res$cadd.gene, FUN = function(x) x$genename)
#       genenames[sapply(genenames, function(x) is.null(x))] <- "none"
#     }
#     
#     genenames_counts <- as.numeric(table(unlist(genenames)))
#     genenames_names <- names(table(unlist(genenames)))
#     d1 <- data.frame(matrix(genenames_counts, nrow = 1))
#     colnames(d1) <-  genenames_names
#     d1 %>% c3() %>% c3_pie()
#   }
# })
# 
# output$cdts_scores <- renderC3PieChart({
#   if(!is.null(values$res) && nrow(values$res) > 0){
#     d1 <- values$res
#     d1 %>% c3(y = 'cdts_score', x = 'query') %>% c3_bar() %>% legend(hide= T) %>% grid('y', show = F, lines = data.frame(value = 0)) %>% tickAxis("x", values = c(1))
#   }
# })
# 
# output$cdts_percentile <- renderC3PieChart({
#   if(!is.null(values$res) && nrow(values$res) > 0){
#     d1 <- values$res
#     d1 %>% c3(y = 'cdts_percentile', x = 'query') %>% c3_bar() %>% legend(hide= T) %>% grid('y', show = F, lines = data.frame(value = 0)) %>% tickAxis("x", values = c(1))
#   }
# })
# getSelectedNodeID <- eventReactive(input$current_node_id, {
#   cat("2",input$current_node_id,"\n")
#   if(is.null(values$my_res)){
#     color_code = "#FFFFFF"
#     snv_score_details <- data.frame("predictors" = "None", "value" = 0)
#     return(DT::datatable(snv_score_details,
#                          rownames = FALSE,
#                          options = list(scrollX = TRUE,
#                                         lengthChange = FALSE,
#                                         searching = FALSE)))
#     
#   }
#   
#   if(!is.null(values$selected_node) && values$selected_node == '') {
#     values$selected_node <- as.character(values$my_res$query)[1]
#   }
#   
#   snv_score_details <- values$scores_data$nodes[values$scores_data$nodes$group == paste0("Scores_",values$selected_node), c("id","color","label")]
#   #print(snv_score_details)
#   color_code <- as.character(snv_score_details$color)
#   snv_score_details$id <- gsub(x = snv_score_details$id, pattern = paste0("_", values$selected_node), replacement = "")
#   snv_score_details <- snv_score_details[,-2] #suppression de la colonne 'color' 
#   colnames(snv_score_details) <- c("predictors", "value")
#   return(snv_score_details)
# })
# 
# #### SVN SCORES DETAIL ####
# output$snv_score_details <- DT::renderDataTable({
#   DT::datatable(getSelectedNodeID(),
#                 rownames = FALSE,
#                 options = list(scrollX = TRUE,
#                                lengthChange = FALSE,
#                                searching = FALSE))
# }%>% DT::formatStyle(
#   'value',
#   backgroundColor = DT::styleEqual(as.character(snv_score_details$value), color_code)
# ))
# 
# output$snv_score_details_id <- renderText({
#   if(is.null(values$my_res))
#     return(NULL)
#   #print(values$selected_node)
#   if(!is.null(values$selected_node) && values$selected_node == '') 
#     values$selected_node <- as.character(values$my_res$query)[1]
#   
#   values$selected_node
# })
#### VARIANT SELECTION FROM PLOTLY ####
# output$selection <- renderUI({
#   event.data <- event_data("plotly_selected", source = "subset")
#   
#   if(is.null(event.data) == T ) return("Please select an area containing variants...")
#   if(!"curveNumber" %in% colnames(event.data)) return("Please select an area containing variants...")
#   
#   x <- values$my_data[subset(event.data, curveNumber == 0)$pointNumber + 1,]
#   x <- x[, c("seqnames","start","query","type","cdts_score","consequences","non_synonymous")] 
#   
#   if(nrow(x) > MAX_VAR){
#     plotted_variants <- list(
#       paste0(nrow(x), " variants selected. The first ", MAX_VAR, " will be plotted."),
#       br(),
#       DT::datatable(x[1:MAX_VAR,],
#                     rownames = TRUE,
#                     options = list(scrollX = TRUE,
#                                    #pageLength = MAX_VAR,
#                                    lengthChange = FALSE,
#                                    searching = FALSE)
#       )
#     )
#   } else {
#     plotted_variants <- list(
#       paste0(nrow(x), " variants selected and plotted."),
#       br(),
#       DT::datatable(x,
#                     rownames = TRUE,
#                     options = list(scrollX = TRUE,
#                                    #pageLength = MAX_VAR,
#                                    lengthChange = FALSE,
#                                    searching = FALSE))
#     )
#   }
#   
#   network_choices <- c()
#   if(sum(x$non_synonymous) > 0){
#     network_choices <- c(network_choices, 
#                          "non synonymous variants" = "non_syn")
#   }
#   
#   if(sum(!x$non_synonymous) > 0){
#     network_choices <- c(network_choices, 
#                          "synonymous and non-coding variants" = "regul")
#   }
#   
#   updateSelectInput(session, "network_type",
#                     choices = network_choices,
#                     selected = network_choices[1]
#   )
#   
#   return(plotted_variants)
# })



