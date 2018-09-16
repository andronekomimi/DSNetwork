require(shiny)
require(d3heatmap)
require(visNetwork)
require(myvariant)
require(shinyBS)
require(plotly)
require(shinyjs)
require(tableHTML)
require(magrittr)
require(shinyalert)

options(shiny.trace = FALSE)

server <- function(input, output, session) {
  
  tmpDir <- tempdir() 
  #dir.create(path = paste0(tmpDir, "ld_figures"), showWarnings = F)
  #dir.create(path = paste0(tmpDir, "objects"), showWarnings = F)
  
  source('helper.R', local = TRUE)
  
  if(is.local()){
    preload <- list(
      locus_80 = "rs2660776\nrs2660774\nrs2575792\nrs2660773\nrs2660772\nrs1370042\nrs1160076\nrs2018075\nrs62259569\nrs2575799\nrs2253915\nrs2253916\nrs2253921\nrs2575802\nrs138876277\nrs1370044\nrs2575805\nrs144055118\nrs35711034\nrs35733956\nrs74520234\nrs147080240\nrs4859022\nrs2262594\nrs2254904\nrs729942\nrs1370045\nrs35870472\nrs34339032\nrs2660785\nrs2575815\nrs2575817\nrs2660787\nrs2575818\nrs2660788\nrs2118082\nrs2244125\nrs2660790\nrs2660791\nrs2660792\nrs2660794\nrs61340205\nrs1821751\nrs2660796\nrs2196523\nrs1370040\nrs2660797\nrs2575783\nrs2660798\nrs2575780\n3:87044690:A:AAT\nrs2575775\nrs2660799\nrs2575774\nrs2575772\nrs2043663\nrs13066793\nrs1960268\nrs1025631\nrs4859107\nrs201414876",
      locus_78 = "rs73167030\nrs12159970\nrs73167031\nrs56150793\nrs6001911\nrs74767555\nrs111843590\nrs6001912\nrs10483203\nrs56283550\nrs12158872\nrs17001907\n22:40836155:CAAAA:C\nrs56215843\nrs61675795\nrs73167042\nrs6001915\nrs17001915\nrs73167045\nrs73167052\nrs5995856\nrs17001920\nrs6001920\nrs73167053\nrs5995860\nrs28419341\nrs28552449\nrs112375848\nrs55864500\nrs183438976\nrs73167058\nrs5995862\nrs17001943\nrs58778028\nrs73167063\nrs5995864\nrs75523053\nrs74486969\nrs61211547\nrs12159787\nrs10483204\nrs145014115\nrs4299422\nrs2899337\nrs201540223\nrs73167066\nrs71777635\nrs73167067\nrs6001930\nrs17001974\nrs6001931\nrs6001932\nrs73167069\nrs73167072\nrs73167073\nrs73167076\nrs17001977\nrs3827381\nrs3827382\nrs73167079\nrs10483205\nrs112880707\nrs73167080\nrs73167082\nrs113409089\nrs74278065\nrs6001935\nrs56135013\nrs6001937\nrs73167089\nrs113798157\nrs5995867\nrs73167090\nrs77426923\nrs56182212\nrs73167092\nrs199614224\nrs6001939\nrs73167093\nrs6001942\nrs138175438\nrs73167096\nrs73167097\nrs73167098\nrs373996442\nrs113966362\nrs73167101\nrs17001993\nrs17001994\nrs73169003\nrs10483206\nrs6001946\nrs6001949\nrs66987842\nrs17001997\nrs6001950\nrs56124626\nrs141580207\nrs112034576\nrs55993771\nrs142270009\nrs6001954\nrs5995870\nrs57693796\nrs5995871\nrs139903991\nrs55849114\nrs55722862\nrs55960299\nrs56334449\nrs55662398\nrs6001962\nrs150091203\nrs932379\nrs73169026\nrs73169028\nrs73169032\nrs55674424\nrs5995875\nrs6001965\nrs56092708\nrs16985899\nrs17002019\nrs17002020\nrs201432472\nrs73169036\nrs12160047\nrs56047425\nrs56411183\nrs73169040\nrs117423948\nrs5995876\nrs138476044\nrs146829921\nrs73169043\nrs116973859\nrs28444678\nrs375665020\nrs17002024\nrs57705902\nrs79799751\nrs6001973\nrs6001974\nrs73169051\nrs143586053\nrs17002026\nrs17002027\nrs28567076\nrs73169056\nrs55853923\nrs17002030\nrs73169057\nrs73169058\nrs73169060\nrs141597790\nrs142127948\nrs17002034\nrs17002036\nrs6001979\nrs17002038\nrs6001980\nrs55719110\nrs55997921\nrs55775031\nrs139645469\nrs73169065\nrs59047419\nrs55669535\nrs192603101\nrs73169071\nrs5995881\nrs73169072\nrs6001981\nrs145351698\nrs73169077\nrs73169078\nrs6001982\nrs73169083\nrs73169084\nrs112497956\nrs6001983\nrs718193\nrs56108505\nrs73169087\nrs73169089\nrs6001984\nrs71785733\nrs73169091\nrs117669949\nrs78694459\nrs73169096\nrs73169097\nrs145005064\nrs73169098\nrs145624734\nrs56051217",
      locus_70 = "rs4889891\nrs9905914\nrs9896202\nrs745571\nrs745570\nrs2587505\nrs139427260\nrs71161686\nrs2587507",
      locus_60 = "13:32866927:CAATAAATAAATA:CAATAAATA\nrs56084662\nrs11571815\nrs11571818\nrs11571833",
      locus_38 = "rs62485509\nrs720475",
      locus_6 = "rs181595184\nrs12129763\nrs6686987\nrs55899544\nrs6427943\nrs6678914\nrs6703244\nrs12028423\nrs12131882\nrs12132085\nrs12129456\nrs12129536\nrs10920365\nrs2167588\nrs4950774\nrs12032424\nrs4950775\nrs4950836\nrs6143572\nrs896548\nrs4245706\nrs12032080\nrs3795598\nrs12143329\nrs4950837\nrs12021815\nrs60573451\nrs140473199",
      locus_2 = "rs11102694\nrs2358994\nrs2358995\nrs7547478\nrs7513707\nrs12022378\nrs3761936\nrs11102701\nrs11102702\nrs12046289\nrs112974454",
      locus_0 = "rs6864776\n5:44527739:A:ATACT\nrs4634356\nrs1905192\nrs4866905\nrs1482663\n5:44496660:A:AG\nrs7710996\nrs6451763\n5:44527050:C:A\nrs1351633\nrs1384453\nrs1482665\nrs983940\nrs6897963\nrs1384454\nrs10079222\nrs7736427\nrs10512860\nrs4866776\nrs1482690\nrs12516346\nrs1482684\n5:44496659:T:TA\nrs1482691\nrs7724859\nrs2128430\nrs7707044\nrs1905191\nrs1120718\nrs4866899\nrs7712213\nrs6451762\nrs7703171\nrs6879342")
    
    #### CONF ####
    appDir <- "/Users/nekomimi/Workspace/dsnetwork/DSNetwork/"
    dataDir <- "/Users/nekomimi/Workspace/dsnetwork/DSNetwork/data/"
    path_to_victor <- paste0(appDir, "softs/VICTOR/")
    python_path <- "/Users/nekomimi/anaconda/bin/python"
    tabix_path <- '/usr/local/bin/tabix'
  } else {
    preload <- list(
      locus_0 = "rs6864776\n5:44527739:A:ATACT\nrs4634356\nrs1905192\nrs4866905\nrs1482663\n5:44496660:A:AG\nrs7710996\nrs6451763\n5:44527050:C:A\nrs1351633\nrs1384453\nrs1482665\nrs983940\nrs6897963\nrs1384454\nrs10079222\nrs7736427\nrs10512860\nrs4866776\nrs1482690\nrs12516346\nrs1482684\n5:44496659:T:TA\nrs1482691\nrs7724859\nrs2128430\nrs7707044\nrs1905191\nrs1120718\nrs4866899\nrs7712213\nrs6451762\nrs7703171\nrs6879342")
    
    #### CONF ####
    appDir <- "/srv/shiny-server/dsnetwork/"
    dataDir <- "/mnt/apps_data/dsnetwork/"
    path_to_victor <- "/mnt/apps_softs/dsnetwork/VICTOR/"
    python_path <- "/Users/nekomimi/anaconda/bin/python" ###
    tabix_path <- '/mnt/apps_softs/dsnetwork/TABIX/tabix' ###
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
  
  # A notification ID
  id <- NULL
  
  cat("Session start\n")
  
  observeEvent(c(input$query_file, input$query), ({
    
    if(!is.null(input$query_file) && input$query_file$size > 0){
      shinyBS::updateButton(session = session, inputId = "fetch_annotations",
                            disabled = FALSE)
    } else {
      if(nchar(input$query) > 0){
        shinyBS::updateButton(session = session, inputId = "fetch_annotations",
                              disabled = FALSE)
      } else {
        shinyBS::updateButton(session = session, inputId = "fetch_annotations",
                              disabled = TRUE)
      }
    }
  }))
  
  #### LOAD PRESET QUERY ####
  observeEvent(input$preload_loci, {
    updateTextAreaInput(session = session, inputId = "query", value = preload[[input$preload]])
  })
  
  #### LOAD DEMO QUERY ####
  observeEvent(input$load_demo, {
    updateTextAreaInput(session = session, inputId = "query", value = preload[["locus_0"]])
  })
  
  #### INPUT DATA MANAGEMENT ####
  transform_query <- function(string){
    if(nchar(input$query) > 0){
      string <- unlist(strsplit(x = string, split = "\n"))
      string <- gsub(x = string, pattern = " +$", replacement = "")
      
      modstring <- sapply(X = string,  FUN = function(x) {
        if(base::startsWith(x = x, prefix = "rs")){
          return(x)
        }  else { 
          x <- unlist(strsplit(x = x, split = ":"))
          if(length(x) != 4){
            return('FAIL')
          } else {
            formatSingleHgvs(x[1], as.numeric(x[2]), x[3], x[4])
          }
        }
      })
      return(modstring)
    }
  }
  
  transform_query_file <- function(filepath){
    if(file.exists(filepath)){
      string <- read.csv(file = filepath, header = F, stringsAsFactors = F)$V1
      string <- gsub(x = string, pattern = " +$", replacement = "")
      
      modstring <- sapply(X = string,  FUN = function(x) {
        if(base::startsWith(x = x, prefix = "rs")){
          return(x)
        }  else { 
          x <- unlist(strsplit(x = x, split = ":"))
          if(length(x) != 4){
            return('FAIL')
          } else {
            formatSingleHgvs(x[1], as.numeric(x[2]), x[3], x[4])
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
    
    shinyBS::updateButton(session = session, inputId = "fetch_annotations", 
                          disabled = TRUE)
    
    withProgress(message = 'Fetching annotations', value = 0, {
      n <- 10
      modstring <- c()
      if(!is.null(input$query_file))
        modstring <- c(modstring, transform_query_file(input$query_file$datapath))
      
      if(nchar(input$query) > 0)
        modstring <- c(modstring, transform_query(input$query))
      
      valid_transfo <- unique(modstring[!grepl(x = modstring, pattern = 'FAIL')])
      print(valid_transfo)
      save(modstring, file = paste0(tmpDir, 'modstring.rda'))
      
      fail_transfo <- names(modstring[grep(x = modstring, pattern = 'FAIL')])
      if(length(modstring) > 0){
        if(length(fail_transfo) > 0){
          res <- paste0('Id recognition fails for: ', paste(fail_transfo, collapse = ","), '.')
          shinyBS::createAlert(session = session, anchorId = "alert_res",
                               alertId = "alert1", title = "Id recognition",
                               content = res , append = TRUE, style = "warning")
        }
      }
      
      if(length(valid_transfo) > 0){
        #### run myvariant ####
        incProgress(1/n, detail = "Interrogating MyVariant.info...")
        res <- as.data.frame(getVariants(hgvsids = valid_transfo,
                                         verbose = F, return.as = "DataFrame",
                                         fields = c("dbsnp","cadd","dbnsfp")))
        
        
        
        #### add metascores columns ####
        res$linsight <- res$cdts_score <- res$cdts_percentile <- NA
        res$bayesdel <- res$iwscoring_known <- res$iwscoring_novel <- NA
        
        values$res <- res
        remove(res)
        if(!is.null(values$res) && nrow(values$res) > 0){
          
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
          
          #### fetch linsight ####
          incProgress(1/n, detail = "Extracting LINSIGHT scores...")
          requested_chromosomes <- seqlevelsInUse(global_ranges)
          
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
            
            #### fetch cdts ####
            incProgress(1/n, detail = "Extracting CDTS data...")
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
            incProgress(4/n, detail = "Interrogating SNPNexus platform...")
            path_to_snpnexus <- paste0(appDir, "scripts/")
            
            # create temp snps file
            filename <- tempfile(tmpdir = tmpDir, fileext = ".txt")
            print(filename)
            createSNPnexusInput(session_values = values, filename = filename)
            value <- runSNPnexus(python_path = python_path,
                                 path_to_snpnexus = path_to_snpnexus, 
                                 filename = filename, 
                                 waiting_time = input$waiting)
            
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
              converted_id <- sapply(X = snpnexus_res$id, FUN = function(x) snpnexusIDconversion(x))
              snpnexus_res$converted_id <- converted_id
              # merge results based on common id
              pre_res <- merge(y = snpnexus_res, x = pre_res, by.y = "converted_id", by.x = "X_id", all = T)
              # metascores
              pre_res$iwscoring_novel <- pre_res$iwscoren6
              pre_res$iwscoring_known <- pre_res$iwscorek10
              # non-coding scores
              colnames(pre_res)[colnames(pre_res) == "deepsea.x"] <- "deepseq_sig_log2@" #deppsea : Lower score indicates higher likelihood of functional significance of the variant. But log2 is in the "good" direction
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
              shinyBS::createAlert(session = session, anchorId = "alert_res",
                                   alertId = "alert4", title = "Data retrieval",
                                   content = "SNPnexus taking too long to answer, sorry..." , 
                                   append = TRUE, style = "warning")
            }
          } else {
            incProgress(4/n, detail = "Skipping SNPNexus...")
          }
          
          res <- values$res
          save(res, file = paste0(tmpDir, '/res.rda'))
          
          incProgress(2/n, detail = "Aggregating data...")
          #### create metascore pies ####
          compute_absolute_metascore(values)
          
          
          common_fields <- colnames(values$res)[grepl(x = colnames(values$res), pattern = "query|X_id|cadd.ref|cadd.alt")]
          
          common_scores <- c("linsight", 
                             "bayesdel",
                             "cdts_score",
                             "cdts_percentile")
          
          snpnexus_scores <- c("gwava_region","gwava_tss","gwava_unmatched",
                               "eigen_nc","eigen_pc_nc",
                               "deepseq_sig_log2",
                               "fathmm_nc","fitcons_nc","funseq","remm",
                               "iwscoring_known","iwscoring_novel")
          
          dbnsfp_rankscores <- colnames(values$res)[grepl(x = colnames(values$res), pattern = "dbnsfp.*.rankscore")]
          
          cadd_raw_scores <- read.csv(file = paste0(dataDir, 'CADD_scores_from_myvariant.info.tsv'), header = T, sep = "\t", stringsAsFactors = F)
          cadd_raw_scores <- cadd_raw_scores[cadd_raw_scores$is_included == "x",]$field
          
          #### split res in 2 -> non-syn vs other
          non_syn_res <- sapply(values$res$cadd.consequence, function(x) "NON_SYNONYMOUS" %in% x)
          values$all_nonsyn_res <- values$res[non_syn_res, (colnames(values$res) %in% c(common_fields, common_scores, dbnsfp_rankscores)) ]
          values$all_regul_res <- values$res[!non_syn_res, (colnames(values$res) %in% c(common_fields, common_scores, snpnexus_scores, cadd_raw_scores)) ]
          
          if(nrow(values$all_nonsyn_res) > 0){
            adjusted_scores <- sapply(X = c(dbnsfp_rankscores, common_scores), FUN = function(x) extract_score_and_convert(values$all_nonsyn_res, score_name = x))
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
              values$all_nonsyn_res[[s]] <- adjusted_scores[[s]]
            }
            
          } else {
            adjusted_scores <- NULL
          }     
          
          
          if(nrow(values$all_regul_res) > 0){
            raw_scores <- sapply(X = c(cadd_raw_scores, snpnexus_scores, common_scores), FUN = function(x) extract_score_and_convert(values$all_regul_res, score_name = x))
            switch (class(raw_scores),
                    matrix = {
                      rownames(raw_scores) <- NULL
                      raw_scores <- as.list(data.frame(raw_scores))
                    },
                    numeric = {
                      names(raw_scores) <- gsub(x = names(raw_scores), pattern = "(.*[a-z]).\\d+.?\\d+$", replacement = "\\1")
                      raw_scores <- lapply(split(raw_scores, names(raw_scores)), unname)
                    }
            )
            
            raw_scores <- raw_scores[sapply(X = raw_scores, FUN = function(x) sum(is.na(x)) != length(x))] 
            
            for(s in names(raw_scores)){
              #print(raw_scores[[s]])
              values$all_regul_res[[s]] <- raw_scores[[s]]
            }
            
          } else {
            raw_scores <- NULL
          }
          
          
          if(nrow(values$all_nonsyn_res) > MAX_VAR){
            #print("faut trier all_nonsyn_res")
            values$nonsyn_res <- values$all_nonsyn_res[1:MAX_VAR,]
          } else {
            #print("pas besoin de trier all_nonsyn_res")
            values$nonsyn_res <- values$all_nonsyn_res
          }
          
          if(nrow(values$all_regul_res) > MAX_VAR){
            #print("faut trier all_regul_res")
            values$regul_res <- values$all_regul_res[1:MAX_VAR,]
          } else {
            #print("pas besoin de trier all_regul_res")
            values$regul_res <- values$all_regul_res
          }
          
          
          save(raw_scores, adjusted_scores, file = paste0(tmpDir, "/scores.rda"))
          values$raw_scores <- names(raw_scores)
          values$adjusted_scores <- names(adjusted_scores)
          
          if(is.null(values$notfound_id)){
            shinyBS::createAlert(session = session, anchorId = "alert_res",
                                 alertId = "alert2", title = "Data retrieval",
                                 content = "Annotations found for all variants" , append = TRUE, style = "success")
          } else {
            shinyBS::createAlert(session = session, anchorId = "alert_res",
                                 alertId = "alert2", title = "Data retrieval",
                                 content = paste0("No Annotations for the following variants: ", 
                                                  values$notfound_id,".") , append = TRUE, style = "warning")
          }
        } else {
          shinyBS::createAlert(session = session, anchorId = "alert_res",
                               alertId = "alert3", title = "Data retrieval",
                               content = "ERROR : No Annotations found !", 
                               append = TRUE, style = "danger")
        }
      }
      
      local({
        #### DISPLAY RESULTS TABLE FOR VARIANTS SELECTION ####
        annotations_fields <- c("query","X_id", "dbsnp.chrom", "dbsnp.hg19.start", "cadd.consequence", common_scores)
        annotations_infos <- values$res[,annotations_fields]
        annotations_infos$cadd.consequence <- sapply(annotations_infos$cadd.consequence, FUN = function(x) paste(x, collapse = ","))
        # annotations_infos$more <- shinyInput(actionButton, nrow(annotations_infos), 
        #                                      'button_', label = "Fire", 
        #                                      onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
        values$annotations <- annotations_infos
        save(annotations_infos, file = paste0(tmpDir, "/annotations.rda"))
        
        output$raw_data <- DT::renderDataTable({
          if(input$network_type == "regul"){
            n <- seq(min(sum(!non_syn_res), MAX_VAR))
            dt = values$annotations[!non_syn_res,]
          } else {
            n <- min(sum(non_syn_res), MAX_VAR)
            dt = values$annotations[non_syn_res,]
          }
          
          DT::datatable(dt, 
                        escape = FALSE,
                        rownames = FALSE,
                        extensions = c('Scroller','Buttons'),
                        selection = list(mode = 'multiple', selected = n),
                        options = list(dom = 'Bfrtip',
                                       buttons = list('copy', 'print', list(
                                         extend = 'collection',
                                         buttons = c('csv', 'excel', 'pdf'),
                                         text = 'Download'
                                       )),
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
      if(sum(!non_syn_res) > 0) available_network_types <- c(available_network_types, "synonymous and non-coding variants" = "regul")
      if(sum(non_syn_res) > 0) available_network_types <- c(available_network_types, "non synonymous variants" = "non_syn")
      
      updateSelectInput(session = session, inputId = "network_type", choices = available_network_types)
      
      #### BUILDING DATA FOR FIRST DEFAULT PLOT ####
      my_data <- data.frame(values$global_ranges, stringsAsFactors = F)
      my_data$cdts_score <- values$res$cdts_score
      my_data$non_synonymous <- sapply(values$res$cadd.consequence, function(x) "NON_SYNONYMOUS" %in% x)
      my_data$consequences <- sapply(values$res$cadd.consequence, FUN = function(x) paste(x, collapse = ","))
      my_data$no_cdts <- is.na(my_data$cdts_score)
      my_data$selected <- TRUE
      
      
      print(input$network_type)
      
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
      
      
    })
    
    
    if(!is.null(values$res) && nrow(values$res) > 0){
      shinyBS::updateButton(session = session, inputId = "buildNetwork", 
                            disabled = FALSE, style = "success")
    } else {
      shinyBS::updateButton(session = session, inputId = "buildNetwork",
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
    
    if (!is.null(id))
      removeNotification(id)
    
    notfound <- do.call("c", values$ld[1,])
    print(notfound)
    
    #### map LD edges ####
    snv_edges <- build_snv_edges(values, "1", NULL, network_type = input$network_type)
    values$all_edges <- snv_edges
    values$current_edges <- snv_edges
    
    visNetworkProxy("my_network") %>%
      visUpdateEdges(edges = snv_edges)
    
    if(length(notfound) > 0){
      shinyBS::createAlert(session = session, anchorId = "alert_ld",
                           alertId = "alert5", title = "Data retrieval",
                           content = paste0('No LD data for the following variants: ', 
                                            paste(notfound, collapse = ', ')) , 
                           append = TRUE, style = "warning")
    } else {
      shinyBS::createAlert(session = session, anchorId = "alert_ld",
                           alertId = "alert5", title = "Data retrieval",
                           content = "LD computation succeed for all variants", 
                           append = TRUE, style = "success")
    }
    
    if(!is.null(ld_results)){
      shinyBS::updateButton(session = session, inputId = "update_ld", 
                            disabled = FALSE)
    } else {
      shinyBS::updateButton(session = session, inputId = "update_ld", 
                            disabled = TRUE)
    } 
  })
  
  #### remove LD infos ####
  observeEvent(input$removeLD, {
    
    shinyBS::updateButton(session = session, inputId = "update_ld", 
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
    shinyBS::updateButton(session = session, inputId = "buildNetwork", disabled = TRUE)
  })
  
  observeEvent(input$raw_data_rows_selected,{
    
    #setdiff(input$tbl_rows_selected, input$tbl_row_last_clicked)
    # shinyBS::updateButton(session = session, inputId = "buildNetwork", 
    #                       disabled = FALSE, style = "warning", label = "Update Network")
    
    if(length(input$raw_data_rows_selected) > MAX_VAR){
      shinyalert(title = "Selection limit", html = TRUE, type = "warning",
                 text = as.character(tags$div(style = "text-align:-webkit-center",
                                              paste0("Variants selection from ",
                                                     (length(input$raw_data_rows_selected) - MAX_VAR),
                                                     " the Network visualisation is limit (",MAX_VAR, ")")
                 )))
      shinyBS::updateButton(session = session, inputId = "buildNetwork", 
                            disabled = TRUE, style = "danger")
    } else {
      shinyBS::updateButton(session = session, inputId = "buildNetwork", 
                            disabled = FALSE, style = "success")
    }
  })
  
  buildNetwork <- eventReactive(input$buildNetwork,{
    print("buildNetwork")
    
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
    
    id <<- showNotification(paste("Building your wonderful network..."), duration = 0, type = "message")
    
    
    
    #     event.data <- event_data("plotly_selected", source = "subset")
    # 
    #     if(is.null(event.data) == T){
    #       if(nrow(my_res) > MAX_VAR){
    #         my_res <- my_res[1:MAX_VAR,]
    #       } else {
    #         my_res <- my_res
    #       }
    #     } else {
    #       my_selection <- values$my_data[subset(event.data, curveNumber == 0)$pointNumber + 1,]
    #       if(nrow(my_selection) > MAX_VAR){
    #         my_res <- my_res[my_res$query %in% my_selection$query[1:MAX_VAR],]
    #       } else {
    #         my_res <- my_res[my_res$query %in% my_selection$query,]
    #       }
    #     }
    
    values$my_res <- my_res
    
    snv_edges <- build_snv_edges(values, "1", NULL, network_type = input$network_type) #create edges without real ld info
    save(snv_edges, file = paste0(tmpDir, "/snv_edges.rda"))
    snv_nodes <- build_snv_nodes(session_values = values, network_type = input$network_type)
    save(snv_nodes, file = paste0(tmpDir, "/snv_nodes.rda"))
    
    non_null_raw_scores <- values$raw_scores
    non_null_adj_scores <- values$adjusted_scores
    
    scores_data <- build_score_nodes(session_values = values, 
                                     selected_adj_scores = non_null_adj_scores, 
                                     selected_raw_scores = non_null_raw_scores, 
                                     network_type = input$network_type)
    
    snv_nodes$infos <- scores_data$nodes_titles
    basic_ranking()
    
    values$scores_data <- scores_data
    values$current_edges <- snv_edges
    values$current_nodes <- snv_nodes
    
    meta_values_mapped <- "blue"
    
    print(paste0("network_type : ",input$network_type, " ; current_nodes : ", nrow(values$current_nodes), " ; current_edges : ", nrow(values$current_edges)))
    
    # vn <- visNetwork(values$current_nodes, 
    #                  values$current_edges) %>%
    #   visEvents(doubleClick = "function(nodes) {
    #             Shiny.onInputChange('current_node_id', nodes.nodes);
    #             ;}") %>%
    #   visInteraction(tooltipDelay = 0, hideEdgesOnDrag = T, navigationButtons = F) %>%
    #   visOptions(highlightNearest = F, clickToUse = T, manipulation = F, nodesIdSelection = TRUE) %>%
    #   #visClusteringByGroup(groups = paste0("Scores_",values$annotations$query), label = "", color = meta_values_mapped, force = T) %>%
    #   #visPhysics(solver = "forceAtlas2Based", maxVelocity = 20, forceAtlas2Based = list(gravitationalConstant = -300)) %>%
    #   visNodes(shapeProperties = list(useBorderWithImage = TRUE)) %>%
    #   visLayout(randomSeed = 2018) %>% 
    #   visPhysics(solver = "forceAtlas2Based", forceAtlas2Based = list(gravitationalConstant = -500))
    #visPhysics(solver = "forceAtlas2Based", maxVelocity = 10, forceAtlas2Based = list(gravitationalConstant = -300))
    
    # for (g in values$annotations$query){
    #   vn <- vn %>% visGroups(groupname = paste0("Scores_",g), background = "#97C2FC",
    #                          color = "#2B7CE9", shape = "square")
    # }
    
    if (!is.null(id))
      removeNotification(id)
    
    
    shinyBS::updateButton(session = session, inputId = "buildNetwork", disabled = TRUE)
    
    if(!is.null(values$current_edges) && nrow(values$current_edges) > 0){
      shinyBS::updateButton(session = session, inputId = "runLD", 
                            disabled = FALSE)
      shinyBS::updateButton(session = session, inputId = "removeLD", 
                            disabled = FALSE)
    } else {
      shinyBS::updateButton(session = session, inputId = "runLD",
                            disabled = TRUE)
      shinyBS::updateButton(session = session, inputId = "removeLD", 
                            disabled = TRUE)
    } 
    
    if(!is.null(values$current_nodes) && nrow(values$current_nodes) > 0){
      shinyBS::updateButton(session = session, inputId = "update_metascore", 
                            disabled = FALSE)
    } else {
      shinyBS::updateButton(session = session, inputId = "update_metascore",
                            disabled = TRUE)
    }
    
    return(list(nodes = values$current_nodes, edges = values$current_edges))
  })
  
  output$my_network <- renderVisNetwork({
    vn_components <- buildNetwork()
    
    if(is.null(vn_components))
      return(NULL)
    
    save(vn_components, file = paste0(tmpDir, "/vn_components.rda"))
    visNetwork(vn_components$nodes, 
               vn_components$edges) %>%
      visEvents(doubleClick = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}") %>%
      visInteraction(tooltipDelay = 0, hideEdgesOnDrag = TRUE, navigationButtons = FALSE) %>%
      visOptions(highlightNearest = TRUE, clickToUse = FALSE, manipulation = FALSE, nodesIdSelection = TRUE) %>%
      visNodes(shapeProperties = list(useBorderWithImage = TRUE)) %>%
      visLayout(randomSeed = 2018) %>% 
      visPhysics("repulsion", repulsion = list(nodeDistance = 1000))
  })
  
  #### UPDATE PLOT ####
  buildPlot <- eventReactive(c(input$fetch_annotations, input$raw_data_rows_selected), {
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
  
  buildPlot_d <- buildPlot %>% debounce(2000)
  
  output$my_plot <- renderPlotly({
    pdf(NULL) # to avoid the production of Rplots.pdf
    
    my_data <- buildPlot_d()
    if(is.null(my_data) || nrow(my_data) == 0)
      return(NULL)
    
    cdts_region_line <- values$cdts_region
    cdts_region_line$color = "gray"
    cdts_region_line$size = 1
    cdts_region_line$shape = "line-ew"
    
    text_snp <- ifelse(test = my_data$no_cdts, yes = paste0(my_data$query," (no CDTS data)"), no = my_data$query)
    text_snp <- paste(text_snp, paste0("\nConsequences: ", my_data$consequences))
    xa <- list(title = levels(my_data$seqnames),
               zeroline = FALSE,
               showline = TRUE,
               showticklabels = TRUE,
               showgrid = TRUE)
    
    plot_ly(data = my_data,
            type = "scatter",
            mode = "markers",
            x = ~start,
            y = ~cdts_score,
            symbol = ~I(shape), color = ~I(color), size = ~I(size),
            marker = list(
              opacity = 0.5
            ),
            showlegend = F, source = "subset") %>%
      add_trace(data = cdts_region_line,
                type = 'scatter',
                mode = 'lines',
                x = ~start,
                y = ~CDTS,
                hoverinfo = "none",
                line = list(color = 'gray'), 
                opacity = 0.3,
                showlegend = F) %>%
      add_trace(data = my_data,
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
      layout(xaxis = xa, dragmode =  "select")
  })
  
  #### NODES FIGURES ####
  observe({
    
    if(!is.null(values$current_nodes) && nrow(values$current_nodes) > 0){
      svn_nodes <- values$current_nodes[values$current_nodes$group == "Variants",]
      
      absolute_scores <- c('pie_bayesdel', 'pie_linsight','pie_iwscoring_known','pie_iwscoring_novel','pie_all')
      
      if(!input$snv_nodes_type %in% absolute_scores){
        if(input$update_metascore > 0)
          svn_nodes$image <- paste0("scores_figures/",input$snv_nodes_type, "_", svn_nodes$id,input$update_metascore,".png")
        else 
          svn_nodes$image <- paste0("scores_figures/",input$snv_nodes_type, "_", svn_nodes$id,".png")
      } else {
        svn_nodes$image <- paste0("scores_figures/",input$snv_nodes_type, "_", svn_nodes$id,".png")
      }
      
      values$current_nodes$image <- as.character(values$current_nodes$image)
      values$current_nodes[values$current_nodes$group == "Variants",] <- svn_nodes #update current_nodes
      
      visNetworkProxy("my_network") %>%
        visUpdateNodes(nodes = svn_nodes)
    }
    
  })
  
  #### SCORES ####
  observe({
    # non_null_raw_scores <- names(values$raw_scores[!sapply(values$raw_scores, is.null)])
    # non_null_adj_scores <- names(values$adjusted_scores[!sapply(values$adjusted_scores, is.null)])
    non_null_raw_scores <- values$raw_scores
    non_null_adj_scores <- values$adjusted_scores
    non_null_raw_scores_pretty_names <- gsub(x = non_null_raw_scores, 
                                             pattern = "score|rawscore", replacement = "")
    non_null_raw_scores_pretty_names <- gsub(x = non_null_raw_scores_pretty_names,
                                             replacement = "", pattern = "\\.$")
    non_null_raw_scores_pretty_names <- gsub(x = non_null_raw_scores_pretty_names, 
                                             replacement = " ", pattern = "\\.")
    non_null_raw_scores_pretty_names <- gsub(x = non_null_raw_scores_pretty_names, 
                                             pattern = "_(.*)_", replacement = " (\\1)")
    
    non_null_adj_scores_pretty_names <- gsub(x = non_null_adj_scores, 
                                             pattern = "score|rankscore", replacement = "")
    non_null_adj_scores_pretty_names <- gsub(x = non_null_adj_scores_pretty_names,
                                             replacement = "", pattern = "\\.$")
    non_null_adj_scores_pretty_names <- gsub(x = non_null_adj_scores_pretty_names, 
                                             replacement = " ", pattern = "\\.")
    non_null_adj_scores_pretty_names <- gsub(x = non_null_adj_scores_pretty_names, 
                                             pattern = "_(.*)_", replacement = " (\\1)")
    
    switch(input$network_type,
           non_syn = {
             predictors <- non_null_adj_scores
           },
           regul = {
             predictors <- non_null_raw_scores
           }
    )
    
    if(length(non_null_raw_scores) > 0 || length(non_null_adj_scores) > 0){
      updateSelectizeInput(session, "selected_scores",
                           choice = predictors,
                           selected = predictors)
    }
  })
  
  
  #### METASCORES ####
  observeEvent(input$update_metascore, {
    
    scores_data <- build_score_nodes(values,
                                     selected_adj_scores = input$selected_scores, 
                                     selected_raw_scores = input$selected_scores, 
                                     inc = input$update_metascore, network_type = input$network_type)
    save(scores_data, file = paste0(tmpDir, "/scores_data.rda"))
    values$scores_data <- scores_data
    
    basic_ranking(inc = input$update_metascore)
    
    svn_nodes <- values$current_nodes[values$current_nodes$group == "Variants",]
    svn_nodes$image <- paste0("scores_figures/",input$snv_nodes_type, "_", 
                              svn_nodes$id,input$update_metascore,".png")
    
    visNetworkProxy("my_network") %>% 
      visUpdateNodes(nodes = svn_nodes)
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
  
  #### SCORES SUBGRAPH ####
  observeEvent(input$current_node_id, {
    cat("1",input$current_node_id,"\n")
    values$selected_node <- input$current_node_id
    snv_score_details <- values$scores_data$nodes[values$scores_data$nodes$group == paste0("Scores_",values$selected_node), c("id","color","label")]
    if(nrow(snv_score_details) > 0){
      color_code <- as.character(snv_score_details$color)
      text_color_code <- sapply(X = color_code, FUN = contrasting_text_color)
      snv_score_details$id <- gsub(x = snv_score_details$id, pattern = paste0("_", values$selected_node), replacement = "")
      snv_score_details <- snv_score_details[,-2] #suppression de la colonne 'color' 
      colnames(snv_score_details) <- c("predictors", "value")
      rownames(snv_score_details) <- NULL
      snv_score_details$value <- round(x = as.numeric(as.character(snv_score_details$value)),digits = 4)
      
      tab <- tableHTML(snv_score_details, theme = 'scientific', rownames = FALSE)
      for (i in 1:nrow(snv_score_details)) { 
        tab <- tab %>% add_css_row(css = list(c('background-color','font-weight','color'), 
                                              c(color_code[i],'bold', text_color_code[i])), rows = i+1) # first column is header
      }
      tab <- tab %>% add_css_column(css = list(c('background-color','text-align','font-weight','color'),
                                               c('white','left','normal','gray')), columns = 1)
    } else {
      tab <- "No data for the selected predictors"
    }
    shinyalert(title = input$current_node_id, html = TRUE,
               text = as.character(tags$div(style = "text-align:-webkit-center",
                                            tab)))
  })
  
  
  observeEvent(input$select_button, {
    selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    print(paste('click on',selectedRow))
  })
  
  output$scores_description <- renderUI({
    X <- readr::read_tsv(file = paste0(dataDir, 'scores_description.tsv'))
    X <- split(X, X$group_name)
    
    table_content <- list()
    for(n in names(X)){
      score_infos <- X[[n]]
      score_infos <- score_infos[order(score_infos$name),]
      integration_list <- list()
      
      for(r in 1:nrow(score_infos)){ 
        integration_list[[length(integration_list) + 1]] <- tags$li(
          paste0(score_infos[r,]$name,
                 " (id: ", score_infos[r,]$id,
                 "; scope:" , score_infos[r,]$scope, 
                 "; orientation: ", score_infos[r,]$orientation, ")")
        )
      }
      
      table_content[[length(table_content) + 1]] <- p(
        strong(n),
        br(),
        helpText(paste0(gsub(x = score_infos$description[1], pattern = ".$", replacement = ""), " (",score_infos$reference[1],").")), 
        "Integration",
        tags$div(style = "margin-left: 30px;",
                 tags$ul(
                   integration_list
                 )
        ),
        hr()
      )
    }
    
    tagList(
      
      table_content
      
      
    )
  })
  
  #### ON STOP ####
  onStop(function() { 
    old_figures <- dir(path = "www/scores_figures/", full.names = T)
    file.remove(old_figures)
    temp_files <- dir(path = tmpDir, full.names = T, recursive = T)
    file.remove(temp_files)
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



