library(shiny)
library(d3heatmap)
library(visNetwork)
library(myvariant)
library(shinyBS)

options(shiny.trace = FALSE)

server <- function(input, output, session) {
  source('helper.R', local = TRUE)
  preload <- list(
    locus_80 = "rs2660776\nrs2660774\nrs2575792\nrs2660773\nrs2660772\nrs1370042\nrs1160076\nrs2018075\nrs62259569\nrs2575799\nrs2253915\nrs2253916\nrs2253921\nrs2575802\nrs138876277\nrs1370044\nrs2575805\nrs144055118\nrs35711034\nrs35733956\nrs74520234\nrs147080240\nrs4859022\nrs2262594\nrs2254904\nrs729942\nrs1370045\nrs35870472\nrs34339032\nrs2660785\nrs2575815\nrs2575817\nrs2660787\nrs2575818\nrs2660788\nrs2118082\nrs2244125\nrs2660790\nrs2660791\nrs2660792\nrs2660794\nrs61340205\nrs1821751\nrs2660796\nrs2196523\nrs1370040\nrs2660797\nrs2575783\nrs2660798\nrs2575780\n3:87044690:A:AAT\nrs2575775\nrs2660799\nrs2575774\nrs2575772\nrs2043663\nrs13066793\nrs1960268\nrs1025631\nrs4859107\nrs201414876",
    locus_78 = "rs73167030\nrs12159970\nrs73167031\nrs56150793\nrs6001911\nrs74767555\nrs111843590\nrs6001912\nrs10483203\nrs56283550\nrs12158872\nrs17001907\n22:40836155:CAAAA:C\nrs56215843\nrs61675795\nrs73167042\nrs6001915\nrs17001915\nrs73167045\nrs73167052\nrs5995856\nrs17001920\nrs6001920\nrs73167053\nrs5995860\nrs28419341\nrs28552449\nrs112375848\nrs55864500\nrs183438976\nrs73167058\nrs5995862\nrs17001943\nrs58778028\nrs73167063\nrs5995864\nrs75523053\nrs74486969\nrs61211547\nrs12159787\nrs10483204\nrs145014115\nrs4299422\nrs2899337\nrs201540223\nrs73167066\nrs71777635\nrs73167067\nrs6001930\nrs17001974\nrs6001931\nrs6001932\nrs73167069\nrs73167072\nrs73167073\nrs73167076\nrs17001977\nrs3827381\nrs3827382\nrs73167079\nrs10483205\nrs112880707\nrs73167080\nrs73167082\nrs113409089\nrs74278065\nrs6001935\nrs56135013\nrs6001937\nrs73167089\nrs113798157\nrs5995867\nrs73167090\nrs77426923\nrs56182212\nrs73167092\nrs199614224\nrs6001939\nrs73167093\nrs6001942\nrs138175438\nrs73167096\nrs73167097\nrs73167098\nrs373996442\nrs113966362\nrs73167101\nrs17001993\nrs17001994\nrs73169003\nrs10483206\nrs6001946\nrs6001949\nrs66987842\nrs17001997\nrs6001950\nrs56124626\nrs141580207\nrs112034576\nrs55993771\nrs142270009\nrs6001954\nrs5995870\nrs57693796\nrs5995871\nrs139903991\nrs55849114\nrs55722862\nrs55960299\nrs56334449\nrs55662398\nrs6001962\nrs150091203\nrs932379\nrs73169026\nrs73169028\nrs73169032\nrs55674424\nrs5995875\nrs6001965\nrs56092708\nrs16985899\nrs17002019\nrs17002020\nrs201432472\nrs73169036\nrs12160047\nrs56047425\nrs56411183\nrs73169040\nrs117423948\nrs5995876\nrs138476044\nrs146829921\nrs73169043\nrs116973859\nrs28444678\nrs375665020\nrs17002024\nrs57705902\nrs79799751\nrs6001973\nrs6001974\nrs73169051\nrs143586053\nrs17002026\nrs17002027\nrs28567076\nrs73169056\nrs55853923\nrs17002030\nrs73169057\nrs73169058\nrs73169060\nrs141597790\nrs142127948\nrs17002034\nrs17002036\nrs6001979\nrs17002038\nrs6001980\nrs55719110\nrs55997921\nrs55775031\nrs139645469\nrs73169065\nrs59047419\nrs55669535\nrs192603101\nrs73169071\nrs5995881\nrs73169072\nrs6001981\nrs145351698\nrs73169077\nrs73169078\nrs6001982\nrs73169083\nrs73169084\nrs112497956\nrs6001983\nrs718193\nrs56108505\nrs73169087\nrs73169089\nrs6001984\nrs71785733\nrs73169091\nrs117669949\nrs78694459\nrs73169096\nrs73169097\nrs145005064\nrs73169098\nrs145624734\nrs56051217",
    locus_70 = "rs4889891\nrs9905914\nrs9896202\nrs745571\nrs745570\nrs2587505\nrs139427260\nrs71161686\nrs2587507",
    locus_60 = "13:32866927:CAATAAATAAATA:CAATAAATA\nrs56084662\nrs11571815\nrs11571818\nrs11571833",
    locus_38 = "rs62485509\nrs720475",
    locus_6 = "rs181595184\nrs12129763\nrs6686987\nrs55899544\nrs6427943\nrs6678914\nrs6703244\nrs12028423\nrs12131882\nrs12132085\nrs12129456\nrs12129536\nrs10920365\nrs2167588\nrs4950774\nrs12032424\nrs4950775\nrs4950836\nrs6143572\nrs896548\nrs4245706\nrs12032080\nrs3795598\nrs12143329\nrs4950837\nrs12021815\nrs60573451\nrs140473199",
    locus_2 = "rs11102694\nrs2358994\nrs2358995\nrs7547478\nrs7513707\nrs12022378\nrs3761936\nrs11102701\nrs11102702\nrs12046289\nrs112974454")
  
  #### CONF ####
  tmpDir <- 'temp/'
  resultDir <- 'temp/'
  
  app.conf <- list(TABIX = '/usr/local/bin/tabix',
                   VCF = '/Users/nekomimi/Workspace/vexor/vexor/data/1000Genomes/')
  path_to_images <- "~/Workspace/DSNetwork/www/scores_figures/"
  
  values <- reactiveValues()
  values$maf_infos <- as.matrix(data.frame(waiting = ""))
  values$annotations <- as.matrix(data.frame(waiting = ""))
  values$all_edges <- NULL
  values$all_nodes <- NULL
  values$current_edges <- NULL
  values$current_nodes <- NULL
  values$scores_data <- NULL
  values$selected_node <- ''
  values$ld_regions <- NULL
  load('data/scores_correlation_matrice.rda')
  
  # observeEvent(input$query, ({
  #   
  #   if(nchar(input$query) > 0){
  #     shinyBS::updateButton(session = session, inputId = "fetch_annotations", 
  #                           disabled = FALSE)
  #   } else {
  #     shinyBS::updateButton(session = session, inputId = "fetch_annotations", 
  #                           disabled = TRUE)
  #   }
  # }))
  
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
            formatSingleHgvs(as.numeric(x[1]), as.numeric(x[2]), x[3], x[4])
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
            formatSingleHgvs(as.numeric(x[1]), as.numeric(x[2]), x[3], x[4])
          }
        }
      })
      return(modstring)
    }
  }
  
  query_control <- eventReactive(input$fetch_annotations, {
    shinyBS::updateButton(session = session, inputId = "fetch_annotations", 
                          disabled = TRUE)
    modstring <- c()
    if(!is.null(input$query_file))
      modstring <- c(modstring, transform_query_file(input$query_file$datapath))
    
    if(nchar(input$query) > 0)
      modstring <- c(modstring, transform_query(input$query))
    
    fail_transfo <- names(modstring[grep(x = modstring, pattern = 'FAIL')])
    if(length(modstring) > 0){
      if(length(fail_transfo) > 0){
        res <- paste0('Id recognition fails for: ', paste(fail_transfo, collapse = ","), '.')
      }
    }            
    
  })
  
  output$transform_res <- renderText({ query_control() })
  
  #### FETCH ANNOTATIONS ####
  observeEvent(input$fetch_annotations, {
    values$maf_infos <- as.matrix(data.frame(waiting = ""))
    values$annotations <- as.matrix(data.frame(waiting = ""))
    
    modstring <- c()
    if(!is.null(input$query_file))
      modstring <- c(modstring, transform_query_file(input$query_file$datapath))
    
    if(nchar(input$query) > 0)
      modstring <- c(modstring, transform_query(input$query))
    
    valid_transfo <- unique(modstring[!grepl(x = modstring, pattern = 'FAIL')])
    
    if(length(valid_transfo) > 0){
      #### run myvariant ####
      res <- as.data.frame(getVariants(hgvsids = valid_transfo,
                                       verbose = F, return.as = "DataFrame",
                                       fields = c("dbsnp","cadd","dbnsfp")), 
                           stringsAsFactors = TRUE)
      res$linsight <- NA
      #save(res, file = 'objects/res.rda')
      values$res <- res
      
      if(!is.null(values$res) && nrow(values$res) > 0){
        
        #### suppression des doublons ####
        duplicated_entries <- values$res$query[duplicated(res$query)]
        row_2_delete <- c() 
        if(length(duplicated_entries) > 0){
          for(duplicated_entry in duplicated_entries){
            b <- values$res[values$res$query == duplicated_entry,]
            minimum_na_row <- which.min(apply(b, 1, function(x) sum(is.na(x))))
            row_2_delete <- c(row_2_delete, 
                              as.numeric(row.names(b)[! row.names(b) %in% names(minimum_na_row)]))
          }
        }
        if(length(row_2_delete) > 0){
          res <- data.frame(values$res[-row_2_delete,], row.names = NULL)
          values$res <- res
          #save(res, file = 'objets/res.rda')
        }
        
        
        #### ranges ####
        my_ranges <- apply(X = values$res, MARGIN = 1, 
                           FUN = function(x) hgvsToGRange(hgvs_id = x[names(x) == "X_id"]$X_id, 
                                                          query_id = x[names(x) == "query"]$query))
        my_ranges <- do.call("c", my_ranges) 
        my_ranges <- sort(unique(my_ranges))
        #save(my_ranges, file = "objets/my_ranges.rda")
        values$my_ranges <- my_ranges
        
        #### run linsight ####
        requested_chromosomes <- seqlevelsInUse(my_ranges)
        
        for(requested_chr in requested_chromosomes){
          if(file.exists(paste0('data/LINSIGHT_',requested_chr,'.rda'))){
            load(paste0('data/LINSIGHT_',requested_chr,'.rda'))
            hits <- findOverlaps(query = my_ranges, subject = gr)
            for(i in seq_along(hits)){ 
              hit <- hits[i]
              query <- my_ranges[queryHits(hit)]$query
              value <- gr[subjectHits(hit)]$score
              values$res[values$res$query == query,"linsight"] <- value
            }
          }
        }
        if(sum(is.na(values$res$linsight)) == length(values$res$linsight)){ # no lINSIGHT RESULTS
          values$res$linsight <- NULL
        }
        
        #values$res <- res
        #save(res, file = 'objets/res.rda')
        
        #### population freq ####
        maf_infos <- values$res[,grepl(x = colnames(values$res), pattern = "query|cadd.ref|cadd.alt|cadd.1000g|maf")] 
        colnames(maf_infos) <- gsub(x = colnames(maf_infos), pattern = "cadd.1000g.(.*)", replacement = "MAF_\\1_1000G")
        row.names(maf_infos) <- NULL
        values$maf_infos <- as.matrix(maf_infos)
        
        #### raw annotations ####
        annotations_infos <- data.frame(values$res, stringsAsFactors = FALSE)
        
        # supprimer not found
        if("notfound" %in% colnames(annotations_infos)){
          annotations_infos <- annotations_infos[is.na(annotations_infos$notfound),]
          annotations_infos$notfound <- NULL
        }
        
        # supprimer licence
        annotations_infos <- annotations_infos[,!grepl(x = colnames(annotations_infos), pattern = "license")]
        
        # suppriner un niveau de nested
        annotations_infos <- data.frame(lapply(annotations_infos, 
                                               function(x) unlist(lapply(x, paste, collapse = ","))))
        # suppriner NA column
        # annotations_infos <- annotations_infos[apply(X = annotations_infos, MARGIN = 1, FUN = function(x) sum(is.na(x)) != (length(x) - 1)),]
        
        values$annotations <- annotations_infos
        #save(annotations_infos, file = "objets/annotations.rda")
        
        adjusted_scores <- list(
          Polyphen2.HDIV.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.polyphen2.hdiv.rankscore"),
          Polyphen2.HVAR.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.polyphen2.hvar.rankscore"),
          Sift.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.sift.converted_rankscore"),
          MutationTaster.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.mutationtaster.converted_rankscore"),
          LRT.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.lrt.converted_rankscore"),
          MutationAssessor.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.mutationassessor.rankscore"),
          FATHMM.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.fathmm.rankscore"),
          FATHMM.MKL.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.fathmm.mkl.coding_rankscore"),
          Provean.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.provean.rankscore"),
          VEST3.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.vest3.rankscore"),
          MetaSVM.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.metasvm.rankscore"),
          MetaLR.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.metalr.rankscore"),
          CADD.phredscore = extract_score_and_convert(annotations_infos, "cadd.phred"),
          DANN.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.dann.rankscore"),
          fitCons.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.integrated.fitcons_rankscore"),
          SiPhy.logodds.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.siphy_29way.logodds_rankscore"),
          GERP.element.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.gerp...rs_rankscore"),
          PhyloP.20way.mammalian.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.phylo.p20way.mammalian_rankscore"),
          PhyloP.100way.vertebrate.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.phylo.p100way.vertebrate_rankscore"),
          phastCons.20way.mammalian.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.phastcons.20way.mammalian_rankscore"),
          phastCons.100way.vertebrate.rankscore = extract_score_and_convert(annotations_infos, 
                                                                            "dbnsfp.phastcons.100way.vertebrate_rankscore"),
          Eigen.phredscore = extract_score_and_convert(annotations_infos, "dbnsfp.eigen.phred"),
          Eigen.PC.phredscore = extract_score_and_convert(annotations_infos, "dbnsfp.eigen.pc.phred"),
          Eigen.PC.rankscore = extract_score_and_convert(annotations_infos, "dbnsfp.eigen.pc.raw_rankscore")
        )
        
        raw_scores <- list(
          Polyphen2.HDIV.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.polyphen2.hdiv.score"),
          Polyphen2.HVAR.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.polyphen2.hvar.score"),
          Polyphen2.rawscore = extract_score_and_convert(annotations_infos, "cadd.polyphen.val"),
          Sift.rawscore = extract_score_and_convert(annotations_infos, "cadd.sift.val"),
          MutationTaster.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.mutationtaster.score"),
          LRT.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.lrt.score"),
          MutationAssessor.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.mutationassessor.score"),
          FATHMM.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.fathmm.score"),
          FATHMM.MKL.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.fathmm.mkl.coding_score"),
          Provean.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.provean.score"),
          VEST3.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.vest3.score"),
          MetaSVM.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.metasvm.score"),
          metaLR.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.metalr.score"),
          CADD.rawscore = extract_score_and_convert(annotations_infos, "cadd.rawscore"),
          DANN.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.dann.score"),
          fitCons.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.integrated.fitcons_score", "cadd.fitcons"), # get cadd.fitcons if missing
          SiPhy.logodds.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.siphy_29way.logodds"),
          GERP.neutral.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.gerp...nr", "cadd.gerp.n"), # get cadd.gerp.n if missing
          GERP.element.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.gerp...rs", "cadd.gerp.rs"), # get cadd.gerp.rs if missing
          
          PhyloP.mammalian_non.human_.rawscore = extract_score_and_convert(annotations_infos, "cadd.phylop.mammalian"),
          PhyloP.primate_non.human_.rawscore = extract_score_and_convert(annotations_infos, "cadd.phylop.primate"),
          PhyloP.vertebrate_non.human_.rawscore = extract_score_and_convert(annotations_infos, "cadd.phylop.vertebrate"),
          PhyloP.20way.mammalian.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.phylo.p20way.mammalian"),
          PhyloP.100way.vertebrate.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.phylo.p100way.vertebrate"),
          phastCons.mammalian_non.human_.rawscore = extract_score_and_convert(annotations_infos, "cadd.phast_cons.mammalian"),
          phastCons.primate_non.human_.rawscore = extract_score_and_convert(annotations_infos, "cadd.phast_cons.primate"),
          phastCons.vertebrate_non.human_.rawscore = extract_score_and_convert(annotations_infos, "cadd.phast_cons.vertebrate"),
          phastCons.20way.mammalian.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.phastcons.20way.mammalian"),
          phastCons.100way.vertebrate.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.phastcons.100way.vertebrate"),
          Eigen.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.eigen.raw"),
          Eigen.PC.rawscore = extract_score_and_convert(annotations_infos, "dbnsfp.eigen.pc.raw"),
          LINSIGHT.rawscore = extract_score_and_convert(annotations_infos, "linsight")
        )
        # save(raw_scores, adjusted_scores, file = "objets/scores.rda")
        values$raw_scores <- raw_scores
        values$adjusted_scores <- adjusted_scores
      }
    }
    
    if(!is.null(values$res)){
      shinyBS::updateButton(session = session, inputId = "runLD", 
                            disabled = FALSE)
    } else {
      shinyBS::updateButton(session = session, inputId = "runLD", 
                            disabled = TRUE)
    }
  })
  
  output$query_res <- renderText({ 
    if(!is.null(values$res) && nrow(values$res) > 0){
      
      if("notfound" %in% colnames(values$res)){
        notfound_id <- values$res[!is.na(values$res$notfound),"query"]
      } else {
        notfound_id <- NULL
      }
      
      ifelse(test = length(notfound_id) > 0, 
             yes = paste0("No Annotations for the following variants: ", 
                          paste(notfound_id, collapse = ","),"."), 
             no = "Annotations found for all variants")  
    }
  })
  
  output$populations <- DT::renderDataTable({
    DT::datatable(values$maf_infos,
                  options = list(scrollX = TRUE))
  })
  
  output$raw_data <- DT::renderDataTable({
    DT::datatable(as.matrix(values$annotations),
                  options = list(scrollX = TRUE))
  })
  
  output$downloadRawTable <- downloadHandler(
    filename = "annotations_table.csv",
    content = function(file) {
      write.csv(as.matrix(values$annotations), file)
    }
  )
  
  output$downloadFreqTable <- downloadHandler(
    filename = "frequencies_table.csv",
    content = function(file) {
      write.csv(values$maf_infos, file)
    }
  )
  
  #### FETCH LD INFOS ####
  observeEvent(input$runLD, {
    updateTabsetPanel(session, "results_tabset",
                      selected = "ld_results"
    )
  })
  
  runLD <- eventReactive(input$runLD,{
    
    shinyBS::updateButton(session = session, inputId = "runLD", 
                          disabled = TRUE)
    
    my_ranges <- values$my_ranges
    
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
                                 results_dir = resultDir, 
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
    
    regions <- unlist(sapply(ilets, function(selection){
      current_range  <- setStudyRange(my_ranges, selection)
      region <- setStudyRegion(current_range)
      region <- gsub(x = region, pattern = 'chr', replacement = '')
      region <- gsub(x = region, pattern = ':', replacement = '_')
    }))
    
    values$ld <- ld_results
    values$ld_regions <- regions
    
    #save(regions, file = "objets/regions.rda")
    
    if(!is.null(values$ld)){
      shinyBS::updateButton(session = session, inputId = "buildNetwork", 
                            disabled = FALSE)
      updateSelectInput(session, "ld_regions", choices = as.character(regions))
    } else {
      shinyBS::updateButton(session = session, inputId = "buildNetwork", 
                            disabled = TRUE)
      updateSelectInput(session, "ld_regions",choices = c())
    }
  })
  
  output$runld_res <- renderText({ 
    if(!is.null(values$res) && nrow(values$res) > 0){
      runLD()
      notfound <- do.call("c", values$ld[1,])
      
      return(ifelse(test = length(notfound) > 0, 
                    yes = paste0('No LD data for the following variants: ', paste(notfound, collapse = ', ')),
                    no = "LD computation succeed for all variants"))  
    }
  })
  
  output$ld_plot <- renderImage({  
    filename <- "processing.gif"
    if(!is.null(values$ld)){
      filename <- paste0(resultDir,"/ld_figures/LD_plot_",input$ld_regions,".png")
    }
    
    list(src = filename,
         alt = "This is alternate text")
    
  },deleteFile = FALSE)
  
  #### BUILDING NETWORK ####
  observeEvent(input$buildNetwork, {
    updateTabsetPanel(session, "results_tabset",
                      selected = "network_results"
    )
  })
  
  buildNetwork <- eventReactive(input$buildNetwork,{
    
    shinyBS::updateButton(session = session, inputId = "buildNetwork", 
                          disabled = TRUE)
    
    snv_dist_edges <- build_snv_edges(values, "0", 1000) #default
    snv_ld_edges <- build_snv_edges(values, "1", NULL) #default
    snv_edges <- rbind(snv_dist_edges, snv_ld_edges)
    snv_nodes <- build_snv_nodes(values)
    non_null_raw_scores <- names(values$raw_scores[!sapply(values$raw_scores, is.null)])
    non_null_adj_scores <- names(values$adjusted_scores[!sapply(values$adjusted_scores, is.null)])
    scores_data <- build_score_nodes(values,selected_adj_scores = non_null_adj_scores, 
                                     selected_raw_scores = non_null_raw_scores)
    basic_ranking()
    
    values$scores_data <- scores_data
    values$all_nodes <- snv_nodes
    values$all_edges <- snv_edges
    
    all_ne <- list(nodes = values$all_nodes,
                   edges = values$all_edges)
    
    # suppress unwanted edges
    current_edges <- values$all_edges
    dist_edges_idx <- which(current_edges$type == "snv_dist_edges")
    current_edges <- current_edges[-dist_edges_idx,]
    
    min_ld <- 0.8
    ld_edges_idx <- which((current_edges$type == "snv_ld_edges") & (as.numeric(current_edges$xvalue) < min_ld))
    current_edges <- current_edges[-ld_edges_idx,]
    
    values$current_edges <- current_edges
    values$current_nodes <- values$all_nodes
    
    cur_ne <- list(nodes = values$current_nodes,
                   edges = values$current_edges)
    
    #save(all_ne, cur_ne, file = "objets/values.rda")
    
    meta_values_mapped <- "blue"
    vn <- visNetwork(values$current_nodes, 
                     values$current_edges, width = "auto", height = "500px") %>%
      visEvents(doubleClick = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}") %>%
      visInteraction(tooltipDelay = 0, hideEdgesOnDrag = T) %>%
      # visOptions(highlightNearest = F, clickToUse = F, manipulation = F) %>%
      visClusteringByGroup(groups = paste0("Scores_",values$annotations$query), label = "", color = meta_values_mapped, force = T) %>%
      visPhysics(solver = "forceAtlas2Based", maxVelocity = 20, forceAtlas2Based = list(gravitationalConstant = -300)) %>%
      visNodes(shapeProperties = list(useBorderWithImage = TRUE)) %>%
      visLayout(randomSeed = 2018)
    # visPhysics(solver = "forceAtlas2Based", maxVelocity = 10, forceAtlas2Based = list(gravitationalConstant = -300))
    
    for (g in values$annotations$query){
      vn <- vn %>% visGroups(groupname = paste0("Scores_",g), background = "#97C2FC",
                             color = "#2B7CE9", shape = "square")
    }
    
    saveRDS(object = vn, file = "objets/init_vn.rds")
    return(vn)
  })
  
  output$my_network <- renderVisNetwork({
    buildNetwork()
  })
  
  #### UPDATE NETWORK ####
  
  #### NODES FIGURES ####
  observe({
    
    if(!is.null(values$current_nodes)){
      svn_nodes <- values$current_nodes[values$current_nodes$group == "Variants",]
      if(input$update_metascore > 0)
        svn_nodes$image <- paste0("scores_figures/",input$snv_nodes_type, "_", svn_nodes$id,input$update_metascore,".png")
      else 
        svn_nodes$image <- paste0("scores_figures/",input$snv_nodes_type, "_", svn_nodes$id,".png")
      
      values$current_nodes$image <- as.character(values$current_nodes$image)
      values$current_nodes[values$current_nodes$group == "Variants",] <- svn_nodes #update current_nodes
      
      visNetworkProxy("my_network") %>%
        visUpdateNodes(nodes = svn_nodes)
    }
    
  })
  
  
  #### LD ####
  observeEvent(input$update_ld, {
    # supprimer les edges de type dist
    dist_edges_id_2_remove <- values$current_edges[values$current_edges$type == "snv_dist_edges",]$id
    
    # extraire les bon edges
    max_ld <- input$ld_range[2]
    min_ld <- input$ld_range[1]
    
    # supprimer les edges de type ld qui ne sont pas dans le range
    ld_edges_id_2_remove <- values$current_edges[(values$current_edges$type == "snv_ld_edges") & 
                                                   ((as.numeric(values$current_edges$xvalue) < min_ld | as.numeric(values$current_edges$xvalue) > max_ld)),]$id
    
    # ajouter les edges de type ld qui sont dans le range mais qui n'étaient pas dans le current_edges
    ld_edges <- values$all_edges[values$all_edges$type == "snv_ld_edges",]
    ld_edges_id_2_add <- ld_edges[(as.numeric(ld_edges$xvalue) >= min_ld) & (as.numeric(ld_edges$xvalue) <= max_ld),]
    
    edges_id_2_remove <- c(dist_edges_id_2_remove, ld_edges_id_2_remove)
    
    edges_2_keep <- values$current_edges[! values$current_edges$id %in% edges_id_2_remove,]
    edges_2_keep <- rbind(edges_2_keep, ld_edges_id_2_add)
    
    values$current_edges <- edges_2_keep
    
    visNetworkProxy("my_network") %>%
      visUpdateEdges(edges = edges_2_keep)
    
    visNetworkProxy("my_network") %>%
      visRemoveEdges(id = edges_id_2_remove)
  })
  
  #### DISTANCE ####
  observeEvent(input$update_dist, {
    # supprimer les edges de type ld
    ld_edges_id_2_remove <- values$current_edges[values$current_edges$type == "snv_ld_edges",]$id
    
    # supprimer les edges entre les nodes trop eloignés
    max_dist <- 1000 * input$dist_range
    dist_edges_id_2_remove <- values$current_edges[(values$current_edges$type == "snv_dist_edges") & (as.numeric(values$current_edges$xvalue) > max_dist),]$id
    
    #dist_edges_id_2_remove <- dist_edges_id_2_remove[as.numeric(dist_edges_id_2_remove$xvalue) > max_dist,]$id
    
    # ajouter les edges de type dist qui sont dans le range mais qui n'étaient pas dans le current_edges
    dist_edges <- values$all_edges[values$all_edges$type == "snv_dist_edges",]
    dist_edges_2_add <- dist_edges[(as.numeric(dist_edges$xvalue) <= max_dist),]
    
    edges_id_2_remove <- c(ld_edges_id_2_remove, dist_edges_id_2_remove)
    
    # edges a garder
    edges_2_keep <- values$current_edges[! values$current_edges$id %in% edges_id_2_remove,]
    edges_2_keep <- rbind(edges_2_keep, dist_edges_2_add)
    
    values$current_edges <- edges_2_keep
    
    visNetworkProxy("my_network") %>%
      visUpdateEdges(edges = edges_2_keep)
    
    visNetworkProxy("my_network") %>%
      visRemoveEdges(id = edges_id_2_remove)
  })
  
  #### FOCUS ####
  observe({
    updateSelectInput(session, "focus",
                      choices = c("None", as.character(values$current_nodes$id)),
                      selected = "None"
    )
  })
  
  observe({
    if(input$focus != "None"){
      visNetworkProxy("my_network") %>%
        visFocus(id = input$focus, scale = 1)
    } else {
      visNetworkProxy("my_network") %>%
        visFit()
    }
  })
  
  #### HIGHLIGHT ####
  observe({
    if("cadd.consdetail" %in% colnames(values$annotations)){
      choices <- unique(unlist(strsplit(x = as.character(values$annotations$cadd.consdetail), split = ","))) 
      updateSelectInput(session, "highlight",
                        choices = c("None", choices),
                        selected = "None")
    }
  })
  
  observe({
    svn_nodes <- values$current_nodes[values$current_nodes$group == "Variants",]
    
    # remise à zero
    svn_nodes$shape <- "circularImage"
    
    if(input$highlight != "None"){
      # nodes to highlight
      svn_id_2_highlight <- as.character(values$annotations$query)[grepl(x = as.character(values$annotations$cadd.consdetail),
                                                                         pattern = input$highlight)]
      svn_nodes[svn_nodes$id %in% svn_id_2_highlight,]$shape <- "image"
      values$current_nodes$shape <- as.character(values$current_nodes$shape)
      values$current_nodes[values$current_nodes$group == "Variants",] <- svn_nodes #update current_nodes
    }
    
    visNetworkProxy("my_network") %>%
      visUpdateNodes(nodes = svn_nodes)
  })
  
  #### SCORES ####
  observe({
    non_null_raw_scores <- names(values$raw_scores[!sapply(values$raw_scores, is.null)])
    non_null_adj_scores <- names(values$adjusted_scores[!sapply(values$adjusted_scores, is.null)])
    
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
    
    
    
    if(length(non_null_raw_scores) > 0){
      names(non_null_raw_scores) <- non_null_raw_scores_pretty_names
      updateSelectizeInput(session, "raw_scores",
                           choice = non_null_raw_scores,
                           selected = non_null_raw_scores)
    }
    
    if(length(non_null_adj_scores) > 0){
      names(non_null_adj_scores) <- non_null_adj_scores_pretty_names
      updateSelectizeInput(session, "adj_scores",
                           choice = non_null_adj_scores,
                           selected = non_null_adj_scores)
    }
    
    if(length(non_null_raw_scores) > 0 || length(non_null_adj_scores) > 0){
      updateSelectizeInput(session, "predictors",
                           choice = list("Raw scores" = non_null_raw_scores,
                                         "Adjusted scores" = non_null_adj_scores),
                           selected = c())
    }
    
  })
  
  #### SCORES STATS ####
  
  compute_scores_missing_data <- eventReactive(input$get_predictors_info, {
    if(length(input$predictors) > 0){
      non_null_raw_scores <- values$raw_scores[!sapply(values$raw_scores, is.null)]
      non_null_adj_scores <- values$adjusted_scores[!sapply(values$adjusted_scores, is.null)]
      
      all_scores <- c(non_null_adj_scores, non_null_raw_scores)
      all_scores <- all_scores[input$predictors]
      all_scores <- data.frame(all_scores)
      
      scores_stats <- do.call("rbind",
                              lapply(X = all_scores,
                                     FUN = function(x) {
                                       data.frame(Missing_values = sum(is.na(x)),
                                                  Min = signif(x = min(x, na.rm = T), digits = 3),
                                                  #First_Qt = quantile(x, probs = 0.25, names = F, na.rm = T),
                                                  Mean = signif(x = mean(x, na.rm = TRUE), digits = 3),
                                                  #Median = median(x, na.rm = TRUE),
                                                  #Third_Qt =  quantile(x, probs = 0.75, names = F, na.rm = T),
                                                  Max = signif(x = max(x, na.rm = TRUE), digits = 3))
                                     }
                              )
      )
      
      pretty_rownames <- row.names(scores_stats)
      pretty_rownames <- gsub(x = pretty_rownames, pattern = "score|rawscore|rankscore", replacement = "")
      pretty_rownames <- gsub(x = pretty_rownames, replacement = "", pattern = "\\.$")
      pretty_rownames <- gsub(x = pretty_rownames, replacement = " ", pattern = "\\.")
      pretty_rownames <- gsub(x = pretty_rownames, pattern = "_(.*)_", replacement = " (\\1)")
      
      row.names(scores_stats) <- pretty_rownames
      
      return(scores_stats)
    }
  })
  
  compute_scores_stats <- eventReactive(input$get_predictors_info,{
    
    non_null_raw_scores <- values$raw_scores[!sapply(values$raw_scores, is.null)]
    non_null_adj_scores <- values$adjusted_scores[!sapply(values$adjusted_scores, is.null)]
    
    all_scores <- c(non_null_adj_scores, non_null_raw_scores)
    all_scores <- all_scores[input$predictors]
    all_scores <- data.frame(all_scores)
    
    df.m <- reshape2::melt(all_scores)
    p <- ggplot(data = df.m, aes(x=variable, y=value)) + geom_boxplot()
    p <- p + facet_wrap( ~ variable, scales="free", strip.position = "left")
    p <- p + theme_bw() + 
      theme(axis.title.y = element_blank(), 
            axis.title.x = element_blank(), 
            axis.text.x = element_blank())
    
    return(p)
  })
  
  output$scores_stats <- renderPlot({
    suppressWarnings(compute_scores_stats())
  })
  
  output$scores_missing_data <- DT::renderDataTable({
    details <- as.matrix(compute_scores_missing_data())
    # if(is.null(detail)){
    #   detail <- as.matrix(data.frame(waiting = ""))
    # }
    DT::datatable(details,
                  rownames = TRUE,
                  options = list(scrollX = TRUE,
                                 lengthChange = FALSE,
                                 searching = FALSE))
  })

  
  #### SCORES CORRELATION MATRICE ####
  compute_scores_corr <- eventReactive(input$get_predictors_info, {
    sub_matrice <- diag(nrow = 10, ncol = 10)
    
    if(length(input$predictors) > 0){
      sub_matrice <- scores_correlation_matrice[rownames(scores_correlation_matrice) %in% input$predictors,
                                                colnames(scores_correlation_matrice) %in% input$predictors]
    }
    
    return(sub_matrice)
  })
  
  output$scores_corr <- renderD3heatmap({
    m <- compute_scores_corr()
    if(!is.null(m)){
      d3heatmap(x = m, Colv = "Rowv", 
                Rowv = NULL, colors = cor_color_breaks, 
                cexCol = .7, na.rm = F, cexRow = .7)
    }
    
  })
  
  
  #### METASCORES ####
  observeEvent(input$update_metascore, {
    old_figures <- dir(path = "www/scores_figures/", full.names = T)
    file.remove(old_figures)
    
    scores_data <- build_score_nodes(values,
                                     selected_adj_scores = input$adj_scores, 
                                     selected_raw_scores = input$raw_scores, 
                                     inc = input$update_metascore)
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
    values$selected_node <- input$current_node_id
  })
  
  
  #### SVN SCORES DETAIL ####
  output$snv_score_details <- DT::renderDataTable({
    if(values$selected_node == '') values$selected_node <- as.character(values$annotations$query)[1]
    snv_score_details <- values$scores_data$nodes[values$scores_data$nodes$group == paste0("Scores_",values$selected_node), c("id","color","label")]
    color_code <- as.character(snv_score_details$color)
    snv_score_details$id <- gsub(x = snv_score_details$id, pattern = paste0("_", values$selected_node), replacement = "")
    snv_score_details <- snv_score_details[,-2] #suppression de la colonne 'color' 
    
    DT::datatable(snv_score_details,
                  rownames = FALSE,
                  options = list(scrollX = TRUE,
                                 lengthChange = FALSE,
                                 searching = FALSE))
  }%>% DT::formatStyle(
    'label',
    target = 'row',
    backgroundColor = DT::styleEqual(as.character(snv_score_details$label), color_code)
  ))
  
  output$snv_score_details_id <- renderText({
    if(values$selected_node == '') values$selected_node <- as.character(values$annotations$query)[1]
    values$selected_node
  })
}



