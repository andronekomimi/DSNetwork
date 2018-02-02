library(shiny)
library(plot3D)
library(visNetwork)
options(shiny.trace = TRUE)
#### load demo data ####
### load(file = 'demo/demo.RData')
source('helper.R', local = TRUE)

server <- function(input, output, session) {
  
  #### CONF ####
  tmpDir <- 'temp/'
  app.conf <- list(TABIX = '/usr/local/bin/tabix',
                   VCF = '/Users/nekomimi/Workspace/vexor/vexor/data/1000Genomes/')
  
  values <- reactiveValues()
  values$maf_infos <- as.matrix(data.frame(waiting = ""))
  values$annotations <- as.matrix(data.frame(waiting = ""))
  values$edges <- NULL
  values$nodes <- NULL
  
  
  #### INPUT DATA MANAGEMENT ####
  transform_query <- function(string){
    string <- unlist(strsplit(x = string, split = "\n"))
    string <- gsub(x = string, pattern = " +$", replacement = "")
    print(string)
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
  
  query_control <- eventReactive(input$fetch_annotations, {
    modstring <- transform_query(input$query)
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
    
    modstring <- transform_query(input$query)
    valid_transfo <- modstring[!grepl(x = modstring, pattern = 'FAIL')]
    if(length(valid_transfo) > 0){
      #### run myvariant ####
      res <- as.data.frame(getVariants(hgvsids = valid_transfo,
                                       verbose = F, return.as = "DataFrame",
                                       fields = c("dbsnp","cadd","dbnsfp")), 
                           stringsAsFactors = TRUE)
      res$linsight <- NA
      save(res, file = 'res.rda')
      values$res <- res
      
      if(!is.null(values$res) && nrow(values$res) > 0){
        
        #### ranges ####
        my_ranges <- apply(X = values$res, MARGIN = 1, 
                           FUN = function(x) hgvsToGRange(hgvs_id = x[names(x) == "X_id"]$X_id, 
                                                          query_id = x[names(x) == "query"]$query))
        my_ranges <- do.call("c", my_ranges) 
        my_ranges <- sort(unique(my_ranges))
        save(my_ranges, file = "my_ranges.rda")
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
              res[res$query == query,"linsight"] <- value
            }
          }
        }
        if(sum(is.na(res$linsight)) == length(res$linsight)){ # no lINSIGHT RESULTS
          res$linsight <- NULL
        }
        
        values$res <- res
        save(res, file = 'res.rda')
        
        #### population freq ####
        maf_infos <- values$res[,grepl(x = colnames( values$res), pattern = "query|cadd.ref|cadd.alt|cadd.1000g|maf")] 
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
        save(annotations_infos, file = "annotations.rda")
        
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
        
        save(raw_scores, adjusted_scores, file = "scores.rda")
        values$raw_scores <- raw_scores
        values$adjusted_scores <- adjusted_scores
      }
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
  runLD <- eventReactive(input$runLD,{
    
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
    })
    
    values$ld <- ld_results
    save(ld_results, file = "ld_results.rda")
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
  
  #### BUILDING NETWORK ####
  buildNetwork <- eventReactive(input$buildNetwork,{
    snv_edges_range = switch (input$snv_edges_type,
                              "0" = input$dist_range,
                              "1" = input$ld_range
    )
    
    values$edges <- build_snv_edges(values, input$snv_edges_type, snv_edges_range)
    values$nodes <- build_snv_nodes(values)
    values$roots_score <- build_score_nodes(values, 
                                            selected_adj_scores = input$adj_scores, 
                                            selected_raw_scores = input$raw_scores)
    
    vn <- visNetwork(rbind(values$nodes,values$roots_score$nodes), 
                     rbind(values$edges,values$roots_score$edges), width = "100%", height = "1000px") %>%
      visOptions(highlightNearest = TRUE) %>%
      visPhysics(solver = "forceAtlas2Based", maxVelocity = 10,
                 forceAtlas2Based = list(gravitationalConstant = -300))
    
    return(vn)
  })
  
  output$my_network <- renderVisNetwork({
    buildNetwork()
  })
  
  #### UPDATE NETWORK ####
  
  #### LD ####
  observeEvent(input$update_ld, {
    # data.frame edges, use colomapping
    edges_2_keep <- values$edges[(values$edges$xvalue <= input$ld_range[2]) & (values$edges$xvalue >= input$ld_range[1]),]
    edges_id_2_remove <- values$edges[! values$edges$id %in% edges_2_keep$id,]$id
    
    visNetworkProxy("my_network") %>%
      visUpdateEdges(edges = edges_2_keep)
    
    visNetworkProxy("my_network") %>%
      visRemoveEdges(id = edges_id_2_remove)
  })
  
  #### DISTANCE ####
  observeEvent(input$update_dist, {
    # data.frame edges, use colomapping
    edges_2_keep <- values$edges[(values$edges$xvalue <= (1000 * input$dist_range)),]
    edges_id_2_remove <- values$edges[! values$edges$id %in% edges_2_keep$id,]$id
    
    visNetworkProxy("my_network") %>%
      visUpdateEdges(edges = edges_2_keep)
    
    visNetworkProxy("my_network") %>%
      visRemoveEdges(id = edges_id_2_remove)
  })
  
  #### FOCUS ####
  observe({
    updateSelectInput(session, "focus",
                      choices = c("None", as.character(values$nodes$id)),
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
      updateCheckboxGroupInput(session, "raw_scores",
                               choiceValues = non_null_raw_scores,
                               choiceNames = non_null_raw_scores_pretty_names
      )}
    
    if(length(non_null_adj_scores) > 0){
      updateCheckboxGroupInput(session, "adj_scores",
                               choiceValues = non_null_adj_scores, 
                               selected = non_null_adj_scores,
                               choiceNames = non_null_adj_scores_pretty_names
      )} 
    updateCheckboxGroupInput(session, "predictors",
                             choices = c(non_null_adj_scores, non_null_raw_scores)
    ) 
  })
  
  # Meta-score
  observeEvent(input$update_metascore, {
    scores_data <- build_score_nodes(values, 
                                     selected_adj_scores = input$adj_scores, 
                                     selected_raw_scores = input$raw_scores)
    
    visNetworkProxy("my_network") %>%
      visUpdateNodes(rbind(values$nodes, scores_data$nodes)) %>%
      visUpdateEdges(rbind(values$edges, scores_data$edges))
    
    
    # meta_scores <- rep(0, times = nrow(candidate_SNP))
    # b <- split(subnodes, f = subnodes$group)
    # 
    # for(i in seq_along(b)){
    #   meta_score <- 0
    #   
    #   for(j in as.numeric(input$metascore)){
    #     classement <- which(as.character(b[[i]]$color[j]) == all_palettes[[j]])
    #     if(length(classement) == 0){
    #       classement <- nrow(candidate_SNP)
    #     }
    #     meta_score <- meta_score + classement
    #   }
    #   meta_scores[i] <- meta_score
    # }
    # 
    # print(meta_scores)
    # 
    # meta_colpalette <- colfunc(n = length(unique(meta_scores)))
    # names(meta_colpalette) <- sort(unique(meta_scores), decreasing = F)
    # 
    # meta_values_mapped <- meta_scores
    # for(i in unique(meta_scores)){
    #   meta_values_mapped[meta_values_mapped == i] <- meta_colpalette[names(meta_colpalette) == i]
    # }
    # 
    # values$nodes$color <- meta_values_mapped
    # values$score_nodes$color <- meta_values_mapped
    # 
    # visNetworkProxy("network_hello") %>%
    #   visUpdateNodes(rbind(values$nodes, values$score_nodes, values$subnodes))
  })
  
  # values$nodes <- nodes
  # values$score_nodes <- score_nodes
  # values$subnodes <- subnodes
  # 
  # isolate(my_nodes <- rbind(values$nodes, values$score_nodes, values$subnodes))
  # my_edges <- rbind(edges, score_edges, subedges,corr_edges)
  # 
  # current_nodes <- my_nodes
  # current_edges <- my_edges
  # 
  # 
  # # MAJ des wigdets
  # observe({ 
  #   updateSelectInput(session, "focus",
  #                     choices = c("None", candidate_SNP$RsID),
  #                     selected = "None"
  #   )
  #   updateSliderInput(session, inputId = "ld_range", min = min(ld_values), max = max(ld_values))
  # })
  # 
  # # LD range
  # observeEvent(input$update_ld, { 
  #   # data.frame edges, use colomapping
  #   color_2_keep <- colormapping[(names(colormapping) <=  input$ld_range[2] & names(colormapping) >= input$ld_range[1])]
  #   edges_2_remove <- as.character(edges[!edges$color %in% color_2_keep, ]$id)
  #   
  #   current_edges <<- my_edges[!my_edges$id %in% edges_2_remove,]
  #   
  #   visNetworkProxy("network_hello") %>%
  #     visUpdateEdges(my_edges)
  #   
  #   visNetworkProxy("network_hello") %>%
  #     visRemoveEdges(id = edges_2_remove)
  # })
  # 
  # 
  # # Annotations
  # observeEvent(input$update_annotations, {
  #   wanted_annotations <- selected_annotations[as.numeric(input$annotations)]
  #   wanted_annotations <- paste(wanted_annotations, collapse = "|")
  #   wanted_annotations <- gsub(x = wanted_annotations, pattern = "\\+\\+gt2", 
  #                              replacement = "")
  #   nodes_2_remove <- as.character(subnodes[!grepl(x = subnodes$id, 
  #                                                  pattern = paste(wanted_annotations, 
  #                                                                  collapse = "|"), 
  #                                                  perl = T),]$id)
  #   
  #   current_nodes <<- my_nodes[!my_nodes$id %in% nodes_2_remove,]
  #   # print(current_edges)
  #   # print(current_nodes)
  #   # # visNetworkProxy("network_hello") %>%
  #   # #   visUpdateNodes(current_nodes) %>%
  #   # #   visUpdateEdges(current_edges)
  #   # 
  #   # visNetworkProxy("network_hello") %>%
  #   #   visUpdateEdges(current_edges) %>%
  #   #   visRemoveNodes(id = nodes_2_remove)
  # })
  # 
  # # Meta-score
  # observeEvent(input$update_metascore, {
  #   
  #   meta_scores <- rep(0, times = nrow(candidate_SNP))
  #   b <- split(subnodes, f = subnodes$group)
  #   
  #   for(i in seq_along(b)){
  #     meta_score <- 0
  #     
  #     for(j in as.numeric(input$metascore)){
  #       classement <- which(as.character(b[[i]]$color[j]) == all_palettes[[j]])
  #       if(length(classement) == 0){
  #         classement <- nrow(candidate_SNP)
  #       }
  #       meta_score <- meta_score + classement
  #     }
  #     meta_scores[i] <- meta_score
  #   }
  #   
  #   print(meta_scores)
  #   
  #   meta_colpalette <- colfunc(n = length(unique(meta_scores)))
  #   names(meta_colpalette) <- sort(unique(meta_scores), decreasing = F)
  #   
  #   meta_values_mapped <- meta_scores
  #   for(i in unique(meta_scores)){
  #     meta_values_mapped[meta_values_mapped == i] <- meta_colpalette[names(meta_colpalette) == i]
  #   }
  #   
  #   values$nodes$color <- meta_values_mapped
  #   values$score_nodes$color <- meta_values_mapped
  #   
  #   visNetworkProxy("network_hello") %>%
  #     visUpdateNodes(rbind(values$nodes, values$score_nodes, values$subnodes))
  # })
  # 
  # createNetwork <- function(){
  #   vn <- visNetwork(my_nodes, my_edges) %>%
  #     # visLegend(addEdges = ledges, addNodes = lnodes, useGroups = F) %>%
  #     # visOptions(highlightNearest = TRUE, selectedBy = "group") %>%
  #     visInteraction(navigationButtons = TRUE) %>%
  #     visClusteringByGroup(groups = paste0("Scores_",candidate_SNP$RsID), label = "") %>%
  #     visPhysics(solver = "forceAtlas2Based", maxVelocity = 20,
  #                forceAtlas2Based = list(gravitationalConstant = -300))
  #   
  #   for (snp in candidate_SNP$RsID){
  #     vn <- vn %>% visGroups(groupname = paste0("Scores_",snp), background = "#97C2FC", 
  #                            color = "#2B7CE9", shape = "square")
  #   }
  #   return(vn)
  # }
  # 
  # output$network_hello <- renderVisNetwork({
  #   createNetwork()
  # })
  # 
  # 
  # 
  # observe({ 
  #   if(input$focus != "None"){
  #     visNetworkProxy("network_hello") %>%
  #       visFocus(id = input$focus, scale = 1)
  #   } else {
  #     visNetworkProxy("network_hello") %>%
  #       visFit()
  #   }
  # })
  # 
  # 
  # #### LEGENDS ####
  # drawLegendPlot <- function(){
  #   
  #   op <- par(mfrow = c(2,4),
  #             oma = c(0,0,0,0) + 0.1,
  #             mar = c(0,0,1,1) + 0.1) 
  #   
  #   colkey(side = 1, col = rev(colormapping), 
  #          clim = range(as.numeric(names(colormapping))), add = FALSE, clab = "LD",
  #          col.clab = "black", adj.clab = 0)
  #   
  #   colkey(side = 1, col = rev(cadd_colpalette), 
  #          clim = range(as.numeric(names(cadd_colpalette))), add = FALSE, clab = "CADD",
  #          col.clab = "black", adj.clab = 0)
  #   
  #   colkey(side = 1, col = rev(eigen_colpalette),
  #          clim = range(as.numeric(names(eigen_colpalette))), add = FALSE, clab = "EIGEN",
  #          col.clab = "black", adj.clab = 0)
  #   
  #   colkey(side = 1, col = rev(fathmm_colpalette),
  #          clim = range(as.numeric(names(fathmm_colpalette))), add = FALSE, clab = "FATHMM",
  #          col.clab = "black", adj.clab = 0)
  #   
  #   colkey(side = 1, col = rev(gerp_colpalette),
  #          clim = range(as.numeric(names(gerp_colpalette)), na.rm = T), add = FALSE, clab = "GERP",
  #          col.clab = "black", adj.clab = 0)
  #   
  #   colkey(side = 1, col = rev(gwava_colpalette),
  #          clim = range(as.numeric(names(gwava_colpalette))), add = FALSE, clab = "GWAVA",
  #          col.clab = "black", adj.clab = 0)
  #   
  #   colkey(side = 1, col = lin_colpalette,
  #          clim = range(as.numeric(names(lin_colpalette)),na.rm = T), add = FALSE, clab = "LINSIGHT",
  #          col.clab = "black", adj.clab = 0)
  #   
  #   par(op)
  # }
  # 
  # output$color_key <- renderPlot(drawLegendPlot())
  
}