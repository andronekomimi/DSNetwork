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
  
  
  query_control <- eventReactive(input$transform_query, {
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
  observeEvent(input$transform_query, {
    values$maf_infos <- as.matrix(data.frame(waiting = ""))
    modstring <- transform_query(input$query)
    valid_transfo <- modstring[!grepl(x = modstring, pattern = 'FAIL')]
    if(length(valid_transfo) > 0){
      # run myvariant
      
      res <- as.data.frame(getVariants(hgvsids = valid_transfo,
                                       verbose = F, return.as = "DataFrame",
                                       fields = c("dbsnp","cadd","dbnsfp")))
      save(res, file = 'res.rda')
      values$res <- res
      
      if(!is.null(values$res) && nrow(values$res) > 0){
        maf_infos <- values$res[,grepl(x = colnames( values$res), pattern = "query|cadd.ref|cadd.alt|cadd.1000g|maf")] 
        #maf_infos <- maf_infos[apply(X = maf_infos, MARGIN = 1, FUN = function(x) sum(is.na(x)) != (length(x) - 1)),] # suppriner NA rows
        colnames(maf_infos) <- gsub(x = colnames(maf_infos), pattern = "cadd.1000g.(.*)", replacement = "MAF_\\1_1000G")
        row.names(maf_infos) <- NULL
        #dbsnp_infos <- data.frame(lapply(dbsnp_infos, function(x) unlist(lapply(x, paste, collapse = ","))))
        values$maf_infos <- as.matrix(maf_infos)
      }
    } 
  })
  
  output$query_res <- renderText({ 
    if(!is.null(values$res) && nrow(values$res) > 0){
      notfound_id <- values$res[!is.na(values$res$notfound),"query"]
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
  
  #### FETCH LD INFOS ####
  runLD <- eventReactive(input$runLD,{
    print('Yo')
    
    hgvsToGRange <- function(hgvs_id, query_id){
      if(is.na(hgvs_id))
        return(GRanges())
      
      chr <- gsub(x = hgvs_id, pattern = "(.*):.*", replacement = "\\1")
      start_pos <- as.numeric(gsub(x = hgvs_id, pattern = ".*:.*[a-z]\\.([0-9]+).*", replacement = "\\1"))
      hgvs_id <- unlist(strsplit(x = hgvs_id, split = ">"))
      if(length(hgvs_id) == 2){
        svn_type <- "snp"
        end_pos <- start_pos
        ref <- ""
        alt <- hgvs_id[2]
      } else {
        hgvs_id <- unlist(strsplit(x = hgvs_id, split = "ins"))
        if(length(hgvs_id) == 2){
          svn_type <- "ins"
          end_pos <- start_pos+1
          ref <- ""
          alt <- hgvs_id[2]
        } else {
          svn_type <- "del"
          end_pos <- as.numeric(gsub(x = hgvs_id, pattern = ".*_([0-9]+).*", replacement = "\\1"))
          ref <- ""
          alt <- ""
        }
      }
      
      GRanges(seqnames = chr, ranges = IRanges(start = start_pos, end = end_pos), 
              query= query_id, type = svn_type)
    }
    
    setStudyRange <- function(granges, selection){
      granges[unlist(strsplitAsListOfIntegerVectors(x = selection, sep = "_"))]
    }
    
    setStudyRegion <- function(studyrange){
      studyrange <- range(studyrange)
      return(paste0(as.character(seqnames(studyrange)), ":", start(studyrange), "-", end(studyrange)))
    }
    
    granges <- apply(X = values$res, MARGIN = 1, FUN = function(x) hgvsToGRange(hgvs_id = x[[2]], query_id = x[[1]]))
    granges <- do.call("c", granges) 
    granges <- sort(granges)
    values$granges <- granges
    
    # build ilets
    hits <- findOverlaps(query = granges, subject = granges, maxgap = 1000000)
    ilets <- c()
    
    for(i in queryHits(hits)){
      ilets <- c(ilets, paste(subjectHits(hits[queryHits(hits) == i]), collapse = "_"))
    }
    
    ilets <- unique(ilets)
    
    ld_results <- sapply(ilets, function(selection){
      current_range  <- setStudyRange(granges, selection)
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
    vn <- build_network(annotations = values$res, ld_results = values$ld)
    return(vn)
  })
  
  output$my_network <- renderVisNetwork({
    buildNetwork()
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