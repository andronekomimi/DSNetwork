require(shiny)
require(plot3D)

#### load demo data ####
load(file = 'demo/demo.RData')


server <- function(input, output, session) {
  
  values <- reactiveValues()
  values$nodes <- nodes
  values$score_nodes <- score_nodes
  values$subnodes <- subnodes
  
  isolate(my_nodes <- rbind(values$nodes, values$score_nodes, values$subnodes))
  my_edges <- rbind(edges, score_edges, subedges,corr_edges)
  
  current_nodes <- my_nodes
  current_edges <- my_edges
  
  
  # MAJ des wigdets
  observe({ 
    updateSelectInput(session, "focus",
                      choices = c("None", candidate_SNP$RsID),
                      selected = "None"
    )
    updateSliderInput(session, inputId = "ld_range", min = min(ld_values), max = max(ld_values))
  })
  
  # LD range
  observeEvent(input$update_ld, { 
    # data.frame edges, use colomapping
    color_2_keep <- colormapping[(names(colormapping) <=  input$ld_range[2] & names(colormapping) >= input$ld_range[1])]
    edges_2_remove <- as.character(edges[!edges$color %in% color_2_keep, ]$id)
    
    current_edges <<- my_edges[!my_edges$id %in% edges_2_remove,]
      
    visNetworkProxy("network_hello") %>%
      visUpdateEdges(my_edges)
    
    visNetworkProxy("network_hello") %>%
      visRemoveEdges(id = edges_2_remove)
  })
  
  
  # Annotations
  observeEvent(input$update_annotations, {
    wanted_annotations <- selected_annotations[as.numeric(input$annotations)]
    wanted_annotations <- paste(wanted_annotations, collapse = "|")
    wanted_annotations <- gsub(x = wanted_annotations, pattern = "\\+\\+gt2", replacement = "")
    nodes_2_remove <- as.character(subnodes[!grepl(x = subnodes$id, pattern = paste(wanted_annotations, collapse = "|"), perl = T),]$id)
    
    current_nodes <<- my_nodes[!my_nodes$id %in% nodes_2_remove,]
    # print(current_edges)
    # print(current_nodes)
    # # visNetworkProxy("network_hello") %>%
    # #   visUpdateNodes(current_nodes) %>%
    # #   visUpdateEdges(current_edges)
    # 
    # visNetworkProxy("network_hello") %>%
    #   visUpdateEdges(current_edges) %>%
    #   visRemoveNodes(id = nodes_2_remove)
  })
  
  # Meta-score
  observeEvent(input$update_metascore, {
    
    meta_scores <- rep(0, times = nrow(candidate_SNP))
    b <- split(subnodes, f = subnodes$group)

    for(i in seq_along(b)){
      meta_score <- 0

      for(j in as.numeric(input$metascore)){
        classement <- which(as.character(b[[i]]$color[j]) == all_palettes[[j]])
        if(length(classement) == 0){
          classement <- nrow(candidate_SNP)
        }
        meta_score <- meta_score + classement
      }
      meta_scores[i] <- meta_score
    }
    
    print(meta_scores)

    meta_colpalette <- colfunc(n = length(unique(meta_scores)))
    names(meta_colpalette) <- sort(unique(meta_scores), decreasing = F)

    meta_values_mapped <- meta_scores
    for(i in unique(meta_scores)){
      meta_values_mapped[meta_values_mapped == i] <- meta_colpalette[names(meta_colpalette) == i]
    }
    
    values$nodes$color <- meta_values_mapped
    values$score_nodes$color <- meta_values_mapped
    
    visNetworkProxy("network_hello") %>%
      visUpdateNodes(rbind(values$nodes, values$score_nodes, values$subnodes))
  })
  
  createNetwork <- function(){
    vn <- visNetwork(my_nodes, my_edges) %>%
      # visLegend(addEdges = ledges, addNodes = lnodes, useGroups = F) %>%
      # visOptions(highlightNearest = TRUE, selectedBy = "group") %>%
      visInteraction(navigationButtons = TRUE) %>%
      visClusteringByGroup(groups = paste0("Scores_",candidate_SNP$RsID), label = "") %>%
      visPhysics(solver = "forceAtlas2Based", maxVelocity = 20,
                 forceAtlas2Based = list(gravitationalConstant = -300))
    
    for (snp in candidate_SNP$RsID){
      vn <- vn %>% visGroups(groupname = paste0("Scores_",snp), background = "#97C2FC", 
                             color = "#2B7CE9", shape = "square")
    }
    return(vn)
  }
  
  output$network_hello <- renderVisNetwork({
    createNetwork()
  })
  
  
  
  observe({ 
    if(input$focus != "None"){
      visNetworkProxy("network_hello") %>%
        visFocus(id = input$focus, scale = 1)
    } else {
      visNetworkProxy("network_hello") %>%
        visFit()
    }
  })
  
  
  #### LEGENDS ####
  drawLegendPlot <- function(){
    
    op <- par(mfrow = c(2,4),
              oma = c(0,0,0,0) + 0.1,
              mar = c(0,0,1,1) + 0.1) 
    
    colkey(side = 1, col = rev(colormapping), 
           clim = range(as.numeric(names(colormapping))), add = FALSE, clab = "LD",
           col.clab = "black", adj.clab = 0)
    
    colkey(side = 1, col = rev(cadd_colpalette), 
           clim = range(as.numeric(names(cadd_colpalette))), add = FALSE, clab = "CADD",
           col.clab = "black", adj.clab = 0)
    
    colkey(side = 1, col = rev(eigen_colpalette),
           clim = range(as.numeric(names(eigen_colpalette))), add = FALSE, clab = "EIGEN",
           col.clab = "black", adj.clab = 0)
    
    colkey(side = 1, col = rev(fathmm_colpalette),
           clim = range(as.numeric(names(fathmm_colpalette))), add = FALSE, clab = "FATHMM",
           col.clab = "black", adj.clab = 0)
    
    colkey(side = 1, col = rev(gerp_colpalette),
           clim = range(as.numeric(names(gerp_colpalette)), na.rm = T), add = FALSE, clab = "GERP",
           col.clab = "black", adj.clab = 0)
    
    colkey(side = 1, col = rev(gwava_colpalette),
           clim = range(as.numeric(names(gwava_colpalette))), add = FALSE, clab = "GWAVA",
           col.clab = "black", adj.clab = 0)
    
    colkey(side = 1, col = lin_colpalette,
           clim = range(as.numeric(names(lin_colpalette)),na.rm = T), add = FALSE, clab = "LINSIGHT",
           col.clab = "black", adj.clab = 0)
    
    par(op)
  }
  
  output$color_key <- renderPlot(drawLegendPlot())
  
}