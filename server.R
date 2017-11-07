require(shiny)
require(plot3D)

#### load demo data ####
load(file = 'demo/demo.RData')


server <- function(input, output, session) {
  
  my_nodes <- rbind(nodes, score_nodes, subnodes)
  my_edges <- rbind(edges, score_edges, subedges,corr_edges)
  
  # MAJ des wigdets
  observe({ 
    updateSelectInput(session, "focus",
                      choices = c("None", candidate_SNP$RsID),
                      selected = "None"
    )
    updateSliderInput(session, inputId = "ld_range", min = min(ld_values), max = max(ld_values))
  })
  
  # LD range
  observe({ 
    # data.frame edges, use colomapping
    color_2_keep <- colormapping[(names(colormapping) <=  input$ld_range[2] & names(colormapping) >= input$ld_range[1])]
    edges_2_remove <- as.character(edges[!edges$color %in% color_2_keep, ]$id)
    
    visNetworkProxy("network_hello") %>%
      visUpdateEdges(my_edges)
    
    visNetworkProxy("network_hello") %>%
      visRemoveEdges(id = edges_2_remove)
  })
  
  
  # Annotations
  # observe({ 
  #   wanted_annotations <- selected_annotations[as.numeric(input$annotations)]
  #   wanted_annotations <- paste(wanted_annotations, collapse = "|")
  #   wanted_annotations <- gsub(x = wanted_annotations, pattern = "\\+\\+gt2", replacement = "")
  #   nodes_2_remove <- as.character(subnodes[!grepl(x = subnodes$id, pattern = paste(wanted_annotations, collapse = "|"), perl = T),]$id)
  #   visNetworkProxy("network_hello") %>%
  #     visUpdateNodes(my_nodes) %>%
  #     visUpdateEdges(my_edges)
  #   
  #   visNetworkProxy("network_hello") %>%
  #     visRemoveNodes(id = nodes_2_remove)
  # })
  
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
    
    gerp_colpalette <- gerp_colpalette[names(gerp_colpalette) != "."]
    lin_colpalette <- lin_colpalette[names(lin_colpalette) != "."]
    
    op <- par(mfrow = c(2,4),
              oma = c(0,0,0,0) + 0.1,
              mar = c(0,0,1,1) + 0.1) 
    
    colkey(side = 1, col = colormapping, 
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