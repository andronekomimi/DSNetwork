#DSNETWORK

require(GenomicRanges)
require(ggplot2)
require(grDevices) 
require(grid)
require(yaml)
require(c3)

#devtools::install_github("mrjoh3/c3")
#devtools::install_github("viking/r-yaml")

path_to_images <- "~/Workspace/DSNetwork/www/scores_figures/"
MAX_VAR <- 20

ld_breaks <- seq(0,1, by = 0.01)
colfunc <- colorRampPalette(c("yellow","red"))
ld_color_breaks <- colfunc(length(seq(0,1, by = 0.01)))
names(ld_color_breaks) <- sort(ld_breaks, decreasing = F)

cor_breaks <- signif(x = seq(-1,1, by = 0.01), digits = 2)
cor_colfunc <- colorRampPalette(c("blue","white","red"))
cor_color_breaks <- cor_colfunc(length(cor_breaks))
names(cor_color_breaks) <- sort(cor_breaks, decreasing = F)

# output$values <- renderPrint({
#   list(x4 = input$x4)
# })


hgvsToGRange <- function(hgvs_id, query_id){
  
  if(is.na(hgvs_id))
    return(GRanges())
  
  chr <- gsub(x = hgvs_id, pattern = "(.*):.*", replacement = "\\1")
  start_pos <- as.numeric(gsub(x = hgvs_id, pattern = ".*:.*[a-z]\\.([0-9]+).*", replacement = "\\1"))
  
  if(grepl(x = hgvs_id, pattern = ">")){ # snp
    svn_type <- "snp"
    hgvs_id <- unlist(strsplit(x = hgvs_id, split = ">"))
    end_pos <- start_pos
    return(GRanges(seqnames = chr, ranges = IRanges(start = start_pos, end = end_pos), 
                   query= query_id, type = svn_type))
  }
  
  if(grepl(x = hgvs_id, pattern = "del")){ #deletion
    svn_type <- "del"
    hgvs_id <- unlist(strsplit(x = hgvs_id, split = "_"))
    if(length(hgvs_id) == 2){
      end_pos <- as.numeric(gsub(x = hgvs_id[2], pattern = "([0-9]+)del", replacement = "\\1"))
    } else {
      end_pos <- start_pos
    }
    return(GRanges(seqnames = chr, ranges = IRanges(start = start_pos, end = end_pos), 
                   query= query_id, type = svn_type))
  }
  
  if(grepl(x = hgvs_id, pattern = "ins")){
    svn_type <- "ins"
    hgvs_id <- unlist(strsplit(x = hgvs_id, split = "ins"))
    end_pos <- start_pos+1
    return(GRanges(seqnames = chr, ranges = IRanges(start = start_pos, end = end_pos), 
                   query= query_id, type = svn_type))
  }
}

setStudyRange <- function(my_ranges, selection){
  my_ranges[unlist(strsplitAsListOfIntegerVectors(x = selection, sep = "_"))]
}

setStudyRegion <- function(studyrange){
  studyrange <- range(studyrange)
  return(paste0(as.character(seqnames(studyrange)), ":", start(studyrange), "-", end(studyrange)))
}

#### Modified from LDheatmap::LDheatmap functions #####
LDheatmapCustom <- function (gdat, genetic.distances = NULL, distances = "physical", 
                             LDmeasure = "r", title = "Pairwise LD", add.map = TRUE, add.key = TRUE, 
                             geneMapLocation = 0.15, geneMapLabelX = NULL, geneMapLabelY = NULL, 
                             SNP.name = NULL, color = NULL, newpage = TRUE, name = "ldheatmap", 
                             vp.name = NULL, pop = FALSE, flip = NULL, text = FALSE, updateProgress = NULL, results_file = "") {
  makeImageRect <- function(nrow, ncol, cols, name, byrow = TRUE) {
    xx <- (1:ncol)/ncol
    yy <- (1:nrow)/nrow
    
    if (byrow) {
      right <- rep(xx, nrow)
      top <- rep(yy, each = ncol)
    }
    else {
      right <- rep(xx, each = nrow)
      top <- rep(yy, ncol)
    }
    rectGrob(x = right, y = top, width = 1/ncol, height = 1/nrow, 
             just = c("right", "top"), gp = gpar(col = NA, fill = cols), 
             name = name)
  }
  makeImageText <- function(nrow, ncol, cols, name) {
    cols <- as.character(cols)
    cols[is.na(cols)] <- ""
    cols <- paste("   ", cols)
    xx <- (1:ncol)/ncol
    yy <- (1:nrow)/nrow
    right <- rep(xx, nrow)
    top <- rep(yy, each = ncol)
    textGrob(cols, x = right, y = top, gp = gpar(cex = 0.3), 
             just = c("right", "top"), name = name)
  }
  if (is.null(color)) {
    if (inherits(gdat, "LDheatmap")) 
      color <- gdat$color
    else color <- grey.colors(20)
  }
  LDheatmap.Legend.add <- function(color, vp = heatmapVP) {
    ImageRect <- makeImageRect(2, length(color), col = c(rep(NA, 
                                                             length(color)), color[length(color):1]), "colorKey")
    keyVP <- viewport(x = 1.1, y = -0.1, height = 0.1, width = 0.5, 
                      just = c("right", "bottom"), name = "keyVP")
    if (LDmeasure == "r") {
      ttt <- expression(paste(R^2, " Color Key"))
    }
    else {
      ttt <- "D' Color Key"
    }
    title <- textGrob(ttt, x = 0.5, y = 1.25, name = "title", 
                      gp = gpar(cex = 0.8))
    labels <- textGrob(paste(0.2 * 0:5), x = 0.2 * 0:5, y = 0.25, 
                       gp = gpar(cex = 0.6), name = "labels")
    ticks <- segmentsGrob(x0 = c(0:5) * 0.2, y0 = rep(0.4, 
                                                      6), x1 = c(0:5) * 0.2, y1 = rep(0.5, 6), name = "ticks")
    box <- linesGrob(x = c(0, 0, 1, 1, 0), y = c(0.5, 1, 
                                                 1, 0.5, 0.5), name = "box")
    key <- gTree(children = gList(ImageRect, title, labels, 
                                  ticks, box), name = "Key", vp = keyVP)
    key
  }
  if (is.null(flip)) {
    if (inherits(gdat, "LDheatmap") && !is.null(gdat$flipVP)) 
      flip <- TRUE
    else flip <- FALSE
  }
  LDheatmap.Map.add <- function(nsnps, add.map, genetic.distances, 
                                geneMapLocation = 0.15, geneMapLabelX = NULL, geneMapLabelY = NULL, 
                                distances = "physical", vp = NULL, SNP.name = NULL, ind = 0, 
                                flip = FALSE) {
    snp <- ((1:nsnps - 1) + 0.5)/nsnps
    if (add.map) {
      min.dist <- min(genetic.distances)
      max.dist <- max(genetic.distances)
      total.dist <- max.dist - min.dist
      if (flip) 
        geneMapLocation <- (-geneMapLocation)
      seq.x <- c(0.5 * geneMapLocation + 1/(nsnps * 2), 
                 1 + 0.5 * geneMapLocation - 1/(nsnps * 2))
      seq.y <- c(-0.5 * geneMapLocation + 1/(nsnps * 2), 
                 1 - 0.5 * geneMapLocation - 1/(nsnps * 2))
      diagonal <- linesGrob(seq.x, seq.y, gp = gpar(lty = 1), 
                            name = "diagonal", vp = vp)
      regionx <- seq.x[1] + ((genetic.distances - min.dist)/total.dist) * 
        (seq.x[2] - seq.x[1])
      regiony <- seq.y[1] + ((genetic.distances - min.dist)/total.dist) * 
        (seq.y[2] - seq.y[1])
      segments <- segmentsGrob(snp, snp, regionx, regiony, 
                               name = "segments", vp = vp)
      if (distances == "physical") 
        mapLabel <- paste("Physical Length:", round((total.dist/1000), 
                                                    1), "kb", sep = "")
      else mapLabel <- paste("Genetic Map Length:", round(total.dist, 
                                                          1), "cM", sep = "")
      if (!flip) {
        if (is.null(geneMapLabelY)) 
          geneMapLabelY <- 0.3
        if (is.null(geneMapLabelX)) 
          geneMapLabelX <- 0.5
      }
      else {
        if (is.null(geneMapLabelY)) 
          geneMapLabelY <- 0.8
        if (is.null(geneMapLabelX)) 
          geneMapLabelX <- 0.4
      }
      title <- textGrob(mapLabel, geneMapLabelX, geneMapLabelY, 
                        gp = gpar(cex = 0.9), just = "left", name = "title")
      geneMap <- gTree(children = gList(diagonal, segments, 
                                        title), name = "geneMap")
      if (!is.null(SNP.name) && (any(ind != 0))) {
        symbols <- pointsGrob(snp[ind], snp[ind], pch = "*", 
                              gp = gpar(cex = 1.25, bg = "blue", col = "blue"), 
                              name = "symbols", vp = vp)
        SNPnames <- textGrob(paste(" ", SNP.name), just = "left", 
                             rot = -45, regionx[ind], regiony[ind], gp = gpar(cex = 0.6, 
                                                                              col = "blue"), name = "SNPnames", vp = vp)
        if (flip) {
          lenght_SNP_name <- max(nchar(SNP.name))
          long_SNP_name <- paste(rep(8, lenght_SNP_name), 
                                 collapse = "")
          name_gap <- convertWidth(grobWidth(textGrob(long_SNP_name)), 
                                   "npc", valueOnly = TRUE)/sqrt(2)
          diagonal <- linesGrob(seq.x, seq.y, gp = gpar(lty = 1), 
                                name = "diagonal", vp = vp)
          segments <- segmentsGrob(snp, snp, regionx, 
                                   regiony, name = "segments", vp = vp)
          symbols <- NULL
          SNPnames <- textGrob(SNP.name, just = "left", 
                               rot = -45, regionx[ind] - name_gap, regiony[ind] + 
                                 name_gap, gp = gpar(cex = 0.6, col = "green"), 
                               name = "SNPnames", vp = vp)
          title <- editGrob(title, y = unit(geneMapLabelY + 
                                              name_gap, "npc"))
        }
        geneMap <- gTree(children = gList(diagonal, segments, 
                                          title, symbols, SNPnames), name = "geneMap")
      }
    }
    else if (!add.map && !is.null(SNP.name) && (any(ind != 
                                                    0))) {
      geneMap <- textGrob(paste(" ", SNP.name), just = "left", 
                          rot = -45, snp[ind], snp[ind], gp = gpar(cex = 0.6, 
                                                                   col = "blue"), name = "SNPnames")
      if (flip) 
        geneMap <- editGrob(geneMap, vp = vp)
    }
    else geneMap <- NULL
    geneMap
  }
  if (is.null(genetic.distances)) {
    if (inherits(gdat, "data.frame")) 
      genetic.distances = 1000 * (1:ncol(gdat))
    else if (inherits(gdat, "matrix")) 
      genetic.distances = 1000 * (1:length(gdat[1, ]))
    else genetic.distances = gdat$genetic.distances
  }
  if (inherits(gdat, "snp.matrix")) {
    require(chopsticks)
    if (!is.vector(genetic.distances)) {
      stop("Distance should be in the form of a vector")
    }
    o <- order(genetic.distances)
    genetic.distances <- genetic.distances[o]
    gdat <- gdat[, o]
    myLD <- ld.snp(gdat, depth = ncol(gdat))
    if (LDmeasure == "r") 
      LDmatrix <- myLD[["rsq2"]]
    else if (LDmeasure == "D") 
      LDmatrix <- myLD[["dprime"]]
    else stop("Invalid LD measurement, choose r or D'.")
    nsnp <- length(genetic.distances)
    tem <- matrix(NA, nrow = nsnp, ncol = nsnp)
    for (i in 1:(nsnp - 1)) {
      tem[i, (i + 1):nsnp] <- LDmatrix[i, 1:(nsnp - i)]
    }
    LDmatrix <- tem
    row.names(LDmatrix) <- attr(myLD, "snp.names")
  }
  else if (inherits(gdat, "data.frame")) {
    for (i in 1:ncol(gdat)) {
      if (!genetics::is.genotype(gdat[, i])) 
        stop("column ", i, " is not a genotype object\n")
    }
    gvars <- unlist(sapply(gdat, function(x) genetics::nallele(x) == 
                             2))
    genetic.distances <- genetic.distances[gvars]
    gdat <- gdat[gvars]
    if (!is.vector(genetic.distances)) {
      stop("Distance should be in the form of a vector")
    }
    o <- order(genetic.distances)
    genetic.distances <- genetic.distances[o]
    gdat <- gdat[, o]
    
    ##### MODIFIED TO ADD PRINT
    myLD <- LD.data.frame.custom(g1 = gdat, updateProgress = updateProgress)
    if (LDmeasure == "r") 
      LDmatrix <- myLD[[LDmeasure]] # R instead of R-squared
    else if (LDmeasure == "D'") 
      LDmatrix <- abs(myLD[[LDmeasure]])
    else stop("Invalid LD measurement, choose r or D'.")
  }
  else if (inherits(gdat, "LDheatmap")) {
    LDmatrix <- gdat$LDmatrix
    distances <- gdat$distances
  }
  else if (inherits(gdat, "matrix")) {
    if (nrow(gdat) != ncol(gdat)) 
      stop("The matrix of linkage disequilibrium measurements must be a square matrix")
    LDmatrix <- gdat
    LDmatrix[lower.tri(LDmatrix, diag = TRUE)] <- NA
  }
  else if (!missing(gdat)) 
    stop(paste("No method for an object of class", class(gdat)))
  else stop("Need to supply LD matrix or genotypes")
  heatmapVP <- viewport(width = unit(0.8, "snpc"), height = unit(0.8, 
                                                                 "snpc"), name = vp.name)
  flipVP <- viewport(width = unit(0.8, "snpc"), height = unit(0.8, 
                                                              "snpc"), y = 0.6, angle = -45, name = "flipVP")
  if (color[1] == "blueToRed") 
    color = rainbow(20, start = 4/6, end = 0, s = 0.7)[20:1]
  png(filename = results_file, width = 500, height = 500, bg = NA)
  if (newpage) 
    grid.newpage()
  mybreak <- 0:length(color)/length(color)
  imgLDmatrix <- LDmatrix
  byrow <- ifelse(flip, FALSE, TRUE)
  colcut <- as.character(cut(1 - imgLDmatrix, mybreak, labels = as.character(color), 
                             include.lowest = TRUE))
  if (is.numeric(color)) 
    colcut <- as.integer(colcut)
  ImageRect <- makeImageRect(dim(LDmatrix)[1], dim(LDmatrix)[2], 
                             colcut, name = "heatmap", byrow)
  ImageText <- NULL
  if (text) 
    ImageText <- makeImageText(dim(LDmatrix)[1], dim(LDmatrix)[2], 
                               round(imgLDmatrix, digits = 2), name = "heatmaptext")
  title <- textGrob(title, 0.5, 1.05, gp = gpar(cex = 1), name = "title")
  if (flip) {
    ImageRect <- editGrob(ImageRect, vp = flipVP)
    if (text) 
      ImageText <- editGrob(ImageText, vp = flipVP, rot = 45, 
                            just = "left")
  }
  heatMap <- gTree(children = gList(ImageRect, ImageText, title), 
                   name = "heatMap")
  nsnps <- ncol(LDmatrix)
  step <- 1/(nsnps - 1)
  ind <- match(SNP.name, row.names(LDmatrix), nomatch = 0)
  geneMapVP <- NULL
  if (flip) 
    geneMapVP <- flipVP
  geneMap <- LDheatmap.Map.add(nsnps, genetic.distances = genetic.distances, 
                               geneMapLocation = geneMapLocation, add.map, geneMapLabelX = geneMapLabelX,
                               geneMapLabelY = geneMapLabelY, distances = distances, 
                               vp = geneMapVP, SNP.name = SNP.name, ind = ind, flip = flip)
  if (add.key) 
    Key <- LDheatmap.Legend.add(color, vp = heatmapVP)
  else Key <- NULL
  LDheatmapGrob <- gTree(children = gList(heatMap, geneMap, 
                                          Key), vp = heatmapVP, name = name, cl = "ldheatmap")
  grid.draw(LDheatmapGrob)
  if (pop) {
    downViewport(heatmapVP$name)
    popViewport()
  }
  dev.off()
  ldheatmap <- list(LDmatrix = LDmatrix, LDheatmapGrob = LDheatmapGrob, 
                    heatmapVP = heatmapVP, flipVP = geneMapVP, genetic.distances = genetic.distances, 
                    distances = distances, color = color)
  class(ldheatmap) <- "LDheatmap"
  invisible(ldheatmap)
}

#### Modified from genetics::LD.data.frame function  ####
LD.data.frame.custom <- function(g1, updateProgress = NULL, ...){
  require(shiny)
  gvars <- sapply( g1, function(x) (is.genotype(x) && nallele(x)==2) )
  if(any(gvars==FALSE))
  {
    warning("Non-genotype variables or genotype variables ",
            "with more or less than two alleles detected. ",
            "These variables will be omitted: ",                
            paste( colnames(g1)[!gvars] , collapse=", " )
    )
    g1 <- g1[,gvars]
  }
  
  
  P <- matrix(nrow=ncol(g1),ncol=ncol(g1))
  rownames(P) <- colnames(g1)
  colnames(P) <- colnames(g1)
  
  P <- D <- Dprime <- nobs <- chisq <- p.value <- corr <- R.2 <- P
  for(i in 1:(ncol(g1)-1) )
    for(j in (i+1):ncol(g1) )
    {
      if (is.function(updateProgress)) {
        text <- paste(rownames(P)[i], 'vs',rownames(P)[j])
        updateProgress(detail = text)
      }
      
      ld <- LD( g1[,i], g1[,j] )
      
      D      [i,j] <- ld$D
      Dprime [i,j] <- ld$"D'"
      corr   [i,j] <- ld$"r"
      R.2    [i,j] <- ld$"R^2"          
      nobs   [i,j] <- ld$"n"
      chisq  [i,j] <- ld$"X^2"
      p.value[i,j] <- ld$"P-value"
    }
  
  retval <- list(
    call=match.call(),
    "D"=D,
    "D'"=Dprime,
    "r" = corr,
    "R^2" = R.2,
    "n"=nobs,
    "X^2"=chisq,
    "P-value"=p.value
  )
  
  class(retval) <- "LD.data.frame"
  
  retval
}

computeLDHeatmap <- function(region, requested_variants, results_dir, 
                             vcf_dir, population, tempDir, tabix_path, i, isShiny = TRUE){
  require(genetics)
  require(LDheatmap)
  
  requested_variants <- unique(gsub(x = requested_variants, pattern = "(.*)_(.*)", replacement = "\\1"))
  region <- gsub(x = region, pattern = 'chr', replacement = '')
  chromosome <- strsplit(x = region, split = ":")[[1]][1]
  
  command_line <- paste0("perl scripts/vcf_to_ped_converter.pl -vcf ", vcf_dir,"ALL.chr",
                         chromosome, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
                         " -sample_panel_file ",vcf_dir,"integrated_call_samples_v3.20130502.ALL.panel",
                         " -region ", region,
                         " -population ", population,
                         " -tabix ",tabix_path,
                         " -output_dir ", results_dir,
                         " -base_format letter")
  
  warning(command_line)
  
  ret <- system(command = command_line, intern = TRUE)
  
  region <- gsub(x = region, pattern = ':', replacement = '_')
  genotypes <- read.csv(paste0(results_dir, "/",region, ".ped"),
                        header = F, sep = '\t',
                        stringsAsFactors = F)
  
  found_variants <- read.table(paste0(results_dir, "/",region, ".info"),
                               stringsAsFactors = F)
  
  header <- c('family', 'individus','pater', 'mater', 'sex', 'pheno', found_variants$V1)
  colnames(genotypes) <- header
  
  variants <- c(requested_variants,found_variants$V1)[duplicated(c(requested_variants,found_variants$V1))]
  notfound = requested_variants[!requested_variants %in% variants]
  
  
  if(length(variants) > 1) {
    genotypes <- subset(x = genotypes, select = variants)
    variants_pos <- subset(found_variants, V1 %in% variants)
    df <- lapply(genotypes, FUN = function(x) gsub(x = x, pattern = " ", replacement = "/"))
    
    my_data <- makeGenotypes(df)
    print('done making genotypes')
    
    if(isShiny) {
      progress <- shiny::Progress$new()
      progress$set(message = "Computing LD...", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / nrow(combn2(x = variants))
        }
        progress$set(value = value, detail = detail)
      }
    } else {
      updateProgress = NULL
    }
    
    
    ldheatmap <- LDheatmapCustom(gdat = my_data, 
                                 genetic.distances = variants_pos$V2,
                                 updateProgress = updateProgress, SNP.name = variants,
                                 results_file = paste0(results_dir,"/ld_figures/LD_plot_", region, ".png"))
    
    ### add lost variants in LDheatmap
    lost_variants <- variants[!variants %in% colnames(ldheatmap$LDmatrix)]
    notfound <- c(notfound, lost_variants)
    
    ### add not found variants
    rounded_LD_results <- round(ldheatmap$LDmatrix, digits = 3)
    empty_test <- matrix(data = NA, nrow = nrow(rounded_LD_results), ncol = length(notfound))
    colnames(empty_test) <- notfound
    rounded_LD_results <- cbind(rounded_LD_results, empty_test)
    empty_test2 <- matrix(data = NA, ncol = ncol(rounded_LD_results), nrow = length(notfound))
    rownames(empty_test2) <- notfound
    rounded_LD_results <- rbind(rounded_LD_results, empty_test2)
    
    rounded_LD_results <- subset(rounded_LD_results, select =  requested_variants)
    rounded_LD_results <- rounded_LD_results[requested_variants,]  
    
    ### add diagonal 1
    diag(rounded_LD_results) <- 1
    lowerTriangle(rounded_LD_results) <- upperTriangle(rounded_LD_results, byrow = T)
    lowerTriangle(rounded_LD_results) <- lowerTriangle(rounded_LD_results) * -1
    
    ### replace NA by 0....
    rounded_LD_results[is.na(rounded_LD_results)] <- 0  
    
    ### save matrix in file
    ld_file = paste0(tempDir,.Platform$file.sep, 'LD_matrix',i,'.csv')
    write.table(x = rounded_LD_results, file = ld_file, quote = FALSE, row.names = F, col.names = F)
    
    print('done')
    
    # R to R squared for the graphic representation
    ldheatmap$LDmatrix <- ldheatmap$LDmatrix^2
    
    return(list( notfound = notfound,
                 filepath = ld_file,
                 data = ldheatmap,
                 SNP.name = variants))
  } else {
    return(list(notfound = requested_variants,
                filepath = '',
                data = NULL,
                SNP.name = NULL))
  }
}

build1000GenomeLDMatrixLDlink <- function(variantsID, population, tempDir, updateProgress = NULL) {
  require(dplyr)
  require(httr)
  rsids = variantsID[grepl(x = variantsID, pattern = 'rs')]
  snpids = variantsID[!grepl(x = variantsID, pattern = 'rs')]
  
  results = NULL
  variants <- paste(rsids, collapse = '%0A')
  tempID <- format(Sys.time(), "VXR%y%m%d%H%M%S")
  LDLink_url <- paste0('http://analysistools.nci.nih.gov/LDlinkRest/ldmatrix?snps=',variants,
                       '&pop=',population,
                       '&reference=',tempID,
                       '&r2_d=r2')
  req <- GET(LDLink_url)
  
  if(req$status_code == 200){
    
    ld_file = paste0(tempDir,.Platform$file.sep, 'LD_matrix.csv')
    
    system(command = paste0('wget -nv -O ',ld_file,' http://analysistools.nci.nih.gov/LDlink/tmp/r2_',tempID,'.txt'))
    results_df <-  NULL
    
    tryCatch({
      results_df <- read.csv(file = ld_file, header = TRUE, sep = '\t')
    }, error = function(err) {
      print(err)
    })
    
    if(!is.null(results_df) && nrow(results_df) > 1) {
      results_df[]
      m <- matrix(nrow = length(variantsID), ncol = length(variantsID), data = -1)
      colnames(m) <- c(rsids, snpids)
      rownames(m) <- c(rsids, snpids)
      
      for(i in seq(1:length(variantsID))) {
        m[i,i] = 1
      }
      
      results_rownames <- as.character(results_df$RS_number)
      found_variants <- names(results_df)[grepl(x = names(results_df), pattern = 'rs[0-9]+$')]
      results_df <- dplyr::select_(results_df, .dots = found_variants)
      
      
      for (rsid in found_variants) {
        
        if (is.function(updateProgress)) {
          text <- rsid
          updateProgress(detail = text)
        }
        
        res = results_df[[rsid]]
        names(res) <- results_rownames
        others = found_variants[found_variants != rsid]
        for(other in others) {
          tryCatch({
            m[rsid,other] = m[other,rsid] = as.numeric(unique(res[names(res) == other]))
          }, error = function(err) {
            print(err)
          })
        }
      }
      
      df = data.frame(m, stringsAsFactors = FALSE)
      df[is.na(df)] <- 0
      df <- df[, colSums(df != "") != 0]
      
      m[m == -1] <- 0
      ld_file = paste0(tempDir,.Platform$file.sep, 'LD_matrix_',tempID,'.csv')
      write.table(x = m, file = ld_file, quote = FALSE, row.names = F, col.names = F)
      
      notfound = variantsID[! variantsID %in% found_variants]
      
      return(list(summary = paste0('No LD data for the following variants : ', paste(notfound, collapse = ', ')),
                  filepath = ld_file,
                  df = df))
    }
  }
  
  return(NULL)
}

#### Network Building ####
build_snv_edges.old <- function(session_values, edges_type, edges_range, network_type){
  my_res <- session_values$my_res
  
  if(nrow(my_res) < 2){
    return(NULL)
  }
  
  edges <- NULL
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
  
  if(edges_type == "0") { #distance
    print("distance")
    hits <- findOverlaps(query = my_ranges, subject = my_ranges, maxgap = 1000 * as.numeric(edges_range)) # range in kb
    ilets <- c()
    
    for(i in unique(queryHits(hits))){
      ilets <- c(ilets, paste(subjectHits(hits[queryHits(hits) == i]), collapse = "_"))
    }
    
    ilets <- strsplitAsListOfIntegerVectors(x = unique(ilets), sep = "_")
    
    if(sum(lengths(ilets) > 1) > 0){ # presence de cluster
      for(i in seq_along(ilets)){
        ilet <- ilets[[i]]
        if(length(ilet) > 1) {
          x <- my_ranges[ilet]
          dist_comb <- combn(x = x$query, m = 2)
          if(is.matrix(dist_comb)){
            snv_dist <- apply(dist_comb, 2, function(y) width(range(x[x$query==y[1]], x[x$query==y[2]])))
            new_row <- data.frame(id = paste0("edge_dist_", i, "_", 1:ncol(dist_comb)),
                                  from = dist_comb[1,], 
                                  to = dist_comb[2,],
                                  width = 1,
                                  dashes = FALSE,
                                  xvalue = snv_dist,
                                  type = "snv_dist_edges",
                                  title = paste0(snv_dist/1000, "Kbp"),
                                  color = "black", stringsAsFactors = FALSE)
          } else {
            snv_dist <-  width(range(x[x$query==dist_comb[1]], x[x$query==dist_comb[2]])) 
            new_row <- data.frame(id = paste0("edge_dist_", i, "_", seq_along(dist_comb)),
                                  from = dist_comb[1], 
                                  to = dist_comb[2],
                                  width = 1,
                                  dashes = FALSE,
                                  xvalue = snv_dist,
                                  type = "snv_dist_edges",
                                  title =  paste0(snv_dist/1000, "Kbp"),
                                  color = "black", stringsAsFactors = FALSE)
          }
          
          if(is.null(edges)){
            edges <- new_row
          } else {
            edges <- rbind(edges, new_row, stringsAsFactors = FALSE)
          }
        }
      }
    }
  }
  
  if(edges_type == "1"){ #linkage
    
    ld_results <- session_values$ld
    nodes_edges <- c()
    
    if(!is.null(ld_results)){
      
      for(i in 1:ncol(ld_results)){
        if("LDmatrix" %in% names(ld_results[,i]$data)){
          ld <- ld_results[,i]$data$LDmatrix
          if(!is.null(ld)){
            ld_comb <- combn(x = colnames(ld), m = 2)
            
            if(is.matrix(ld_comb)){
              ld_values <- rep(x = NA, times = ncol(ld_comb))
              for (j in 1:ncol(ld_comb)) { 
                out_j = ld_comb[,j] ; 
                ld_values[j] <- ld[out_j[1], colnames(ld) == out_j[2]]
              }
              names(ld_values) <- apply(ld_comb, 2, function(x) paste(x, collapse = "_"))
            } else {
              ld_values <- ld[ld_comb[1], colnames(ld) == ld_comb[2]] 
              names(ld_values) <- paste0(ld_comb[1], "_", ld_comb[2])
            }
            
            nodes_edges <- c(nodes_edges, ld_values) 
          }
        }
      }
      
      if(length(nodes_edges) > 0){
        # colfunc<-colorRampPalette(c("red","yellow","springgreen"))
        # colormapping <- colfunc(n = length(unique(nodes_edges)))
        # names(colormapping) <- sort(unique(nodes_edges), decreasing = T)
        
        nodes_edges <- round(x = nodes_edges, digits = 2)
        ld_values_mapped <- nodes_edges
        
        for(i in unique(nodes_edges)){
          ld_values_mapped[ld_values_mapped == i] <- ld_color_breaks[names(ld_color_breaks) == i]
        }
        
        snp_comb <- unlist(strsplit(x = names(ld_values_mapped), split = "_"))
        edges <- data.frame(id = paste0("edge_ld_", 1:length(ld_values_mapped)),
                            from = snp_comb[(1:length(snp_comb) %% 2 != 0)], 
                            to = snp_comb[(1:length(snp_comb) %% 2 == 0)],
                            width = 1,
                            dashes = FALSE,
                            xvalue = nodes_edges,
                            type = "snv_ld_edges",
                            title = nodes_edges,
                            color = ld_values_mapped, stringsAsFactors = FALSE)
        
        # ajout artificiel du LD pour les variant polyalleleic
        variants_2_update <- my_res$query[grepl(x = my_res$query, pattern = "rs.*_.*")]
        corresponding_variants <- gsub(x = variants_2_update, pattern = "(.*)_(.*)", replacement = "\\1")
        
        new_rows <- NULL
        
        for(i in seq_along(corresponding_variants)){
          v <- corresponding_variants[i]
          corresponding_rows <- edges[grepl(x = rownames(edges), pattern = v),]
          if(nrow(corresponding_rows) > 0) {
            corresponding_rows[(corresponding_rows$from == v | corresponding_rows$to == v),]$id <- paste0(corresponding_rows[(corresponding_rows$from == v | corresponding_rows$to == v),]$id, "_",i)
            corresponding_rows[corresponding_rows$from == v,]$from <- corresponding_rows[corresponding_rows$to == v,]$to <- variants_2_update[i]
            
            if(is.null(new_rows)){
              new_rows <- corresponding_rows
            } else {
              new_rows <- rbind(new_rows, corresponding_rows)
            }
          }
        }
        
        edges <- rbind(edges, new_rows)
      }
    } else {
      snp_comb <- t(combn(x = my_ranges$query, m = 2))
      edges <- data.frame(id = paste0("edge_ld_",1:nrow(snp_comb)),
                          from = snp_comb[,1], 
                          to = snp_comb[,2],
                          width = 1,
                          dashes = FALSE,
                          xvalue = 0,
                          type = "snv_ld_edges",
                          title = "NA",
                          color = "#DDDDDD", stringsAsFactors = FALSE)
    }
  }
  
  return(edges)
} 

build_snv_edges <- function(session_values, edges_type, edges_range, network_type){
  my_res <- session_values$my_res
  
  if(nrow(my_res) < 2){
    return(NULL)
  }
  
  edges <- NULL
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
  
  if(edges_type == "0") { #distance
    print("distance")
    hits <- findOverlaps(query = my_ranges, subject = my_ranges, maxgap = 1000 * as.numeric(edges_range)) # range in kb
    ilets <- c()
    
    for(i in unique(queryHits(hits))){
      ilets <- c(ilets, paste(subjectHits(hits[queryHits(hits) == i]), collapse = "_"))
    }
    
    ilets <- strsplitAsListOfIntegerVectors(x = unique(ilets), sep = "_")
    
    if(sum(lengths(ilets) > 1) > 0){ # presence de cluster
      for(i in seq_along(ilets)){
        ilet <- ilets[[i]]
        if(length(ilet) > 1) {
          x <- my_ranges[ilet]
          dist_comb <- combn(x = x$query, m = 2)
          if(is.matrix(dist_comb)){
            snv_dist <- apply(dist_comb, 2, function(y) width(range(x[x$query==y[1]], x[x$query==y[2]])))
            new_row <- data.frame(id = paste0("edge_dist_", i, "_", 1:ncol(dist_comb)),
                                  from = dist_comb[1,], 
                                  to = dist_comb[2,],
                                  width = 1,
                                  dashes = FALSE,
                                  xvalue = snv_dist,
                                  type = "snv_dist_edges",
                                  title = paste0(snv_dist/1000, "Kbp"),
                                  color = "black", stringsAsFactors = FALSE)
          } else {
            snv_dist <-  width(range(x[x$query==dist_comb[1]], x[x$query==dist_comb[2]])) 
            new_row <- data.frame(id = paste0("edge_dist_", i, "_", seq_along(dist_comb)),
                                  from = dist_comb[1], 
                                  to = dist_comb[2],
                                  width = 1,
                                  dashes = FALSE,
                                  xvalue = snv_dist,
                                  type = "snv_dist_edges",
                                  title =  paste0(snv_dist/1000, "Kbp"),
                                  color = "black", stringsAsFactors = FALSE)
          }
          
          if(is.null(edges)){
            edges <- new_row
          } else {
            edges <- rbind(edges, new_row, stringsAsFactors = FALSE)
          }
        }
      }
    }
  }
  
  if(edges_type == "1"){ #linkage
    
    ld_results <- session_values$ld
    nodes_edges <- c()
    
    snp_comb <- t(combn(x = my_ranges$query, m = 2))
    edges <- data.frame(id = paste0("edge_ld_",1:nrow(snp_comb)),
                        from = snp_comb[,1], 
                        to = snp_comb[,2],
                        width = 2,
                        dashes = FALSE,
                        xvalue = 0,
                        type = "snv_ld_edges",
                        title = "NA",
                        color = "#DDDDDD", stringsAsFactors = FALSE)
    
    if(!is.null(ld_results)){
      
      for(i in 1:ncol(ld_results)){
        if("LDmatrix" %in% names(ld_results[,i]$data)){
          ld <- ld_results[,i]$data$LDmatrix
          not_found <- ld_results[,i]$notfound
          if(!is.null(ld)){
            ld[lower.tri(ld)] = t(ld)[lower.tri(ld)]
            
            for(j in 1:nrow(edges)){
              e <- edges[j,]
              e_from <- unique(gsub(x = e$from, pattern = "(.*)_(.*)", replacement = "\\1"))
              e_to <- unique(gsub(x = e$to, pattern = "(.*)_(.*)", replacement = "\\1"))
              
              if(e_from == e_to){
                e$title <- "NA"
                e$xvalue <- 0
              } else {
                if(e_from %in% colnames(ld) && e_to %in% colnames(ld)){
                  e$title <- e$xvalue <- round(x = ld[e_from, e_to], digits = 2)
                  edges[j,] <- e
                }
              }
            }
            
            # color
            ld_values_mapped <- edges$title
            
            for(x in unique(edges$title)){
              if(x == "NA"){
                ld_values_mapped[ld_values_mapped == x] <- "#DDDDDD"
              } else {
                ld_values_mapped[ld_values_mapped == x] <- ld_color_breaks[names(ld_color_breaks) == x]
              }
            }
            
            edges$color <- ld_values_mapped
          }
        }
      }
    }
  }
  
  return(edges)
} 



build_snv_nodes <- function(session_values, network_type){
  
  nodes <- NULL
  my_res <- session_values$my_res
  
  if(nrow(my_res) < 1){
    return(NULL)
  }
  
  node_names <- unique(my_res$query)
  
  nodes <- data.frame(id = node_names, 
                      color = "blue",
                      label = node_names,
                      shape = "circularImage",
                      image = paste0("scores_figures/pie_scores_",node_names,".png"),
                      font.size = 30,
                      size = 50,
                      fixed = F, x = NA, y = NA,
                      group = "Variants"
  )
  
  return(nodes)
}

build_score_nodes <- function(session_values, selected_adj_scores, selected_raw_scores, network_type, inc = NULL){
  
  load('data/scores_correlation_matrice.rda')
  colfunc <- colorRampPalette(c("red","yellow","springgreen"))
  
  my_res <- session_values$my_res
  
  switch(network_type,
         non_syn = {
           selected_scores <- selected_adj_scores
         },
         regul = {
           selected_scores <- selected_raw_scores
         }
  )
  
  nodes <- as.character(my_res$query)
  nodes_data <- data.frame(nodes = nodes)
  a_scores <- my_res[,selected_scores]
  
  # print(selected_scores)
  # print(nodes_data)
  # print(a_scores)
  nodes_data <- cbind(nodes_data, as.data.frame(a_scores))
  
  if(length(selected_scores) == 1){
    colnames(nodes_data) <- c("nodes", selected_scores)
  } 
  
  ### supprimer les lignes vides
  nodes_data <- nodes_data[(apply(X = nodes_data,
                                  MARGIN = 1, 
                                  FUN = function(x) sum(is.na(x))) < ncol(nodes_data) - 1), ]
  nodes_data <- data.frame(nodes_data, stringsAsFactors = FALSE)
  
  save(nodes_data, file = "objects/nodes_data.rda") # keep this one #
  
  if(ncol(nodes_data) >= 1){
    score_nodes <- NULL
    score_edges <- NULL
    
    scores_values_mapped <- list()
    
    # créer les echelles coloriques pour chaque scores
    for(selected_score in selected_scores){
      d <-  nodes_data[[selected_score]]
      colpalette <- colfunc(n = length(unique(d)))
      names(colpalette) <- sort(unique(d), decreasing = T, na.last = T)
      d[is.na(d)] <- names(colpalette)[is.na(names(colpalette))] <- "NA"
      
      values_mapped <- d
      
      save(values_mapped, file = "objects/values_mapped.rda")
      save(colpalette, file = "objects/colpalette.rda")
      
      for(i in unique(d)){
        if(i == "NA"){
          values_mapped[values_mapped == i] <- "#CCCCCC"
        } else {
          values_mapped[values_mapped == i] <- colpalette[names(colpalette) == i]
        }
        
      }
      scores_values_mapped[[selected_score]] <- values_mapped
    }
    
    scores_values_mapped <- data.frame(scores_values_mapped, 
                                       row.names = nodes_data$nodes, 
                                       stringsAsFactors = FALSE)
    
    
    # les scores à proprement parler
    for(n in nodes_data$nodes){
      
      scores_values <- as.numeric(nodes_data[nodes_data$nodes == n,-1])
      new_n_rows <- data.frame(id = paste0(selected_scores, "_",n),
                               color = as.character(scores_values_mapped[row.names(scores_values_mapped) == n,]),
                               label = as.character(scores_values),
                               shape = "square",
                               image = NA,
                               font.size = 10,
                               size = 30,
                               group = paste0("Scores_",n))
      
      if(is.null(score_nodes)){
        score_nodes <- new_n_rows
      } else {
        score_nodes <- rbind(score_nodes, new_n_rows)
      }
      
      #save(score_nodes, file = "objects/score_nodes.rda")
      
      
      ## creation des figures
      new_n_rows$h <- 1
      new_n_rows$id <- gsub(x = new_n_rows$id, pattern = paste0("_",n), replacement = "")
      new_n_rows$id <- as.character(new_n_rows$id)
      new_n_rows$label <- as.character(new_n_rows$label)
      new_n_rows$color <- as.character(new_n_rows$color)
      custom_colors <- new_n_rows$color
      names(custom_colors) <- new_n_rows$label
      
      ## pies
      png(paste0(path_to_images,"pie_scores_",n,inc,".png"), width = 2000, height = 2000,
          units = "px")
      par(lwd = 0.001)
      pie(x = new_n_rows$h, col = new_n_rows$color, labels = "")
      dev.off()
      par(lwd = 1)
      
    }
    
  }
  
  
  return(list(nodes = score_nodes))
}

build_snv_scores_detail_node <- function(id){
  
  node <- data.frame(id = "snv_scores_detail", 
                     color = "black",
                     label = paste0("Score details for : ",id),
                     shape = "image",
                     image = paste0("scores_figures/bar_scores_",id,".png"),
                     font.size = 30,
                     size = 175,
                     fixed = F, x = NA, y = NA,
                     group = "Detail"
  )
  
  return(node)
}

mean_score <- function(x){
  suppressWarnings(mean(x = as.numeric(unlist(strsplit(x = x, split= ","))), na.rm= T))
}

extract_score_and_convert <- function(my_res_infos, score_name){
  if(sum(colnames(my_res_infos) == score_name) == 1){
    my_data <- my_res_infos[,score_name]
    x <- sapply(X = my_data, FUN = is.numeric)
    my_data[x] <- sapply(X = my_data[x], FUN = function(y) mean(x = y, na.rm = T))
    my_data[!x] <- NA
    return(unlist(my_data))
  } else {
    return(NULL)
  }
}

extract_score_and_convert.old <- function(my_res_infos, score_name, sub_score_name = NULL){
  if(sum(colnames(my_res_infos) == score_name) == 1){
    x <- as.character(my_res_infos[,score_name])
    x <- sapply(x, mean_score)
    return(x)
  } else {
    if(!is.null(sub_score_name)){
      if(sum(colnames(my_res_infos) == sub_score_name) == 1){
        x <- as.character(my_res_infos[,sub_score_name])
        x <- sapply(x, mean_score)
        return(x)
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
    return(NULL)
  }
}

#### Basic Ranking ####
basic_ranking <- function(inc = NULL){
  colfunc <- colorRampPalette(c("springgreen","yellow","red"))
  load("objects/nodes_data.rda")
  nodes <- as.character(nodes_data$nodes)
  
  if(length(nodes) < 2){
    return(NULL)
  }
  
  # compute RANK for each classifier 
  if(ncol(nodes_data) > 2){ # nodes + at least 2 scores
    ranked_nodes_data_na_last <- apply(X = nodes_data[,-1], MARGIN = 2, 
                                       FUN = function(x) data.table::frankv(x, na.last = T, order = -1))
    # replace NA by mean and re rank
    nodes_data_na_mean <- apply(X = nodes_data[,-1], MARGIN = 2, 
                                FUN = function(x) {
                                  x[is.na(x)] <- mean(x = x, na.rm =T)
                                  return(x)
                                }) 
    
    ranked_nodes_data_na_mean <- apply(X = nodes_data_na_mean, MARGIN = 2, 
                                       FUN = function(x) data.table::frankv(x, na.last = F, order = -1))
    rownames(ranked_nodes_data_na_last) <- rownames(ranked_nodes_data_na_mean) <- nodes
    
    # compute MEAN RANK accross all the classifier 
    mean.rank_na_last = rowMeans(ranked_nodes_data_na_last)
    mean.rank_na_mean = rowMeans(ranked_nodes_data_na_mean)
    
  } else { # only on predictor
    x <- nodes_data[,-1]
    ranked_nodes_data_na_last <- data.table::frankv(x, na.last = T, order = -1)
    # replace NA by mean and re rank
    nodes_data_na_mean <- mean(x = x, na.rm =T)
    x[is.na(x)] <- nodes_data_na_mean
    ranked_nodes_data_na_mean <- data.table::frankv(x, na.last = T, order = -1)
    names(ranked_nodes_data_na_last) <- names(ranked_nodes_data_na_mean) <- nodes
    
    # compute MEAN RANK accross all the classifier 
    mean.rank_na_last = ranked_nodes_data_na_last
    mean.rank_na_mean = ranked_nodes_data_na_mean
  }
  
  
  
  # compute FINAL RANK
  rank_na_last = as.numeric(rank(mean.rank_na_last, ties.method = "average"))
  rank_na_mean = as.numeric(rank(mean.rank_na_mean, ties.method = "average"))
  
  colpalette_na_last <- colfunc(n = length(unique(rank_na_last)))
  colpalette_na_mean <- colfunc(n = length(unique(rank_na_mean)))
  
  names(colpalette_na_last) <- sort(unique(rank_na_last), decreasing = T) # plus petit ranking = meilleur = plus rouge
  names(colpalette_na_mean) <- sort(unique(rank_na_mean), decreasing = T)
  
  
  values_mapped_na_last <- as.character(rank_na_last)
  for(i in unique(rank_na_last)){
    values_mapped_na_last[values_mapped_na_last == i] <- colpalette_na_last[names(colpalette_na_last) == i]
  }
  
  values_mapped_na_mean <- as.character(rank_na_mean)
  for(i in unique(rank_na_mean)){
    values_mapped_na_mean[values_mapped_na_mean == i] <- colpalette_na_mean[names(colpalette_na_mean) == i]
  }
  
  df = data.frame(id = nodes,
                  mean.rank_na_last,
                  final_rank_na_last = rank_na_last,
                  col_na_last = as.character(values_mapped_na_last),
                  mean.rank_na_mean,
                  final_rank_na_mean = rank_na_mean,
                  col_na_mean = as.character(values_mapped_na_mean),
                  stringsAsFactors = F,
                  row.names = NULL)
  
  apply(X = df, MARGIN = 1, FUN = function(n){
    ## pies
    png(paste0(path_to_images,"pie_rank_na_last_",n[names(n) == "id"],inc,".png"),
        width = 2000, height = 2000,
        units = "px")
    par(lwd = 0.001)
    pie(x = 1, col = n[names(n) == "col_na_last"], labels = n[names(n) == "final_rank_na_last"])
    dev.off()
    par(lwd = 1)
    
    png(paste0(path_to_images,"pie_rank_na_mean_",n[names(n) == "id"],inc,".png"),
        width = 2000, height = 2000,
        units = "px")
    par(lwd = 0.001)
    pie(x = 1, col = n[names(n) == "col_na_mean"], labels = n[names(n) == "final_rank_na_mean"])
    dev.off()
    par(lwd = 1)
  })
  
}

#### Absolute metascores ####
compute_absolute_metascore <- function(session_values){
  absolute_metascore <- list(
    eigen = list(min = 1,
                 max = 99),
    eigen_pc = list(min = 1,
                    max = 99),
    linsight = list(min = 0.03,
                    max = 1),
    iwscoring_known = list(min = -5,
                     max = 6),
    iwscoring_novel = list(min = -5,
                           max = 6),
    bayesdel = list(min = -1.4,
                    max = 0.9)
  )
  
  colfunc <- colorRampPalette(c("springgreen","yellow","red"))
  colpalette_lenght <- 10
  colpalette <- colfunc(n = colpalette_lenght)
  
  absolute_metascore_colpalette <- list(
    eigen = colpalette,
    eigen_pc = colpalette,
    linsight = colpalette,
    iwscoring_known = colpalette,
    iwscoring_novel = colpalette,
    bayesdel = colpalette)
  
  for(n in names(absolute_metascore)){
    neg_direction <- absolute_metascore[[n]]$min > absolute_metascore[[n]]$max
    names(absolute_metascore_colpalette[[n]]) <- sort(round(x = seq(from = absolute_metascore[[n]]$min, 
                                                                    to = absolute_metascore[[n]]$max, 
                                                                    length.out = colpalette_lenght), digits = 2), 
                                                      decreasing = neg_direction) # plus petit ranking = meilleur = plus rouge
  }
  
  where <- function(x, n){
    if(is.na(x))
      return("#CCCCCC")
    neg_direction <- absolute_metascore[[n]]$min > absolute_metascore[[n]]$max
    xrange <- sort(round(x = seq(from = absolute_metascore[[n]]$min, 
                                 to = absolute_metascore[[n]]$max, 
                                 length.out = colpalette_lenght), digits = 2), 
                   decreasing = neg_direction)
    if(!neg_direction){
      if(sum(x >= xrange))
        return(absolute_metascore_colpalette[[n]][max(which(x >= xrange))])
      return(absolute_metascore_colpalette[[n]][1])
    } else {
      if(sum(x <= xrange))
        return(absolute_metascore_colpalette[[n]][max(which(x <= xrange))])
      return(absolute_metascore_colpalette[[n]][1])
    }
  }
  
  res_meta <- session_values$res[, c("query", names(absolute_metascore))]
  for(n in names(absolute_metascore)){
    res_meta[[n]] <- sapply(X = res_meta[[n]], FUN = function(x1) where(x = x1, n = n))
  }
  
  apply(X = res_meta, MARGIN = 1, FUN = function(n){
    ## pies
    for(i in names(absolute_metascore)){
      png(paste0(path_to_images,"pie_",i,"_",n[names(n) == "query"],".png"),
          width = 2000, height = 2000,
          units = "px")
      par(lwd = 0.001)
      pie(x = 1, col = n[names(n) == i], labels = "")
      dev.off()
      par(lwd = 1)
    }
    
    png(paste0(path_to_images,"pie_all_",n[names(n) == "query"],".png"), width = 2000, height = 2000,
        units = "px")
    par(lwd = 0.001)
    pie(x = rep(1, times = length(absolute_metascore)), col = n[names(n) != "query"], labels = "")
    dev.off()
    par(lwd = 1)
  })
  
  
  
}

#### Friedman Test : compute and plot ####
## This function has been copied from mlr/R/plotCritDifferences.R and modified to fit with our purpose
plotCritDifferences <- function(nodes_data, baseline){
  ranked_nodes_data <- apply(X = nodes_data[,-1], MARGIN = 2, FUN = function(x) data.table::frankv(x, na.last = T, order = -1))
  rownames(ranked_nodes_data) <- nodes_data$nodes
  n.tasks <- ncol(ranked_nodes_data)
  n.learners =  nrow(ranked_nodes_data)
  
  p.value = 0.05
  test = "bd"
  mean.rank = rowMeans(ranked_nodes_data)
  
  df = data.frame(mean.rank,
                  learner.id = names(mean.rank),
                  rank = rank(mean.rank, ties.method = "average"))
  
  
  # Orientation of descriptive lines yend(=y-value of horizontal line)
  right = df$rank > median(df$rank)
  # Better learners are ranked ascending
  df$yend[!right] = rank(df$rank[!right], ties.method = "first") - 0.5
  # Worse learners ranked descending
  df$yend[right] = rank(-df$rank[right], ties.method = "first") - 0.5
  # Better half of learner have lines to left / others right.
  df$xend = ifelse(!right, 0L, max(df$rank) + 1L)
  # Save orientation, can be used for vjust of text later on
  df$right = as.numeric(right)
  df$short.name = names(mean.rank)
  
  # Perform nemenyi test
  # require("PMCMR")
  f.test = friedman.test(t(ranked_nodes_data))
  f.rejnull = f.test$p.value < p.value
  if (!is.na(f.test$p.value)) {
    f.rejnull = f.test$p.value < p.value
    if (!f.rejnull)
      warning("Cannot reject null hypothesis of overall Friedman test,
              returning overall Friedman test.")
  } else {
    f.rejnull = FALSE
    warning("P-value not computable. Learner performances might be exactly equal.")
  }
  
  q.nemenyi = qtukey(1 - p.value, n.learners, 1e+06) / sqrt(2L)
  cd.nemenyi = q.nemenyi * sqrt(n.learners * (n.learners + 1L) / (6L * n.tasks))
  q.bd = qtukey(1L - (p.value / (n.learners - 1L)), 2L, 1e+06) / sqrt(2L)
  cd.bd = q.bd * sqrt(n.learners * (n.learners + 1L) / (6L * n.tasks))
  
  if (f.rejnull) {
    nem.test = PMCMR::posthoc.friedman.nemenyi.test(t(ranked_nodes_data))
    nem.test$crit.difference = list("nemenyi" = cd.nemenyi, "bd" = cd.bd)
    nem.test$f.rejnull = f.rejnull
    #return(nem.test)
  } else {
    f.test$f.rejnull = f.rejnull
    f.test$crit.difference = list("nemenyi" = cd.nemenyi, "bd" = cd.bd)
    #return(f.test)
  }
  
  # Store Info for plotting the cricital differences
  cd.info = list("test" = test,
                 "cd" = nem.test$crit.difference[[test]],
                 "x" = df$mean.rank[df$learner.id == baseline],
                 "y" = 0.1)
  
  # Plot the critical difference bars
  cd.x = df$mean.rank[df$learner.id == baseline]
  cd.y = cd.info$y
  cd = cd.info$cd
  
  # Plot descritptive lines and learner names
  p = ggplot(df)
  # Point at mean rank
  p = p + geom_point(aes_string("mean.rank", 0, colour = "learner.id"), size = 3)
  # Horizontal descriptive bar
  p = p + geom_segment(aes_string("mean.rank", 0, xend = "mean.rank", yend = "yend",
                                  color = "learner.id"), size = 1)
  # Vertical descriptive bar
  p = p + geom_segment(aes_string("mean.rank", "yend", xend = "xend",
                                  yend = "yend", color = "learner.id"), size = 1)
  # Plot Learner name
  p = p + geom_text(aes_string("xend", "yend", label = "learner.id", color = "learner.id",
                               hjust = "right"), vjust = -1)
  p = p + xlab("Average Rank")
  # Change appearance
  p = p + scale_x_continuous(breaks = c(0:max(df$xend)))
  p = p + theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "none",
                panel.background = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(size = 1),
                axis.line.y = element_blank(),
                panel.grid.major = element_blank(),
                plot.background = element_blank())
  
  # Add horizontal bar arround baseline
  p = p + annotate("segment", x = cd.x + cd, xend = cd.x - cd, y = cd.y, yend = cd.y,
                   alpha = 0.5, color = "darkgrey", size = 2)
  # Add intervall limiting bar's
  p = p + annotate("segment", x = cd.x + cd, xend = cd.x + cd, y = cd.y - 0.05,
                   yend = cd.y + 0.05, color = "darkgrey", size = 1)
  p = p + annotate("segment", x = cd.x - cd, xend = cd.x - cd, y = cd.y - 0.05,
                   yend = cd.y + 0.05, color = "darkgrey", size = 1)
  # Add point at learner
  p = p + annotate("point", x = cd.x, y = cd.y, alpha = 0.5)
  # Add critical difference text
  p = p + annotate("text", label = paste("Critical Difference =", round(cd, 2), sep = " "),
                   x = cd.x, y = cd.y + 0.05)
  return(p)
}

#### robust system call ####
## This function has been copied from 
## http://stackoverflow.com/questions/7014081/capture-both-exit-status-and-output-from-a-system-call-in-r
robust.system <- function (cmd, stdoutFile = NULL) {
  stderrFile = tempfile(pattern="R_robust.system_stderr", fileext=as.character(Sys.getpid()))
  if (is.null(stdoutFile)) stdoutFile = tempfile(pattern="R_robust.system_stdout", fileext=as.character(Sys.getpid()))
  
  retval = list()
  retval$exitStatus = system(paste0(cmd, " 2> ", shQuote(stderrFile), " > ", shQuote(stdoutFile)))
  #if (is.null(stdoutFile)) retval$stdout = readLines(stdoutFile)
  retval$stderr = readLines(stderrFile)
  
  if (is.null(stdoutFile)) 
    unlink(c(stdoutFile, stderrFile))
  else
    unlink(c(stderrFile))
  
  return(retval)
}

C3PieChartOutput <- function(outputId, width = '100%', height = '200px'){
  htmlwidgets::shinyWidgetOutput(outputId = outputId, name = 'c3', width = width, height = height, package = 'c3')
}

#' @rdname C3PieChart-shiny
#' @export
renderC3PieChart <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, C3PieChartOutput, env, quoted = TRUE)
}


createVCF <- function(session_values, filename){
  res <- session_values$res
  vcf_desc_0 <- "##fileformat=VCFv4.1"
  vcf_desc_1 <- paste0("##fileDate=",gsub(x = Sys.Date(), pattern = "-",replacement = ""))
  vcf_desc_2 <- '##INFO=<ID=MAF,Number=A,Type=Float,Description="MAF">'
  vcf_header <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  vcf_content <- paste(res$dbsnp.chrom, res$dbsnp.hg19.start, res$query, res$dbsnp.ref, res$dbsnp.alt, 100, "PASS", paste0("MAF:",res$dbsnp.gmaf), sep = "\t")
  
  write(x = vcf_desc_0, file = filename)
  write(x = vcf_desc_1, file = filename, append = T)
  write(x = vcf_desc_2, file = filename, append = T)
  write(x = vcf_header, file = filename, append = T)
  write(x = vcf_content, file = filename, append = T)
}


extractBayesDel <- function(path_to_victor, filename){
  cmd <- paste(path_to_victor,filename,"--ann=BayesDel_nsfp33a_noAF -x=3 --min=-1.5 --step=0.01 --indel=max --padding=1 -o", paste0(filename,".gz"))
  system(command = cmd, intern = F)
  if(file.exists(paste0(filename,".gz"))){
    vcf <- vcfR::read.vcfR(paste0(filename,".gz"), verbose = F)
    return(as.numeric(vcf@gt))
  } else {
    return(NULL)
  }
}

extract_LINSIGHT_range <- function(){
  min_range <- +Inf
  max_range <- -Inf
  for(chr in c(1:22,"X")){
    load(paste0('data/LINSIGHT/LINSIGHT_chr',chr,'.rda'))
    local_min <- min(gr$score, na.rm = TRUE)
    local_max <- max(gr$score, na.rm = TRUE)
    print(paste0("chr",chr, " : ", local_min, " <-> ", local_max))
    if(local_max > max_range){
      max_range <- local_max
    }
    
    if(local_min < min_range){
      max_range <- local_min
    }
  }
  return(list(min = min_range, max = max_range))
}

extract_BAYESDEL_range <- function(){
  min_range <- +Inf
  max_range <- -Inf
  #cadd_data_dir <- "/Users/nekomimi/Workspace/vexor/vexor/data/cadd_scores/"
  #path_to_victor <- "/Users/nekomimi/Workspace/Exomes/softs/VICTOR/vAnnBase"
  cadd_data_dir <- "~/Transit/data/cadd_scores/"
  path_to_victor <- "/home/lemaud01/thesis/workspace/exomes/softs/VICTOR/vAnnBase"
  
  RES <- mclapply(X = c(1:22,"X","Y","MT"), mc.cores = 8, FUN = function(chr) {
    filename <- tempfile(tmpdir = "~/Transit/data/bayesdel", fileext = ".vcf")
    gr <- readRDS(file = paste0(cadd_data_dir,"cadd",chr,".RDS"))
    print(length(gr))
    # create a VCF
    vcf_desc_0 <- "##fileformat=VCFv4.1"
    vcf_desc_1 <- paste0("##fileDate=",gsub(x = Sys.Date(), pattern = "-",replacement = ""))
    vcf_desc_2 <- '##INFO=<ID=MAF,Number=A,Type=Float,Description="MAF">'
    vcf_header <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    vcf_content <- paste(rep(seqlevelsInUse(gr), times = length(gr)), start(gr), ".", gr$Ref, gr$Alt, 100, "PASS", paste0("MAF:"), sep = "\t")
    
    write(x = vcf_desc_0, file = filename)
    write(x = vcf_desc_1, file = filename, append = T)
    write(x = vcf_desc_2, file = filename, append = T)
    write(x = vcf_header, file = filename, append = T)
    write(x = vcf_content, file = filename, append = T)
    
    value <- extractBayesDel(path_to_victor, filename)
    
    local_min <- min(value, na.rm = TRUE)
    local_max <- max(value, na.rm = TRUE)
    print(paste0("chr",chr, " : ", local_min, " <-> ", local_max))
    if(local_max > max_range){
      max_range <- local_max
    }
    
    if(local_min < min_range){
      min_range <- local_min
    }
    
    return(list(min = min_range, max = max_range))
  })
  
  return(RES)
  
}
