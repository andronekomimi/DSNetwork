#VEXOR EVOL

require(GenomicRanges)
require(ggplot2)
require(grDevices) 


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
                             vp.name = NULL, pop = FALSE, flip = NULL, text = FALSE, updateProgress = NULL) {
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
  
  requested_variants <- unique(requested_variants)
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
    
    data <- makeGenotypes(df)
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
    
    
    ldheatmap <- LDheatmapCustom(gdat = data, 
                                 genetic.distances = variants_pos$V2,
                                 updateProgress = updateProgress, SNP.name = variants)
    
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

build_snv_edges <- function(session_values, edges_type, edges_range){
  
  edges <- NULL
  
  if(edges_type == "0") { #distance
    print("distance")
    
    my_ranges <- session_values$my_ranges
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
            new_row <- data.frame(id = paste0("edge_dist_", i, "_", 1:ncol(dist_comb)),
                                  from = dist_comb[1,], 
                                  to = dist_comb[2,],
                                  width = 3,
                                  dashes = FALSE,
                                  xvalue = apply(dist_comb, 2, function(y) width(range(x[x$query==y[1]], x[x$query==y[2]]))),
                                  color = "black", stringsAsFactors = FALSE)
          } else {
            new_row <- data.frame(id = paste0("edge_dist_", i, "_", seq_along(dist_comb)),
                                  from = dist_comb[1], 
                                  to = dist_comb[2],
                                  width = 3,
                                  dashes = FALSE,
                                  xvalue = width(range(x[x$query==dist_comb[1]], x[x$query==dist_comb[2]])),
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
    print("linkage")
    ld_results <- session_values$ld
    nodes_edges <- c()
    
    if(!is.null(ld_results)){
      for(i in 1:ncol(ld_results)){
        ld <- ld_results[,i]$data$LDmatrix
        if(!is.null(ld)){
          ld_comb <- combn(x = colnames(ld), m = 2)
          ld_values <- rep(x = NA, times = ncol(ld_comb))
          for (j in 1:ncol(ld_comb)) { 
            out_j = ld_comb[,j] ; 
            ld_values[j] <- ld[out_j[1], colnames(ld) == out_j[2]]
          }
          names(ld_values) <- apply(ld_comb, 2, function(x) paste(x, collapse = "_"))
          nodes_edges <- c(nodes_edges, ld_values)
        }
      }
      
      colfunc<-colorRampPalette(c("red","yellow","springgreen"))
      colormapping <- colfunc(n = length(unique(nodes_edges)))
      names(colormapping) <- sort(unique(nodes_edges), decreasing = T)
      
      ld_values_mapped <- nodes_edges
      for(i in unique(nodes_edges)){
        ld_values_mapped[ld_values_mapped == i] <- colormapping[names(colormapping) == i]
      }
      
      snp_comb <- unlist(strsplit(x = names(ld_values_mapped), split = "_"))
      edges <- data.frame(id = paste0("edge_ld_", 1:length(ld_values_mapped)),
                          from = snp_comb[(1:length(snp_comb) %% 2 != 0)], 
                          to = snp_comb[(1:length(snp_comb) %% 2 == 0)],
                          width = 3,
                          dashes = FALSE,
                          xvalue = nodes_edges,
                          color = ld_values_mapped, stringsAsFactors = FALSE)
    }
  }
  
  return(edges)
} 

build_snv_nodes <- function(session_values){
  
  nodes <- NULL
  
  annotations <- session_values$res
  
  if("notfound" %in% colnames(annotations)){
    annotations <- annotations[is.na(annotations$notfound),]
    annotations$notfound <- NULL
  }
  
  node_names <- unique(annotations$query)
  
  nodes <- data.frame(id = node_names, 
                      color = NA,
                      label = node_names,
                      shape = "dot",
                      font.size = 30,
                      group = "Variants"
                      #color.background = ifelse(test = myanno$Eigen_coding_or_noncoding == "c", yes = "lightgreen", no = "lightblue"),
                      #color.border = ifelse(test = myanno$Eigen_coding_or_noncoding == "c", yes = "lightgreen", no = "lightblue"),
                      #color.highlight.background = ifelse(test = myanno$Eigen_coding_or_noncoding == "c", yes = "darkgreen", no = "darkblue"),
                      #color.highlight.border = ifelse(test =myanno$Eigen_coding_or_noncoding == "c", yes = "darkgreen", no = "darkblue")
  )
  
  return(nodes)
}


build_score_nodes <- function(session_values, selected_adj_scores, selected_raw_scores){
  nodes <- as.character(session_values$annotations$query)
  nodes_data <- data.frame(nodes = nodes)
  selected_scores <- c(selected_adj_scores, selected_raw_scores)
  
  if(length(selected_adj_scores) > 0){
    a_scores <- session_values$adjusted_scores[selected_adj_scores]
    nodes_data <- cbind(nodes_data, as.data.frame(a_scores))
  }
  
  if(length(selected_raw_scores) > 0){
    r_scores <- session_values$raw_scores[selected_raw_scores]
    nodes_data <- cbind(nodes_data, as.data.frame(r_scores))
  }
  
  # supprimer les lignes vides
  nodes_data <- nodes_data[(apply(X = nodes_data,
                                  MARGIN = 1, 
                                  FUN = function(x) sum(is.na(x))) < ncol(nodes_data) - 1), ]
  nodes_data <- data.frame(nodes_data, stringsAsFactors = FALSE)
  
  meta_values_mapped <- "blue"
  root_score_nodes <- data.frame(id = paste0("root_scores_",nodes_data$nodes), 
                            color = meta_values_mapped,
                            label = paste0("scores_",nodes_data$nodes),
                            shape = "database",
                            font.size = 10,
                            group = paste0("Scores_",nodes_data$nodes))
  
  root_score_edges <-  data.frame(id = paste0("edge_score_root_", 1:nrow(nodes_data)),
                                  from = nodes_data$nodes,
                                  to = paste0("root_scores_",nodes_data$nodes),
                                  width = 1,
                                  dashes = FALSE,
                                  xvalue = "",
                                  color = "black")
  
  if(ncol(nodes_data) > 1){
    score_nodes <- NULL
    score_edges <- NULL
    
    # les scores à proprement parler
    for(n in nodes_data$nodes){
      scores_values <- as.numeric(nodes_data[nodes_data$nodes == n,-1])
      new_n_rows <- data.frame(id = paste0(selected_scores, "_",n), 
                             color = "red",
                             label = as.character(scores_values),
                             shape = "triangle",
                             font.size = 10,
                             group = paste0("Scores_",n))
      
      if(is.null(score_nodes)){
        score_nodes <- new_n_rows
      } else {
        score_nodes <- rbind(score_nodes, new_n_rows)
      }
      
      new_e_rows <-  data.frame(id = paste0("edge_scores_", n, "_", seq_along(selected_scores)),
                                from = paste0("root_scores_",n),
                                to = paste0(selected_scores, "_",n),
                                width = 1,
                                dashes = FALSE,
                                xvalue = "",
                                color = "green")
      
      if(is.null(score_edges)){
        score_edges <- new_e_rows
      } else {
        score_edges <- rbind(score_edges, new_e_rows)
      }
    }
    
    root_score_nodes <- rbind(root_score_nodes, score_nodes)
    root_score_edges <- rbind(root_score_edges, score_edges)
  }
  
  
  return(list(nodes = root_score_nodes, edges = root_score_edges))
  
  # if("notfound" %in% colnames(annotations)){
  #   annotations <- annotations[is.na(annotations$notfound),]
  #   annotations$notfound <- NULL
  # }
  # 
  # node_names <- unique(annotations$query)
  # 
  # nodes <- data.frame(id = node_names, 
  #                     color = NA,
  #                     label = node_names,
  #                     #shape = ifelse(test = myanno$Eigen_coding_or_noncoding == "c", yes = "star", no = "dot"),
  #                     font.size = 30,
  #                     group = "Variants"
  #                     #color.background = ifelse(test = myanno$Eigen_coding_or_noncoding == "c", yes = "lightgreen", no = "lightblue"),
  #                     #color.border = ifelse(test = myanno$Eigen_coding_or_noncoding == "c", yes = "lightgreen", no = "lightblue"),
  #                     #color.highlight.background = ifelse(test = myanno$Eigen_coding_or_noncoding == "c", yes = "darkgreen", no = "darkblue"),
  #                     #color.highlight.border = ifelse(test =myanno$Eigen_coding_or_noncoding == "c", yes = "darkgreen", no = "darkblue")
  # )
  # 
  # return(nodes)
}


mean_score <- function(x){
  suppressWarnings(mean(x = as.numeric(unlist(strsplit(x = x, split= ","))), na.rm= T))
}

extract_score_and_convert <- function(annotations_infos, score_name, sub_score_name = NULL){
  if(sum(colnames(annotations_infos) == score_name) == 1){
    x <- as.character(annotations_infos[,score_name])
    x <- sapply(x, mean_score)
    return(x)
  } else {
    if(!is.null(sub_score_name)){
      if(sum(colnames(annotations_infos) == sub_score_name) == 1){
        x <- as.character(annotations_infos[,sub_score_name])
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
