### Get candidate SNPs
require('readr')
require(visNetwork, quietly = TRUE)
require(plot3D)

inputfile <- "~/Workspace/Finemapping/networks/locus2/locus2_chr1_113948389_114948389_summary.csv"
locus_summary <- read_delim(inputfile, "\t", 
                            escape_double = FALSE, trim_ws = TRUE)
candidate_SNP <- locus_summary[(locus_summary$type_snp == "candidate_SNP"),]

# Prepare data for VEXOR
cn <- colnames(x = candidate_SNP) 
cn[cn == "chr"] <- "Chr"
cn[cn == "position"] <- "Pos"
colnames(x = candidate_SNP) <- cn
candidate_SNP$RsID <- gsub(x = candidate_SNP$snp_name_icogs, replacement = "\\1", pattern = "(rs\\d+):.*")
candidate_SNP$End_pos <- candidate_SNP$Pos
candidate_SNP[, colnames(candidate_SNP) %in% c("info_onco","info_icogs","EAF_onco","EAF_icogs")] <- apply (X = candidate_SNP[, colnames(candidate_SNP) %in% c("info_onco","info_icogs","EAF_onco","EAF_icogs")], MARGIN = 2, FUN = function(x) {as.numeric(gsub(pattern = ",", replacement = ".", x = x))})

# outputfile <- gsub(x = inputfile, pattern = "summary", replacement = "candidate")
# write_delim(x = candidate_SNP, path = outputfile, delim = '\t')
# 
# # Prepare data for ANNOVAR
# outputfile <- gsub(x = inputfile, pattern = "summary", replacement = "annovar")
# write_delim(x = candidate_SNP[,c("Chr","Pos","Pos","Ref","Alt")], path = outputfile, delim = '\t', col_names = F)

# ANNOVAR Part
# telecharger le soft : http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
# liste des ressources disponibles : perl annotate_variation.pl -downdb avdblist humandb -buildver hg19
# telecharger les resources necessaires : 
# perl annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg19 humandb
# perl annotate_variation.pl -webfrom annovar -downdb NOM_DB -buildver hg19 humandb
# NOM_DB = eigen, gwava, fathmm, cadd13, gerp++gt2, dbnsfp33a
# perl table_annovar.pl data/research.txt humandb/ -buildver hg19 -out data/myanno -remove -protocol cadd13,dbnsfp33a,eigen,fathmm,gerp++gt2,gwava -operation f,f,f,f,f,f -nastring .

# ANNOVAR Results
myanno <- read_delim("~/Workspace/Finemapping/networks/locus2/myanno.hg19_multianno.csv", "\t", escape_double = FALSE, trim_ws = TRUE)

# ADD manually LINSIGHT score
myanno$LINSIGHT <- c(0.0507033,0.105529,0.130147,0.065669,0.0544197,'.','.',0.0501699,0.0499148,0.0496126,0.050308)

# Récuperer les données de LD
build1000GenomeLDMatrixLDlink <- function(variantsID, population, tempDir) {
  require(dplyr)
  require(httr)
  rsids = variantsID[grepl(x = variantsID, pattern = 'rs')]
  snpids = variantsID[!grepl(x = variantsID, pattern = 'rs')]
  
  results = NULL
  variants <- paste(rsids, collapse = '%0A')
  tempID <- format(Sys.time(), "VXR%y%m%d%H%M%S")
  
  LDLink_url <- paste0('https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldmatrix?snps=',
                       variants,
                       '&pop=',population,
                       '&reference=',tempID,
                       '&r2_d=r2')
  req <- GET(LDLink_url)
  results_df <-  NULL
  
  if(req$status_code == 200){
    
    ld_file = paste0(tempDir, 'LD_matrix.csv')
    
    system(command = paste0('wget -nv -O ',ld_file,' http://analysistools.nci.nih.gov/LDlink/tmp/r2_',tempID,'.txt'))
    
    
    tryCatch({
      results_df <- read.csv(file = ld_file, header = TRUE, sep = '\t')
    }, error = function(err) {
      print(err)
    })
    
  }
  
  return(results_df)
}

ldmat <- build1000GenomeLDMatrixLDlink(variantsID = candidate_SNP$RsID, population = "CEU", tempDir = '~/Workspace/Finemapping/networks/locus2/')
ldmat <- data.frame(ldmat)

#### Visualisation
### step 1 : LD edges
snp_comb <- combn(x = as.character(ldmat[,1]), m = 2)
ld_values <- rep(x = 0, times = ncol(snp_comb))

for (i in 1:ncol(snp_comb)) { 
  out_i = snp_comb[,i] ; 
  ld_values[i] <- ldmat[ldmat$RS_number == out_i[1], colnames(ldmat) == out_i[2]]
}

colfunc<-colorRampPalette(c("red","yellow","springgreen"))
colormapping <- colfunc(n = length(unique(ld_values)))
names(colormapping) <- sort(unique(ld_values), decreasing = T)
ld_values_mapped <- ld_values
for(i in unique(ld_values)){
  ld_values_mapped[ld_values_mapped == i] <- colormapping[names(colormapping) == i]
}


nodes <- data.frame(id = candidate_SNP$RsID, 
                    color = NA,
                    label = candidate_SNP$RsID,
                    shape = ifelse(test = myanno$Eigen_coding_or_noncoding == "c", yes = "star", no = "dot"),
                    font.size = 30,
                    group = "Variants"
                    #color.background = ifelse(test = myanno$Eigen_coding_or_noncoding == "c", yes = "lightgreen", no = "lightblue"),
                    #color.border = ifelse(test = myanno$Eigen_coding_or_noncoding == "c", yes = "lightgreen", no = "lightblue"),
                    #color.highlight.background = ifelse(test = myanno$Eigen_coding_or_noncoding == "c", yes = "darkgreen", no = "darkblue"),
                    #color.highlight.border = ifelse(test =myanno$Eigen_coding_or_noncoding == "c", yes = "darkgreen", no = "darkblue")
)

edges <- data.frame(from = snp_comb[1,], 
                    to = snp_comb[2,],
                    width = 3,
                    dashes = FALSE,
                    color = ld_values_mapped)



# selected annotations
selected_annotations <- c("CADD13_PHRED", "Eigen", "FATHMM_noncoding", "gerp++gt2", "GWAVA_region_score", "LINSIGHT")
sub_myanno <- myanno[, colnames(myanno) %in% selected_annotations]
row.names(sub_myanno) <- candidate_SNP$RsID

# CORRELATION BETWEEN annotations
annot_corr <-  diag(x = 1, length(selected_annotations))
colnames(annot_corr) <- rownames(annot_corr) <- selected_annotations

## CADD <-> EIGEN : .64 lin
## CADD <-> FATHMM : .59 prv
## CADD <-> GERP++ : .36 lin
## CADD <-> GWAVA : .28 prv
## CADD <-> LINSIGHT : .24 lin
annot_corr["CADD13_PHRED",] <- c(1,.64,.59,.36,.28,.24)

## EIGEN <-> FATHMM : ?
## EIGEN <-> GERP++ : .43 lin
## EIGEN <-> GWAVA : ?
## EIGEN <-> LINSIGHT : .43 lin
annot_corr["Eigen",] <- c(0,1,-1,.43,-1,.43)

## FATHMM <-> GERP++ : .76 db
## FATHMM <-> GWAVA : .29 prv
## FATHMM <-> LINSIGHT : ?
annot_corr["FATHMM_noncoding",] <- c(0,0,1,.76,.29,-1)

## GERP++ <-> GWAVA : ?
## GERP++ <-> LINSIGHT : 0.06 lin
annot_corr["gerp++gt2",] <- c(0,0,0,1,-1,0.06)

## GWAVA <-> LINSIGHT : ?
annot_corr["GWAVA_region_score",] <- c(0,0,0,0,1,-1)

## prepare color palettes for annotation

sub_myanno_mapped <- sub_myanno

##CADD plus le score est élevé plus snp deletere
cadd_colpalette <- colfunc(n = length(unique(sub_myanno$CADD13_PHRED)))
names(cadd_colpalette) <- sort(unique(sub_myanno$CADD13_PHRED), decreasing = T)

cadd_values_mapped <- sub_myanno$CADD13_PHRED
for(i in unique(sub_myanno$CADD13_PHRED)){
  cadd_values_mapped[cadd_values_mapped == i] <- cadd_colpalette[names(cadd_colpalette) == i]
}

sub_myanno_mapped$CADD13_PHRED <- cadd_values_mapped

##EIGEN plus le score est élevé plus snp deletere
eigen_colpalette <- colfunc(n = length(unique(sub_myanno$Eigen)))
names(eigen_colpalette) <- sort(unique(sub_myanno$Eigen), decreasing = T)

eigen_values_mapped <- sub_myanno$Eigen
for(i in unique(sub_myanno$Eigen)){
  eigen_values_mapped[eigen_values_mapped == i] <- eigen_colpalette[names(eigen_colpalette) == i]
}

sub_myanno_mapped$Eigen <- eigen_values_mapped

##FATHMM plus le score est élevé plus snp deletere
fathmm_colpalette <- colfunc(n = length(unique(sub_myanno$FATHMM_noncoding)))
names(fathmm_colpalette) <- sort(unique(sub_myanno$FATHMM_noncoding), decreasing = T)

fathmm_values_mapped <- sub_myanno$FATHMM_noncoding
for(i in unique(sub_myanno$FATHMM_noncoding)){
  fathmm_values_mapped[fathmm_values_mapped == i] <- fathmm_colpalette[names(fathmm_colpalette) == i]
}

sub_myanno_mapped$FATHMM_noncoding <- fathmm_values_mapped


##GERP++ plus le score est élevé plus snp deletere
gerp_colpalette <- colfunc(n = length(unique(sub_myanno$`gerp++gt2`)))
names(gerp_colpalette) <- sort(unique(sub_myanno$`gerp++gt2`), decreasing = T)

gerp_values_mapped <- sub_myanno$`gerp++gt2`
for(i in unique(sub_myanno$`gerp++gt2`)){
  gerp_values_mapped[gerp_values_mapped == i] <- gerp_colpalette[names(gerp_colpalette) == i]
}

sub_myanno_mapped$`gerp++gt2` <- gerp_values_mapped

##GWAVA plus le score est élevé plus snp deletere
gwava_colpalette <- colfunc(n = length(unique(sub_myanno$GWAVA_region_score)))
names(gwava_colpalette) <- sort(unique(sub_myanno$GWAVA_region_score), decreasing = T)

gwava_values_mapped <- sub_myanno$GWAVA_region_score
for(i in unique(sub_myanno$GWAVA_region_score)){
  gwava_values_mapped[gwava_values_mapped == i] <- gwava_colpalette[names(gwava_colpalette) == i]
}

sub_myanno_mapped$GWAVA_region_score <- gwava_values_mapped


##LINSIGHT proba donc le plus petit le plus rouge
lin_colpalette <- colfunc(n = length(unique(sub_myanno$LINSIGHT)))
names(lin_colpalette) <- sort(unique(sub_myanno$LINSIGHT), decreasing = F)

lin_values_mapped <- sub_myanno$LINSIGHT
for(i in unique(sub_myanno$LINSIGHT)){
  lin_values_mapped[lin_values_mapped == i] <- lin_colpalette[names(lin_colpalette) == i]
}

sub_myanno_mapped$LINSIGHT <- lin_values_mapped


score_nodes <- data.frame(id = paste0("scores_",candidate_SNP$RsID), 
                          color = NA,
                          label = paste0("scores_",candidate_SNP$RsID),
                          shape = "database",
                          font.size = 10,
                          group = paste0("Scores_",candidate_SNP$RsID))

score_edges <-  data.frame(from = candidate_SNP$RsID, 
                           to = paste0("scores_",candidate_SNP$RsID),
                           width = 1,
                           dashes = FALSE,
                           color = "black")

subnodes <- NULL
subedges <- NULL

for(i in 1:nrow(sub_myanno_mapped)) {
  rs_annot <- sub_myanno_mapped[i,]
  rs_values <- sub_myanno[i,]
  NA_val <- which(as.character(rs_values) == ".")
  rs_annot[NA_val] <- "#DDDDDD"
  
  if(is.null(subnodes)){
    subnodes <- data.frame(id = paste0(colnames(rs_annot),"_",i),
                           color = as.character(rs_annot),
                           label = paste0(colnames(rs_annot),":",as.character(rs_values)),
                           shape = c("square"),
                           font.size = 30,
                           group = paste0("Scores_",rownames(sub_myanno_mapped)[i]))
  } else {
    subnodes <- rbind(subnodes, 
                      data.frame(id = paste0(colnames(rs_values),"_",i),
                                 color = as.character(rs_annot),
                                 label = paste0(colnames(rs_annot),":",as.character(rs_values)),
                                 shape = c("square"),
                                 font.size = 30,
                                 group = paste0("Scores_",rownames(sub_myanno_mapped)[i])))
  }
  
  if(is.null(subedges)){
    subedges <- data.frame(from = paste0("scores_",rownames(sub_myanno_mapped)[i]),
                           to = paste0(colnames(rs_annot),"_",i),
                           width = 1,
                           dashes = FALSE,
                           color = "gray")
  } else {
    subedges <- rbind(subedges, 
                      data.frame(from = paste0("scores_",rownames(sub_myanno_mapped)[i]),
                                 to = paste0(colnames(rs_annot),"_",i),
                                 width = 1,
                                 dashes = FALSE,
                                 color = "gray"))
  }
}


#### EDGES BETWEEN ANNOTATIONS
#### get subedges
b <- split(subedges, subedges$from)
corr_edges <- NULL

for(j in seq_along(b)){
  subedges_ids <- as.character(b[[j]]$to)
  annot_comb <- combn(x = subedges_ids, m = 2)
  for( k in 1:ncol(annot_comb)) {
    # find correl
    l <- annot_comb[,k]
    lcor <- gsub(x = l, pattern = "(.*)_\\d+", replacement = "\\1")
    
    if(!is.null(corr_edges)){
      corr_edges <-  rbind(corr_edges, 
                           data.frame(from = l[1],
                                      to = l[2],
                                      width = 1 + annot_corr[lcor[1],lcor[2]],
                                      dashes = ifelse(test = annot_corr[lcor[1],lcor[2]] == -1, yes = TRUE, no = FALSE),
                                      color = "black"))
    } else {
      corr_edges <- data.frame(from = l[1],
                               to = l[2],
                               width = 1 + annot_corr[lcor[1],lcor[2]],
                               dashes = ifelse(test = annot_corr[lcor[1],lcor[2]] == -1, yes = TRUE, no = FALSE),
                               color = "black")
    }
  }
}



# annot_comb <- combn(x = selected_annotations, m = 2)
# annot_values_mapped <- rep(x = 1, times = ncol(annot_comb))
# for(i in selected_annotations){
#   annot_values_mapped[annot_values_mapped == i] <- 
# }
# 
# annot_edges <- data.frame(from = annot_comb[1,], 
#                     to = annot_comb[2,],
#                     width = annot_values_mapped,
#                     color = NA)





lnodes <- data.frame(shape = c("star","dot"),
                     # color = c("lightgreen","lightblue"),
                     label = c("Coding","Non-coding"), title = "Variants")

ledges <- data.frame(color = colormapping,
                     label = names(colormapping), width = 2, title = "LD")


vn <- visNetwork(rbind(nodes, score_nodes, subnodes), rbind(edges, score_edges, subedges,corr_edges), width = "100%", height = "1000px") %>%
  visLegend(addEdges = ledges, addNodes = lnodes, useGroups = F) %>%
  visOptions(highlightNearest = TRUE, selectedBy = "group") %>%
  visClusteringByGroup(groups = paste0("Scores_",candidate_SNP$RsID), label = "") %>%
  visPhysics(solver = "forceAtlas2Based", maxVelocity = 10,
             forceAtlas2Based = list(gravitationalConstant = -300))

for (snp in candidate_SNP$RsID){
  vn <- vn %>% visGroups(groupname = paste0("Scores_",snp), background = "#97C2FC", 
                         color = "#2B7CE9", shape = "square")
}

# require(plot3D)
# colkey(side = 3, add = F, col = rev(cadd_colpalette), clim = range(as.numeric(names(cadd_colpalette))))

