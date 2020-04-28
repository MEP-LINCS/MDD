makeDBSampleSheet <- function(sampleAnno) {
  # renames the columns of a sample annotation dataframe
  # to create a new dataframe that can be used as input to Diffbind.
  ss <- sampleAnno %>% 
    dplyr::rename("SampleID"   = specimenID,
                  "Collection" = collection,
                  "Time"       = experimentalTimePoint) %>% 
    mutate(Replicate           = as.factor(replicate),
           Factor              = fct_inorder(as.factor(experimentalCondition)),
           PeakCaller          = "narrow",
           Peaks               = paste0("../", Peaks),
           bamReads            = paste0("../", bamReads)
    )
  ss
}

bufferCoords <- function(START, END, BUFFER) {
  START <- as.numeric(START)
  END   <- as.numeric(END)
  
  if (length(BUFFER) == 2) {
    START <- START - BUFFER[1]
    END   <- END   + BUFFER[2]
  } else {
    START <- START - BUFFER
    END   <- END   + BUFFER
  }
  return(c(START, END))
}

makeGranges <- function(CHR, START, END, GENOME = "hg38", BINWIDTH = NULL) {
  
  if (!is.null(BINWIDTH)) {
    END   <- START + BINWIDTH*round((END - START)/BINWIDTH)
    STARTS <- seq(START, END - BINWIDTH, by = BINWIDTH) + 1
    ENDS   <- seq(START + BINWIDTH, END, by = BINWIDTH)
    DF <- data.frame(chr = CHR, start = STARTS, end = ENDS)
  } else {
    DF <- data.frame(chr = CHR, start = START, end = END)
  }
  
  GR <- makeGRangesFromDataFrame(DF)
  genome(GR) <- GENOME
  return(GR)
}

makeGrangesAndCountSplit <- function(dob, chr, start, end, binwidth = 20, score = DBA_SCORE_RPKM) {
  gr_temp   <- makeGranges(CHR = chr, START = start, 
                           END = end, BINWIDTH = binwidth)
  
  which <- as.character(unique(dba.show(dob)$Factor))
  
  dobsByFactor <-
    lapply(which, function(x) {
      print(x)
      dobX <- dba(dob, mask = dob$masks[[x]])
      print(dba.show(dobX))
      return(dobX)
  })
  
  dobGranges     <- lapply(dobsByFactor, dba.count, minOverlap = 1,
                           peaks = gr_temp, score = score)
  grangesCounted <- lapply(dobGranges,
                           dba.peakset, bRetrieve = TRUE)
  
  names(grangesCounted) <- which
  return(grangesCounted)
}

makeGrangesFromGene <- function(symbol, buffer, basicDob, mart, aTrna = NULL) {
  
  if(is.null(aTrna)) {
    aTrna <- getBM(attributes = c("ensembl_gene_id",
                                  "hgnc_symbol",
                                  "chromosome_name",
                                  "start_position",
                                  "end_position"), 
                   mart = mart,
                   filters = "hgnc_symbol",
                   values = symbol) %>% 
      mutate(chromosome_name = paste0("chr", chromosome_name))
  }
  
  if(nrow(aTrna) > 1) {
    stop("biomaRt returned multiple matches to symbol")
  }
  
  range <- bufferCoords(aTrna$start_position, aTrna$end_position, buffer)
  chr   <- aTrna$chromosome_name
  start <- range[1]
  end   <- range[2]
  
  geneGranges <- makeGrangesAndCountSplit(dob = basicDob,
                                         chr = chr,
                                         start = start,
                                         end = end,
                                         binwidth = 20)
  return(geneGranges)
}

makeGrangesFromTranscript <- function(transcript, buffer, basicDob, mart, aTrna = NULL) {
  
  if(is.null(aTrna)) {
    aTrna <- getBM(attributes = c("ensembl_gene_id",
                                  "ensembl_transcript_id",
                                  "hgnc_symbol",
                                  "chromosome_name",
                                  "start_position",
                                  "end_position"), 
                   mart = mart,
                   filters = "ensembl_transcript_id",
                   values = transcript) %>% 
      mutate(chromosome_name = paste0("chr", chromosome_name))
  }
  
  if(nrow(aTrna) > 1) {
    stop("biomaRt returned multiple matches to symbol")
  }
  
  range <- bufferCoords(aTrna$start_position, aTrna$end_position, buffer)
  chr   <- aTrna$chromosome_name
  start <- range[1]
  end   <- range[2]
  
  transcriptGranges <- makeGrangesAndCountSplit(dob = basicDob,
                                          chr = chr,
                                          start = start,
                                          end = end,
                                          binwidth = 20)
  return(transcriptGranges)
}


makeGeneTrackFromGene <- function(symbol, mart, grayCol = "gray70") {
  biomTrack <- BiomartGeneRegionTrack(genome="hg38",
                                      biomart = mart,
                                      filters=list(hgnc_symbol=symbol),
                                      name=symbol,
                                      collapseTranscripts = "longest",
                                      col = "black",
                                      col.line = "black",
                                      fill = grayCol,
                                      protein_coding = grayCol,
                                      "utr3" = grayCol,
                                      "composite" = grayCol,
                                      "utr5" = grayCol)
}

makeGeneTrackFromTranscript <- function(transcript, mart, grayCol = "gray70") {
  biomTrack <- BiomartGeneRegionTrack(genome="hg38",
                                      biomart = mart,
                                      filters=list(ensembl_transcript_id=transcript),
                                      name=transcript,
                                      # collapseTranscripts = "longest",
                                      col = "black",
                                      col.line = "black",
                                      fill = grayCol,
                                      protein_coding = grayCol,
                                      "utr3" = grayCol,
                                      "composite" = grayCol,
                                      "utr5" = grayCol)
}

geneToStacks <- function(symbol, buffer, basicDob, mart, aTrna = NULL) {
  
  colorLookup <- c("ctrl_0"  = "#7A4A2A",
                   "PBS_24" = "#8dd3c7",
                   "PBS_48" = "#8dd3c7",
                   "HGF_24"  = "#80b1d3",
                   "HGF_48"  = "#80b1d3",
                   "OSM_24"  = "#fdb462",
                   "OSM_48"  = "#fdb462",
                   "EGF_24"  = "#FB8072",
                   "EGF_48"  = "#FB8072",
                   "BMP2_24"  = "#B3DE69",
                   "BMP2_48"  = "#B3DE69",
                   "IFNG_24" = "#BEBADA",
                   "IFNG_48" = "#BEBADA",
                   "TGFB_24" = "#FFD92F",
                   "TGFB_48" = "#FFD92F",
                   "ctrl_0"  = "#7A4A2A")
  
  biomTrack <- makeGeneTrackFromGene(symbol, mart = mart)
  grList <- makeGrangesFromGene(symbol, buffer, basicDob, mart, aTrna = aTrna)
  
  dTracks <- mapply(function(x, x_nm) {
    DTx <- DataTrack(x,
                     range             = x,
                     genome            = "hg38",
                     isPaired          = TRUE,
                     name = x_nm,
                     fill              = colorLookup[x_nm],
                     col.histogram     = colorLookup[x_nm]
                     )
    return(DTx)
  } , x = grList, x_nm = names(grList)
  )
  
  return(c(dTracks, biomTrack))
}

transcriptToStacks <- function(transcript, buffer, basicDob, mart, aTrna = NULL) {
  
  colorLookup <- c("ctrl_0"  = "#7A4A2A",
                   "PBS_24" = "#8dd3c7",
                   "PBS_48" = "#8dd3c7",
                   "HGF_24"  = "#80b1d3",
                   "HGF_48"  = "#80b1d3",
                   "OSM_24"  = "#fdb462",
                   "OSM_48"  = "#fdb462",
                   "EGF_24"  = "#FB8072",
                   "EGF_48"  = "#FB8072",
                   "BMP2_24"  = "#B3DE69",
                   "BMP2_48"  = "#B3DE69",
                   "IFNG_24" = "#BEBADA",
                   "IFNG_48" = "#BEBADA",
                   "TGFB_24" = "#FFD92F",
                   "TGFB_48" = "#FFD92F",
                   "ctrl_0"  = "#7A4A2A")
  
  biomTrack <- makeGeneTrackFromTranscript(transcript, mart = mart)
  grList <- makeGrangesFromTranscript(transcript, buffer, basicDob, mart, aTrna = aTrna)
  
  dTracks <- mapply(function(x, x_nm) {
    DTx <- DataTrack(x,
                     range             = x,
                     genome            = "hg38",
                     isPaired          = TRUE,
                     name = x_nm,
                     fill              = colorLookup[x_nm],
                     col.histogram     = colorLookup[x_nm]
    )
    return(DTx)
  } , x = grList, x_nm = names(grList)
  )
  
  return(c(dTracks, biomTrack))
}

geneToPlot <- function(symbol, buffer, basicDob, mart, sA, extraGenes = NULL) {
  aTrna <- getBM(attributes = c("ensembl_gene_id",
                                "hgnc_symbol",
                                "chromosome_name",
                                "start_position",
                                "end_position"), 
                 mart = mart,
                 filters = "hgnc_symbol",
                 values = symbol) %>% 
    mutate(chromosome_name = paste0("chr", chromosome_name))
  
  plotRange <- bufferCoords(aTrna$start_position, aTrna$end_position, 
                            buffer)
  
  stack <- geneToStacks(symbol, buffer, basicDob, mart, aTrna)
  
  axisTrack   <- GenomeAxisTrack()
  
  if(!is.null(extraGenes)) {
    extraGeneTracks <- lapply(extraGenes, 
                              makeGeneTrackFromGene, 
                              mart = mart)
    
    plotTracks(c(stack, extraGeneTracks, axisTrack),
               from = plotRange[1],
               to   = plotRange[2],
               type = "histogram",
               scale = .1)
  } else {
    plotTracks(c(stack, axisTrack),
               from = plotRange[1],
               to   = plotRange[2],
               type = "histogram",
               scale = .1)
  }
}

transcriptToPlot <- function(transcript, buffer, basicDob, mart, sA, extraGenes = NULL) {
  aTrna <- getBM(attributes = c("ensembl_gene_id",
                                "hgnc_symbol",
                                "ensembl_transcript_id",
                                "chromosome_name",
                                "start_position",
                                "end_position"), 
                 mart = mart,
                 filters = "ensembl_transcript_id",
                 values = transcript) %>% 
    mutate(chromosome_name = paste0("chr", chromosome_name))
  
  print(aTrna)
  
  plotRange <- bufferCoords(aTrna$start_position, aTrna$end_position, 
                            buffer)
  
  stack <- transcriptToStacks(transcript, buffer, basicDob, mart, aTrna)
  
  axisTrack   <- GenomeAxisTrack()
  
  if(!is.null(extraGenes)) {
    extraGeneTracks <- lapply(extraGenes, 
                              makeGeneTrackFromGene, 
                              mart = mart)
    
    plotTracks(c(stack, extraGeneTracks, axisTrack),
               from = plotRange[1],
               to   = plotRange[2],
               type = "histogram",
               scale = .1)
  } else {
    plotTracks(c(stack, axisTrack),
               from = plotRange[1],
               to   = plotRange[2],
               type = "histogram",
               scale = .1)
  }
}
