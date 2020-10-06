calculateMotifScores <- function(countsDob, referenceMotifs, colData, 
                                 returnZ = FALSE,
                                 returnObject = FALSE,
                                 expect = TRUE) {
  
  atacseqCounts <- 
    dba.peakset(countsDob, bRetrieve = TRUE) %>% 
    data.frame()
  
  atacseqCounts <- as.matrix(atacseqCounts[, colData$specimenID])
  
  atacseqRowRanges <- 
    dba.peakset(countsDob, bRetrieve = TRUE)
  
  # Creating a RangedSummariedExperiment, adding GC bias, filtering peaks
  atacRSE <-
    SummarizedExperiment(assays = list(counts = atacseqCounts),
                         rowRanges = atacseqRowRanges,
                         colData = colData)
  
  atacRSE <- addGCBias(atacRSE, genome = BSgenome.Hsapiens.UCSC.hg38)
  atacRSE <- filterPeaks(atacRSE, non_overlapping = TRUE)
  
  # Finding motifs in peaks, computing accessibliity deviations
  motifsInPeaks <- matchMotifs(referenceMotifs, atacRSE,
                               genome = BSgenome.Hsapiens.UCSC.hg38)
  
  if(expect) {
    expected <- computeExpectations(atacRSE, norm = TRUE)
    motifDev <-
      computeDeviations(atacRSE,
                        annotations = motifsInPeaks,
                        expectation = expected)
  } else {
    motifDev <-
      computeDeviations(atacRSE,
                        annotations = motifsInPeaks)
  }
  
  if(!returnZ) {
    # getting raw scores of deviations
    motifScores <- deviations(motifDev)
  } else {
    # getting z-scores of deviations
    motifScores <- deviationScores(motifDev)
  }
  
  if(!returnObject) {
    return(motifScores)
  } else {
    return(list("object" = motifDev, 
                "scores" = motifScores))
  }
}

collapseReps <- function(x, sA) {
  mat_medians <- sapply(unique(sA$experimentalCondition), function(X) {
    Xids <-
      sA %>% 
      filter(experimentalCondition == X) %>% 
      pull(specimenID)
    Xmat <- x[, Xids]    # matrix of one condition's replicates
    if (length(Xids) > 1) {
      Xmat_median <- apply(Xmat, 1, median)
    } else {
      Xmat_median <- Xmat
    }
    return(Xmat_median)
  })
  colnames(mat_medians) <- unique(sA$experimentalCondition)
  mat_medians <- data.frame(mat_medians)
  return(mat_medians)
}

getJasparMotifs

jaspar <- function (...) 
  # This function is edited from the getJasparMotifs function in 
  # the R package chromVAR.
  # It hasbeen changed to retrieve JASPAR 2018 motifs, rather than
  # the JASPAR 2016 motifs retrieved by the function in chromVAR
{
  opts <- list()
  opts["species"] <- "Homo sapiens"
  opts["collection"] <- "CORE"
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
    names(out) <- paste(names(out), TFBSTools::name(out), 
                        sep = "_")
  return(out)
}