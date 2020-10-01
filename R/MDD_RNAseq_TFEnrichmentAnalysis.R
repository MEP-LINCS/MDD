# This script uses the CHEA3 API to test for TF enrichment in the RNAseq features from
# the MCF10A integrated clusters.
# Adapted from script from Mark Dane.
library(tidyverse)
# library(xlsx)
library(openxlsx)
library(readxl)
options(java.parameters = "-Xmx8000m")

###############################################################################
# Functions
# loadExcelFile <- function(excelFileName) {
#   wb <- loadWorkbook(excelFileName)
#   sheets <- names(getSheets(wb))
#   
#   myTables <- lapply(sheets, function(x) {
#     X <- read.xlsx2(excelFileName, sheetName = x)
#     X
#   })
#   
#   if (any(grepl("EGF", sheets))) {
#     names(myTables) <- str_remove(sheets, "module_")
#   } else {
#     names(myTables) <- sheets
#   }
#   return(myTables)
# }

loadExcelFile <- function(excelFileName) {
  sheets <- excel_sheets(excelFileName)
  
  myTables <- lapply(1:length(sheets), function(x) {
    X <- read_excel(excelFileName, sheet = x)
    X
  })
  
  if (any(grepl("EGF", sheets))) {
    names(myTables) <- str_remove(sheets, "module_")
  } else {
    names(myTables) <- sheets
  }
  return(myTables)
}

getTFs <- function(genes, qname){
  library(httr)
  library(jsonlite)
  
  url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
  encode = "json"
  payload = list(query_name = qname, gene_set = genes)
  
  #POST to ChEA3 server
  response = POST(url = url, body = payload, encode = encode)
  json = content(response, "text")
  
  #results as list of R dataframes
  results <- try(jsonlite::fromJSON(json))
  return(results)
}

getTFList <- function(geneList) {
  CHEA3Results <- mapply(getTFs,
                         genes = geneList,
                         qname = names(geneList),
                         SIMPLIFY = FALSE)
  return(CHEA3Results)
}

# writeListToExcel <- function(x, fname) {
#   nm <- names(x)
#   
#   write.xlsx2(x[[1]], fname, nm[1], row.names = FALSE)
#   
#   mapply(write.xlsx2, 
#          x = x[-1], 
#          sheetName = nm[-1],
#          file = fname,
#          row.names = FALSE,
#          append = TRUE)
# }

# writeListToExcel <- function(x, fname) {
#   mapply(write.xlsx, 
#          x = x, 
#          sheetName = names(x),
#          file = fname,
#          row.names = FALSE)
# }

writeListToExcel <- function(x, fname) {
  write.xlsx(x = x,
             file = fname,
             row.names = FALSE)
}


subsetListToLibrary <- function(CHEA3List, subset_to) {
  X <- lapply(CHEA3List, function(x) {
    x <- x[[subset_to]]
    return(x)
  })
  X
}

subsetListToLibraryAndWrite <- function(CHEA3List, subset_to, outDirRoot, FEATURE) {
  X <- subsetListToLibrary(CHEA3List = CHEA3List,
                           subset_to = subset_to)
  
  X_CSV <- Reduce(bind_rows, X)
  
  subset_to_clean <- str_split_fixed(subset_to, "--", 2)[1]
  
  outDirLibrary <- sprintf("%s/%s", outDirRoot, subset_to_clean)
  
  if(!dir.exists(outDirLibrary)) {
    dir.create(outDirLibrary, recursive = TRUE)
  }
  
  writeListToExcel(X, 
                   fname = sprintf("%s/MDD_%s_%s.xlsx", outDirLibrary, FEATURE, subset_to_clean))
  
  write_csv(X_CSV, sprintf("%s/MDD_%s_%s.csv", outDirLibrary, FEATURE, subset_to_clean))
}

writeListOfListsToExcel <- function(CHEA3List, outDirRoot, FEATURE) {
  outDirAll <- sprintf("%s/%s", outDirRoot, "All_Libraries")
  
  if(!dir.exists(outDirAll)) {
    dir.create(outDirAll, recursive = TRUE)
  }
  
  mapply(writeListToExcel,
         x = CHEA3List, 
         fname = sprintf("%s/MDD_%s_%s.xlsx", 
                         outDirAll, FEATURE, names(CHEA3List)))
}


combineModules <- function(moduleList, modules_to_combine) {
  
  toCombine <- modules_to_combine[sapply(modules_to_combine, length) > 1]
  toRename  <- modules_to_combine[sapply(modules_to_combine, length) == 1]
  
  combinedModules <- lapply(toCombine, function(x) {
    print(sprintf("combining %s and %s", x[1], x[2]))
    X <- suppressWarnings(
      bind_rows(moduleList[[ x[1] ]],
                moduleList[[ x[2] ]])
    )
    X
  })
  
  moduleListRename <- moduleList[unlist(toRename)]
  names(moduleListRename) <- names(toRename)
  moduleListFull <- c(combinedModules, moduleListRename)
  moduleOrder <- as.numeric(str_extract(names(moduleListFull), "[0-9]+"))
  
  moduleListFull <- moduleListFull[order(moduleOrder)]

  moduleListFull <- 
    mapply(function(x, x_nm) {
    x <- x %>% 
      mutate(clusterCombined = x_nm,
             clusterOriginal  = Cluster) %>% 
      dplyr::select(-Cluster) %>% 
      dplyr::select(feature_type, contains("cluster"), everything())
    x
  }, 
  x = moduleListFull, x_nm = names(moduleListFull),
  SIMPLIFY = FALSE)
  
  return(moduleListFull)
}
  
moduleFileToCHEA3 <- function(excelFileName,
                              outDirRoot,
                              librariesToExport = c("Literature--ChIP-seq",
                                                    "ReMap--ChIP-seq",
                                                    "ENCODE--ChIP-seq",
                                                    "Integrated--meanRank",
                                                    "Integrated--topRank"),
                              FDR_THRESHOLD = 0.25,
                              min_set_size = NULL,
                              FEATURE,
                              FEATURENAME,
                              FEATURETYPE,
                              modulesToCombine = NULL,
                              RmarkdownVersion = NULL) {
  print("loading modules")
  moduleList <- loadExcelFile(excelFileName)
  
  if (!is.null(modulesToCombine)) {
    print("joining modules...")
    moduleList <- combineModules(moduleList, 
                                 modules_to_combine = modulesToCombine)
  }
  
  moduleListLong <- Reduce(bind_rows, moduleList)
  
  if (!dir.exists(outDirRoot)) {
    dir.create(outDirRoot, recursive = TRUE)
  }
  
  print("writing modules")
  
  write_csv(moduleListLong, 
            path = sprintf("%s/MDD_%s_allFeatures.csv", 
                           outDirRoot, 
                           str_remove(FEATURE, "_CHEA3"))
            )

  moduleListLongRNAseq <- 
    moduleListLong %>% 
    filter(Type == "RNAseq")
  
  write_csv(moduleListLongRNAseq, 
            path = sprintf("%s/MDD_%s_RNAseqFeatures.csv", 
                           outDirRoot, 
                           str_remove(FEATURE, "_CHEA3"))
            )
  
  print("selecting RNAseq features")
  moduleRNAseqFeatures <- lapply(moduleList, function(x) {
    x <- x %>%
      filter(Type == "RNAseq") %>%
      distinct(feature) %>%
      pull(feature) %>%
      as.character()
    x
  })

  print("submitting to CHEA3 API")
  moduleCHEA3ResultsList <- getTFList(moduleRNAseqFeatures)

  if (is.null(RmarkdownVersion)) {
    # save results for selected libraries
    print("subsetting results to selected libraries and writing")
    lapply(librariesToExport, function(x) {
      subsetListToLibraryAndWrite(moduleCHEA3ResultsList,
                                  x,
                                  outDirRoot,
                                  FEATURE = FEATURE)
      }
    )

    # Save results for all libraries
    print("writing all libraries")

    writeListOfListsToExcel(moduleCHEA3ResultsList, outDirRoot, FEATURE = FEATURE)

    moduleCHEA3ResultsListLong <- Reduce(bind_rows, moduleCHEA3ResultsList)

    write_csv(moduleCHEA3ResultsListLong, path = sprintf("%s/MDD_%s_all.csv", outDirRoot,
                                                         FEATURE))

    save(moduleList, moduleCHEA3ResultsList, moduleRNAseqFeatures,
         file = sprintf("%s/MDD_RNAseq_%s.Rdata", outDirRoot, FEATURE))

    print("creating Rmarkdown document")

    rmarkdown::render("MDD_RNAseq_showCHEA3.Rmd",
                      output_file = sprintf("%s/MDD_%s_summary.html", outDirRoot, FEATURE))
  } else {
    
    print("creating Rmarkdown document")
    rmarkdown::render("MDD_RNAseq_showCHEA3_geneListLengths.Rmd",
                      output_file = sprintf("%s/MDD_%s_summary.html", outDirRoot, FEATURE))
    
    print(sprintf("filtering for set size length > %s", min_set_size))
    lengthFilt <- sapply(moduleRNAseqFeatures, length) > min_set_size

    print(lengthFilt)

    moduleCHEA3ResultsList <- moduleCHEA3ResultsList[lengthFilt]

    print("subsetting results to selected libraries and writing")
    lapply(librariesToExport, function(x) {
      subsetListToLibraryAndWrite(moduleCHEA3ResultsList,
                                  x,
                                  outDirRoot,
                                  FEATURE = FEATURE)
      }
      )

    # Save results for all libraries
    print("writing all libraries")

    writeListOfListsToExcel(moduleCHEA3ResultsList, outDirRoot, FEATURE = FEATURE)

    moduleCHEA3ResultsListLong <- Reduce(bind_rows, moduleCHEA3ResultsList)

    write_csv(moduleCHEA3ResultsListLong, path = sprintf("%s/MDD_%s_all.csv", outDirRoot,
                                                         FEATURE))

    save(moduleList, moduleCHEA3ResultsList, moduleRNAseqFeatures,
         file = sprintf("%s/MDD_RNAseq_%s.Rdata", outDirRoot, FEATURE))
  }
}


cTableFile <- "../RNAseq/misc/MDD_RNAseq_combinedTable.csv"

combinedTable <- read.csv(cTableFile, stringsAsFactors = FALSE)

combinedTable %>% 
  filter(Induced) %>% 
  distinct(Ligand, hgnc_symbol) %>% 
  split(.$Ligand) %>% 
  lapply(function(x) {x <- x$hgnc_symbol})
#   
# combinedTableToCHEA3 <- function(combinedTableFileName,
#                                  outDirRoot,
#                                  librariesToExport = c("Literature--ChIP-seq",
#                                                        "ReMap--ChIP-seq",
#                                                        "ENCODE--ChIP-seq",
#                                                        "Integrated--meanRank",
#                                                        "Integrated--topRank"),
#                                  FDR_THRESHOLD = 0.20,
#                                  min_set_size = 30,
#                                  FEATURE,
#                                  FEATURENAME,
#                                  FEATURETYPE,
#                                  RmarkdownVersion = NULL) {
#   print("loading combined table")
#   
#   moduleListLong <- read.csv(combinedTableFileName, stringsAsFactors = FALSE) %>% 
#     filter(Induced) %>% 
#     distinct(Ligand, hgnc_symbol, Induced)
#   
#   if (!dir.exists(outDirRoot)) {
#     dir.create(outDirRoot, recursive = TRUE)
#   }
#   
#   # write_csv(moduleListLongRNAseq, 
#   #           path = sprintf("%s/MDD_%s_.csv", 
#   #                          outDirRoot, 
#   #                          str_remove(FEATURE, "_CHEA3"))
#   # )
#   
#   print("selecting RNAseq features")
#   # moduleRNAseqFeatures <- lapply(moduleList, function(x) {
#   #   x <- x %>%
#   #     filter(Type == "RNAseq") %>%
#   #     distinct(feature) %>%
#   #     pull(feature) %>%
#   #     as.character()
#   #   x
#   # })
#   
#   moduleRNAseqFeatures <-
#     combinedTable %>% 
#     filter(Induced) %>% 
#     distinct(Ligand, hgnc_symbol) %>% 
#     split(.$Ligand) %>% 
#     lapply(function(x) {x <- x$hgnc_symbol})
#   
#   
#   lenFilt <- sapply(moduleRNAseqFeatures, length) > min_set_size
#   
#   print("removing ligands with too few genes:")
#   print(names(moduleRNAseqFeatures)[!lenFilt])
#   
#   moduleRNAseqFeatures <- moduleRNAseqFeatures[lenFilt]
#   
#   print("submitting to CHEA3 API")
#   moduleCHEA3ResultsList <- getTFList(moduleRNAseqFeatures)
#   
#   print("creating Rmarkdown document")
#   rmarkdown::render("MDD_RNAseq_showCHEA3_geneListLengths.Rmd",
#                     output_file = sprintf("%s/MDD_%s_summary.html", outDirRoot, FEATURE))
#   
#   print("subsetting results to selected libraries and writing")
#   
#   lapply(librariesToExport, function(x) {
#     subsetListToLibraryAndWrite(moduleCHEA3ResultsList,
#                                 x,
#                                 outDirRoot,
#                                 FEATURE = FEATURE)
#     }
#     )
#     
#     # Save results for all libraries
#     print("writing all libraries")
#     
#     writeListOfListsToExcel(moduleCHEA3ResultsList, outDirRoot, FEATURE = FEATURE)
#     
#     moduleCHEA3ResultsListLong <- Reduce(bind_rows, moduleCHEA3ResultsList)
#     
#     write_csv(moduleCHEA3ResultsListLong, path = sprintf("%s/MDD_%s_all.csv", outDirRoot,
#                                                          FEATURE))
#     
#     save(moduleList, moduleCHEA3ResultsList, moduleRNAseqFeatures,
#          file = sprintf("%s/MDD_RNAseq_%s.Rdata", outDirRoot, FEATURE))
#     }
# 


combinedTableToCHEA3 <- function(combinedTableFileName,
                                 outDirRoot,
                                 librariesToExport = c("Literature--ChIP-seq",
                                                       "ReMap--ChIP-seq",
                                                       "ENCODE--ChIP-seq",
                                                       "Integrated--meanRank",
                                                       "Integrated--topRank"),
                                 FDR_THRESHOLD = 0.20,
                                 min_set_size = 30,
                                 FEATURE,
                                 FEATURENAME,
                                 FEATURETYPE,
                                 RmarkdownVersion = NULL) {
  print("loading combined table")
  
  moduleListLong <- read.csv(combinedTableFileName, stringsAsFactors = FALSE) %>% 
    filter(Induced) %>% 
    distinct(Ligand, hgnc_symbol, Induced)
  
  if (!dir.exists(outDirRoot)) {
    dir.create(outDirRoot, recursive = TRUE)
  }
  
  # write_csv(moduleListLongRNAseq, 
  #           path = sprintf("%s/MDD_%s_.csv", 
  #                          outDirRoot, 
  #                          str_remove(FEATURE, "_CHEA3"))
  # )
  
  print("selecting RNAseq features")
  # moduleRNAseqFeatures <- lapply(moduleList, function(x) {
  #   x <- x %>%
  #     filter(Type == "RNAseq") %>%
  #     distinct(feature) %>%
  #     pull(feature) %>%
  #     as.character()
  #   x
  # })
  
  moduleRNAseqFeatures <-
    combinedTable %>% 
    filter(Induced) %>% 
    distinct(Ligand, hgnc_symbol) %>% 
    split(.$Ligand) %>% 
    lapply(function(x) {x <- x$hgnc_symbol})
  
  
  lenFilt <- sapply(moduleRNAseqFeatures, length) > min_set_size
  
  print("removing ligands with too few genes:")
  print(names(moduleRNAseqFeatures)[!lenFilt])
  
  moduleRNAseqFeatures <- moduleRNAseqFeatures[lenFilt]
  
  print("submitting to CHEA3 API")
  moduleCHEA3ResultsList <- getTFList(moduleRNAseqFeatures)
  # return(moduleCHEA3ResultsList)
  print("creating Rmarkdown document")
  rmarkdown::render("MDD_RNAseq_showCHEA3_geneListLengths.Rmd",
                    output_file = sprintf("%s/MDD_%s_summary.html", outDirRoot, FEATURE))

  print("subsetting results to selected libraries and writing")

  lapply(librariesToExport, function(x) {
    subsetListToLibraryAndWrite(moduleCHEA3ResultsList,
                                x,
                                outDirRoot,
                                FEATURE = FEATURE)
  }
  )

  # Save results for all libraries
  print("writing all libraries")

  writeListOfListsToExcel(moduleCHEA3ResultsList, outDirRoot, FEATURE = FEATURE)

  moduleCHEA3ResultsListLong <- Reduce(bind_rows, moduleCHEA3ResultsList)

  write_csv(moduleCHEA3ResultsListLong, path = sprintf("%s/MDD_%s_all.csv", outDirRoot,
                                                       FEATURE))

  save(moduleListLong, moduleCHEA3ResultsList, moduleRNAseqFeatures,
       file = sprintf("%s/MDD_RNAseq_%s.Rdata", outDirRoot, FEATURE))
}

# outDirRoot = "../RNAseq/Data/Enrichment/ligandUniqueFeatures_CHEA3"
# FEATURE = "ligandUniqueFeatures_CHEA3"
# FEATURENAME = "ligand-unique RNAseq features from integrated modules"
# FEATURETYPE = "ligand"
# FDR_THRESHOLD = 0.25
# 
# moduleFileToCHEA3(excelFileName = "../misc/MDD_int_rr_unique_tables.xlsx",
#                   outDirRoot = "../RNAseq/Data/Enrichment/ligandUniqueFeatures_CHEA3",
#                   FEATURE = "ligandUniqueFeatures_CHEA3",
#                   FEATURENAME = "ligand-unique RNAseq features from integrated modules",
#                   FEATURETYPE = "ligand",
#                   FDR_THRESHOLD = 0.25)

###############################################################################
# Read cluster feature files

if(!grepl("R$", getwd())) { setwd("R") }

# mtc_Combined <- list("module_1" = "module_1",
#                            "module_2" = "module_2",
#                            "module_3" = c("module_3", "module_4"),
#                            "module_4" = "module_5",
#                            "module_5" = c("module_6", "module_16"),
#                            "module_6" = "module_7",
#                            "module_7" = c("module_8", "module_17"),
#                            "module_8" = "module_9",
#                            "module_9" = c("module_10", "module_15"),
#                            "module_10" = "module_11",
#                            "module_11" = "module_12",
#                            "module_12" = "module_13",
#                            "module_13" = "module_14",
#                            "module_14" = "module_18")

moduleFileToCHEA3(excelFileName = "../misc/MDD_int_rr_combined14_tables.xlsx",
                  outDirRoot = "../RNAseq/Data/Enrichment/combined14Module_CHEA3",
                  FEATURE = "combined14Module_CHEA3",
                  FEATURENAME = "combined module RNAseq features",
                  FEATURETYPE = "module",
                  FDR_THRESHOLD = 0.20)

moduleFileToCHEA3(excelFileName = "../misc/MDD_int_rr_unique_tables.xlsx",
                  outDirRoot = "../RNAseq/Data/Enrichment/ligandUniqueFeatures_CHEA3",
                  FEATURE = "ligandUniqueFeatures_CHEA3",
                  FEATURENAME = "ligand-unique RNAseq features from integrated modules",
                  FEATURETYPE = "ligand",
                  FDR_THRESHOLD = 0.20,
                  min_set_size = 30,
                  RmarkdownVersion = "DIFFERENT")

combinedTableToCHEA3(combinedTableFileName = "../RNAseq/misc/MDD_RNAseq_combinedTable.csv",
                     outDirRoot = "../RNAseq/Data/Enrichment/inducedGene_CHEA3",
                     FEATURE = "inducedGene_CHEA3",
                     FEATURENAME = "ligand-induced RNAseq genes",
                     FEATURETYPE = "ligand",
                     min_set_size = 30,
                     FDR_THRESHOLD = 0.20)

FEATURE = "combined14Module_CHEA3"
FEATURENAME = "combined module RNAseq features"
FET_THRESHOLD = .05
outDir <- "../RNAseq/Data/Enrichment/combined14Module_CHEA3/binaryHeatmaps"
CHEA3_in  <- read.csv("../RNAseq/Data/Enrichment/combined14Module_CHEA3/ReMap/MDD_combined14Module_CHEA3_ReMap.csv",
                      stringsAsFactors = FALSE) %>% 
  mutate(significant = FET.p.value < FET_THRESHOLD) %>% 
  dplyr::rename(Module = Query.Name) %>% 
  mutate(Module = as.factor(Module)) %>% 
  mutate(Module = fct_inorder(Module))
size_for_allSig = 20
if(!dir.exists(outDir)) {dir.create(outDir)}
rmarkdown::render("MDD_RNAseq_TF_Enrichment_binaryHeatmaps.Rmd", 
                  output_file = sprintf("%s/MDD_%s_binaryHeatmaps.html", outDir, FEATURE))



FEATURE = "ligandUniqueFeatures_CHEA3"
FEATURETYPE = "ligand"
FEATURENAME = "ligand-unique RNAseq features from integrated modules"
FET_THRESHOLD = .05
outDir <- "../RNAseq/Data/Enrichment/ligandUniqueFeatures_CHEA3/binaryHeatmaps"
CHEA3_in  <- read.csv("../RNAseq/Data/Enrichment/ligandUniqueFeatures_CHEA3/ReMap/MDD_ligandUniqueFeatures_CHEA3_ReMap.csv",
                      stringsAsFactors = FALSE) %>% 
  mutate(significant = FET.p.value < FET_THRESHOLD) %>% 
  dplyr::rename(Module = Query.Name) %>% 
  mutate(Module = as.factor(Module)) %>% 
  mutate(Module = fct_inorder(Module)) %>% 
  mutate(Direction = str_extract(Module, "_[:alnum:]+$")) %>% 
  mutate(Direction = str_remove(Direction, "_")) %>% 
  mutate(Ligand = str_remove(Module, "_[:alnum:]+$")) %>% 
  mutate(annotation = case_when(significant & Direction == "Negative" ~ "Negative",
                                significant & Direction == "Positive" ~ "Positive",
                                !significant ~ ""))
size_for_allSig = 20
if(!dir.exists(outDir)) {dir.create(outDir)}
rmarkdown::render("MDD_RNAseq_TF_Enrichment_binaryHeatmaps_multiDirection.Rmd", 
                  output_file = sprintf("%s/MDD_%s_binaryHeatmaps.html", outDir, FEATURE))

