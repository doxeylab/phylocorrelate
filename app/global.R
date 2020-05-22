#------------------------------------------------------------------------------
# PhyloCorrelations v1.0
# global.R
# Last modified: 2020-05-22 21:57:16 (CEST)
# BJM Tremblay

# library(profvis)

msg <- function(...) {
  time <- format(as.POSIXlt(Sys.time(), tz = "America/Toronto"))
  message(paste0(c(time, ...), collapse = " "))
}
msg("Loading global.R")

msg("Loading packages")

msg("  shiny")
library(shiny)
msg("  shinyWidgets")
library(shinyWidgets)
msg("  DT")
suppressPackageStartupMessages(library(DT))
msg("  shinycssloaders")
library(shinycssloaders)
msg("  shinyjs")
suppressPackageStartupMessages(library(shinyjs))
msg("  igraph")
suppressPackageStartupMessages(library(igraph))
msg("  plotly")
suppressPackageStartupMessages(library(plotly))
msg("  visNetwork")
library(visNetwork)
msg("  topGO")
suppressMessages(suppressPackageStartupMessages(library(topGO)))
msg("  shinythemes")
library(shinythemes)
msg("  ggplot")
library(ggplot2)
msg("  reshape2")
library(reshape2)
msg("  magrittr")
library(magrittr)
msg("  promises")
library(promises)
msg("  future")
suppressPackageStartupMessages(library(future))
plan(multiprocess)

CONFIGS <- readr::read_lines("PhyloCorrConfig.txt") %>%
  Filter(function(x) x != "", .) %>%
  Filter(function(x) !grepl("^\\s+$", x), .) %>%
  paste0(collapse = ",") %>%
  paste0("list(", ., ")") %>%
  parse(text = .) %>%
  eval()

startBlast <- function() {
  if (CONFIGS$UseBlastp) system("Rscript blastp/blastp.R &")
}

BlastpVer <- 4L

startBlast()

#------------------------------------------------------------------------------
# Data - in memory

msg("Loading data")

KOCounts <- readRDS("data/KOCounts.RDS")
KODesc <- readRDS("data/KODescriptions.RDS")
KEGGPathways <- readRDS("data/KEGGPathways.RDS")
KO2Pathway <- readRDS("data/KO2Pathway.RDS")
KOs <- names(KOCounts)
KODesc <- KODesc[names(KODesc) %in% KOs]
KO2Pathway <- KO2Pathway[names(KO2Pathway) %in% KOs]
KO2Pathway <- KO2Pathway[!grepl("ko", KO2Pathway)]
KOAnnoTree <- readRDS("data/KOAnnoTree.RDS")
KOLinks <- readRDS("data/KOLinks.RDS")
KOCombLinks <- readRDS("data/KOLinks2.RDS")

PFAMCounts <- readRDS("data/PFAMCounts.RDS")
PFAMDesc <- readRDS("data/PFAMDescriptions.RDS")
PFAMGO <- readRDS("data/PFAMGOTerms.RDS")
PFAMs <- names(PFAMCounts)
PFAMGO <- PFAMGO[PFAMGO$PFAM %in% PFAMs, ]
PFAMGObp <- PFAMGO[PFAMGO$Ontology == "BP", ]
PFAMGOmf <- PFAMGO[PFAMGO$Ontology == "MF", ]
PFAMGOcc <- PFAMGO[PFAMGO$Ontology == "CC", ]
PFAM2GObp <- readRDS("data/PFAM2GObp.RDS")
PFAM2GOcc <- readRDS("data/PFAM2GOcc.RDS")
PFAM2GOmf <- readRDS("data/PFAM2GOmf.RDS")
PFAMAnnoTree <- readRDS("data/PFAMAnnoTree.RDS")
PFAMLinks <- readRDS("data/PFAMLinks.RDS")
PFAMCombLinks <- readRDS("data/PFAMLinks2.RDS")
PFAMDesc <- PFAMDesc[names(PFAMDesc) %in% names(PFAMCounts)]

TIGRFAMCounts <- readRDS("data/TIGRFAMCounts.RDS")
TIGRFAMDesc <- readRDS("data/TIGRFAMDescriptions.RDS")
TIGRFAMGO <- readRDS("data/TIGRFAMGOTerms.RDS")
TIGRFAMs <- names(TIGRFAMCounts)
TIGRFAMGO <- TIGRFAMGO[TIGRFAMGO$TIGRFAM %in% TIGRFAMs, ]
TIGRFAMGObp <- TIGRFAMGO[TIGRFAMGO$Ontology == "BP", ]
TIGRFAMGOmf <- TIGRFAMGO[TIGRFAMGO$Ontology == "MF", ]
TIGRFAMGOcc <- TIGRFAMGO[TIGRFAMGO$Ontology == "CC", ]
TIGRFAM2GObp <- readRDS("data/TIGRFAM2GObp.RDS")
TIGRFAM2GOcc <- readRDS("data/TIGRFAM2GOcc.RDS")
TIGRFAM2GOmf <- readRDS("data/TIGRFAM2GOmf.RDS")
TIGRFAMAnnoTree <- readRDS("data/TIGRFAMAnnoTree.RDS")
TIGRFAMLinks <- readRDS("data/TIGRFAMLinks.RDS")
TIGRFAMCombLinks <- readRDS("data/TIGRFAMLinks2.RDS")

#------------------------------------------------------------------------------
# Data - on disk

getDat <- function(f, i, j, drp, rn) {
  d <- fst::read_fst(f, j)
  rownames(d) <- rn
  if (!is.null(i)) d <- d[i, , drop = FALSE]
  if (drp && length(d) == 1) structure(d[[1]], names = rn) else d
}

KOOverlaps <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/KOOverlaps.fst", i, j, drop, KOs)
}
KORunHyperP <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/KO_rJC_HyperP.fst", i, j, drop, KOs)
}
KOOccDiff <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/KOOccDiff.fst", i, j, drop, KOs)
}
KORunJaccardCoef <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/KORunJaccardCoef.fst", i, j, drop, KOs)
}
KORunLens <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/KO_RunsLen.fst", i, j, drop, KOs)
}
KOPathTDR <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/KO_PMF_OccDiffvsrHyperP.fst", i, j, drop, KOs)
}
KOJaccardCoef <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/KOJaccardCoef.fst", i, j, drop, KOs)
}
KOTable <- function(j) {
  fst::read_fst("data/KOTable.fst", j)[[1]]
}
KOBestPMF <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/KOBestPMF.fst", i, j, drop, KOs)
}

PFAMOverlaps <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/PFAMOverlaps.fst", i, j, drop, PFAMs)
}
PFAMPagelPearson <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/PFAMPagelPearson.fst", i, j, drop, PFAMs)
}
PFAMRunHyperP <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/PFAM_rJC_HyperP.fst", i, j, drop, PFAMs)
}
PFAMOccDiff <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/PFAMOccDiff.fst", i, j, drop, PFAMs)
}
PFAMRunJaccardCoef <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/PFAMRunJaccardCoef.fst", i, j, drop, PFAMs)
}
PFAMRunLens <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/PFAM_RunsLen.fst", i, j, drop, PFAMs)
}
PFAMGObpTDR <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/PFAM_PMF_OccDiffvsrHyperP.fst", i, j, drop, PFAMs)
}
PFAMJaccardCoef <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/PFAMJaccardCoef.fst", i, j, drop, PFAMs)
}
PFAMTable <- function(j) {
  fst::read_fst("data/PFAMTable.fst", j)[[1]]
}
PFAMBestPMF <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/PFAMBestPMF.fst", i, j, drop, PFAMs)
}

TIGRFAMOverlaps <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/TIGRFAMOverlaps.fst", i, j, drop, TIGRFAMs)
}
TIGRFAMPagelPearson <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/TIGRFAMPagelPearson.fst", i, j, drop, TIGRFAMs)
}
TIGRFAMRunHyperP <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/TIGRFAM_rJC_HyperP.fst", i, j, drop, TIGRFAMs)
}
TIGRFAMOccDiff <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/TIGRFAMOccDiff.fst", i, j, drop, TIGRFAMs)
}
TIGRFAMRunJaccardCoef <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/TIGRFAMRunJaccardCoef.fst", i, j, drop, TIGRFAMs)
}
TIGRFAMRunLens <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/TIGRFAM_RunsLen.fst", i, j, drop, TIGRFAMs)
}
TIGRFAMGObpTDR <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/TIGRFAM_PMF_OccDiffvsrHyperP3.fst", i, j, drop, TIGRFAMs)
}
TIGRFAMJaccardCoef <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/TIGRFAMJaccardCoef.fst", i, j, drop, TIGRFAMs)
}
TIGRFAMTable <- function(j) {
  fst::read_fst("data/TIGRFAMTable.fst", j)[[1]]
}
TIGRFAMBestPMF <- function(i = NULL, j = NULL, drop = TRUE) {
  getDat("data/TIGRFAMBestPMF.fst", i, j, drop, TIGRFAMs)
}

#------------------------------------------------------------------------------
# Server functions: misc

msg("Loading functions")

makeButtons <- function(name, terms) {
  structure(paste0(
    vapply(
      paste0("button_", terms),
      function(x) as.character(
        actionButton(x, label = "", icon = icon("columns"),
          style = "padding:0px; color: grey; background-color: white; border-color: white;",
          onclick = sprintf(
            'Shiny.onInputChange("%s_MAT_BUTTON", this.id)', name
          )
        )
      ),
      character(1)
    )
  ), names = terms)
}
PFAMButtons <- readRDS("data/PFAMButtons.RDS")
TIGRFAMButtons <- readRDS("data/TIGRFAMButtons.RDS")
KOButtons <- readRDS("data/KOButtons.RDS")

getPFAMtable <- function(entry, keepEntry = FALSE, entryOcc = 0) {
  msg("Making PFAM table for", entry)
  o <- data.frame(
    Links = PFAMCombLinks[PFAMs],
    Description = PFAMDesc[PFAMs],
    Compare = PFAMButtons,
    Ov = PFAMOverlaps(, entry),
    Occ = PFAMCounts[PFAMs],
    OccDiff = PFAMOccDiff(, entry),
    JC = round(PFAMJaccardCoef(, entry), 3),
    rJC = round(PFAMRunJaccardCoef(, entry), 3),
    rHyperP = PFAMRunHyperP(, entry),
    CS = PFAMBestPMF(, entry),
    row.names = PFAMs,
    stringsAsFactors = FALSE
  )
  o <- o[order(match(o$CS, c("Very high", "High", "Low", "Very low")), o$rHyperP), ]
  if (keepEntry) o else o[rownames(o) != entry, ]
}

getTIGRFAMtable <- function(entry, keepEntry = FALSE, entryOcc = 0) {
  msg("Making TIGRFAM table for", entry)
  o <- data.frame(
    Links = TIGRFAMCombLinks[TIGRFAMs],
    Description = TIGRFAMDesc[TIGRFAMs],
    Compare = TIGRFAMButtons,
    Ov = TIGRFAMOverlaps(, entry),
    Occ = TIGRFAMCounts[TIGRFAMs],
    OccDiff = TIGRFAMOccDiff(, entry),
    JC = round(TIGRFAMJaccardCoef(, entry), 3),
    rJC = round(TIGRFAMRunJaccardCoef(, entry), 3),
    rHyperP = TIGRFAMRunHyperP(, entry),
    CS = TIGRFAMBestPMF(, entry),
    row.names = TIGRFAMs
  )
  o <- o[order(match(o$CS, c("Very high", "High", "Low", "Very low")), o$rHyperP), ]
  if (keepEntry) o else o[rownames(o) != entry, ]
}

getKOtable <- function(entry, keepEntry = FALSE, entryOcc = 0) {
  msg("Making KO table for", entry)
  o <- data.frame(
    Links = KOCombLinks[KOs],
    Description = KODesc[KOs],
    Compare = KOButtons,
    Ov = KOOverlaps(, entry),
    Occ = KOCounts[KOs],
    OccDiff = KOOccDiff(, entry),
    JC = round(KOJaccardCoef(, entry), 3),
    rJC = round(KORunJaccardCoef(, entry), 3),
    rHyperP = KORunHyperP(, entry),
    CS = KOBestPMF(, entry),
    row.names = KOs
  )
  o <- o[order(match(o$CS, c("Very high", "High", "Low", "Very low")), o$rHyperP), ]
  if (keepEntry) o else o[rownames(o) != entry, ]
}

getPFAMinfo <- function(entry) {
  msg("Making input info for", entry)
  if (!entry %in% PFAMs) return()
  list(
    pfam = PFAMLinks[entry],
    description = PFAMDesc[entry],
    occurrences = PFAMCounts[entry],
    annotree = PFAMAnnoTree[entry]
  )
}

getTIGRFAMinfo <- function(entry) {
  msg("Making input info for", entry)
  if (!entry %in% TIGRFAMs) return()
  list(
    tigrfam = TIGRFAMLinks[entry],
    description = TIGRFAMDesc[entry],
    occurrences = TIGRFAMCounts[entry],
    annotree = TIGRFAMAnnoTree[entry]
  )
}

getKEGGinfo <- function(entry) {
  msg("Making input info for", entry)
  if (!entry %in% KOs) return()
  list(
    kegg = KOLinks[entry],
    description = KODesc[entry],
    occurrences = KOCounts[entry],
    annotree = KOAnnoTree[entry]
  )
}

searchPFAMs <- function(entry) {
  res <- PFAMDesc[grepl(entry, PFAMDesc, ignore.case = TRUE) |
                  grepl(entry, names(PFAMDesc), ignore.case = TRUE)]
  if (length(res) == 0) return()
  data.frame(row.names = names(res), Description = unname(res),
             stringsAsFactors = FALSE)
}

searchTIGRFAMs <- function(entry) {
  res <- TIGRFAMDesc[grepl(entry, TIGRFAMDesc, ignore.case = TRUE) |
                     grepl(entry, names(TIGRFAMDesc), ignore.case = TRUE)]
  if (length(res) == 0) return()
  data.frame(row.names = names(res), Description = unname(res),
             stringsAsFactors = FALSE)
}

searchKEGGs <- function(entry) {
  res <- KODesc[grepl(entry, KODesc, ignore.case = TRUE) |
                grepl(entry, names(KODesc), ignore.case = TRUE)]
  if (length(res) == 0) return()
  data.frame(row.names = names(res), Description = unname(res),
             stringsAsFactors = FALSE)
}

filter_blastp_col <- function(df, tokeep) {
  if (all(is.na(tokeep))) df else df[tokeep, ]
}

getTableFiltered <- function(entry, globals = list(), keepEntry = FALSE, isBlastp = FALSE) {
  if (!isBlastp) {
    what <- substr(entry, 1, 1)
    o <- switch(what, P = getPFAMtable(entry, keepEntry, PFAMCounts[entry]),
                K = getKOtable(entry, keepEntry, TIGRFAMCounts[entry]),
                T = getTIGRFAMtable(entry, keepEntry, KOCounts[entry]))
    o <- o[
           o$Ov >= globals$Ov &
           o$Occ >= globals$Occ &
           o$OccDiff <= globals$OccDiff &
           o$JC >= globals$JC &
           o$rHyperP <= globals$rHyperP &
           o$rJC >= globals$rJC, ]
    o <- o[o$CS %in% globals$CS, ]
    o <- o[!is.na(o$rJC), ]
  } else {
    o <- cleanBlastpTable(entry)
    o <- filter_blastp_col(o, o$Ov >= globals$Ov)
    o <- filter_blastp_col(o, o$Occ >= globals$Occ)
    o <- filter_blastp_col(o, o$OccDiff <= globals$OccDiff)
    o <- filter_blastp_col(o, o$JC >= globals$JC)
    o <- filter_blastp_col(o, o$rJC >= globals$rJC)
    o <- filter_blastp_col(o, o$rHyperP <= globals$rHyperP)
    o <- filter_blastp_col(o, o$CS %in% globals$CS)
  }
  msg("Table size:", nrow(o))
  o
}

getCorrMat <- function(entry, size, globals = list(), isBlastp = FALSE) {
  if (isBlastp) {
    msg("Making correlation matrix for blastp entry")
    d <- getTableFiltered(entry, globals, isBlastp = TRUE)
  } else {
    msg("Making correlation matrix for", entry)
    d <- getTableFiltered(entry, globals, keepEntry = TRUE)
  }
  if (nrow(d) <= 1) {
    msg("No results")
    return(matrix(nrow = 0, ncol = 0))
  }
  if (globals$metric == "CS") {
    CSlevels <- c("Very low", "Low", "High", "Very high")
    d[[globals$metric]] <- factor(
      d[[globals$metric]],
      levels = CSlevels
    )
  }
  d <- d[order(d[[globals$metric]], decreasing = TRUE), ]
  if (nrow(d) > size) d <- d[1:size, ]
  if (nrow(d) <= 1) {
    msg("No results")
    return(matrix(nrow = 0, ncol = 0))
  }
  msg("Correlation table matrix:", nrow(d), "x", ncol(d))
  toget <- rownames(d)
  if (isBlastp) {
    what <- substr(toget[1], 1, 1)
  } else {
    if (!entry %in% toget) toget <- c(toget, entry)
    what <- substr(entry, 1, 1)
  }
  out <- switch(what,
    P = switch(globals$metric,
          JC = as.matrix(PFAMJaccardCoef(toget, toget)),
          rJC = as.matrix(PFAMRunJaccardCoef(toget, toget)),
          CS = as.matrix(PFAMBestPMF(toget, toget))
        ),
    T = switch(globals$metric,
          JC = as.matrix(TIGRFAMJaccardCoef(toget, toget)),
          rJC = as.matrix(TIGRFAMRunJaccardCoef(toget, toget)),
          CS = as.matrix(TIGRFAMBestPMF(toget, toget))
        ),
    K = switch(globals$metric,
          JC = as.matrix(KOJaccardCoef(toget, toget)),
          rJC = as.matrix(KORunJaccardCoef(toget, toget)),
          CS = as.matrix(KOBestPMF(toget, toget))
        )
  )
  if (isBlastp) {
    out <- cbind(out, Query = d[rownames(out), globals$metric])
    out <- rbind(out, Query = c(d[rownames(out), globals$metric], 1))
  }
  out
}

intToCS <- function(x) {
  x <- as.integer(x)
  CSlevels <- c("Very low", "Low", "High", "Very high", "")
  x[x == 0] <- 5
  CSlevels[x]
}

makeCorrNetwork <- function(entry, globals = list(), isBlastp = FALSE) {
  if (!isBlastp) msg("Making correlation network for", entry)
  else msg("Making correlation network for blastp entry")
  m <- getCorrMat(entry, 200, globals, isBlastp = isBlastp)
  if (isBlastp) entry <- "Query"
  if (ncol(m) <= 1) {
    msg("No results")
    return()
  }
  if (isBlastp) {
    what <- substr(rownames(m)[1], 1, 1)
  } else {
    what <- substr(entry, 1, 1)
  }
  CSlevels <- c("Very low", "Low", "High", "Very high")
  if (globals$metric != "CS") {
    m[m < switch(globals$metric, JC = globals$JC, rJC = globals$rJC)] <- 0
    diag(m) <- 0
  } else {
    m <- matrix(
      as.integer(factor(as.character(m), levels = CSlevels)),
      nrow = nrow(m), dimnames = dimnames(m)
    )
    m[!intToCS(m) %in% globals$CS] <- 0
    diag(m) <- 0
  }
  m <- m[apply(m, 2, function(x) any(x > 0)), apply(m, 2, function(x) any(x > 0))]
  col_e <- m[, entry]
  row_e <- m[entry, ]
  if (globals$metric != "CS") {
    m[m < globals$min2nd] <- 0
  } else {
    m[!intToCS(m) %in% globals$cswhich] <- 0
  }
  m[, entry] <- col_e
  m[entry, ] <- row_e
  msg("Network size:", nrow(m))
  network <- graph_from_adjacency_matrix(m, weighted = TRUE, mode = "undirected",
                                         diag = FALSE)
  todel <- which(igraph::degree(network) == 0)
  if (length(todel) > 0) network <- igraph::delete.vertices(network, todel)
  data <- toVisNetworkData(network)
  data$nodes$shape <- "circle"
  data$nodes$color.border <- "black"
  data$nodes$border.width <- 2
  descs <- switch(what, P = PFAMDesc[data$nodes$id],
                  T = TIGRFAMDesc[data$nodes$id], K = KODesc[data$nodes$id])
  data$nodes$title <- paste(data$nodes$id, descs)
  if (isBlastp) {
    data$nodes$title[data$nodes$id == entry] <- entry
  }
  if (globals$metric == "CS") {
    data$edges$value <- data$edges$weight / 4
    data$edges$title <- intToCS(data$edges$weight)
  } else {
    data$edges$value <- data$edges$weight
    data$edges$title <- data$edges$weight
  }
  data$nodes$color.background <- "gold"
  data$nodes$color.background[data$nodes$id == entry] <- "tomato"
  visNetwork(nodes = data$nodes, edges = data$edges) %>%
    visNodes(physics = FALSE, borderWidthSelected = 2) %>%
    visInteraction(selectConnectedEdges = FALSE,
                   tooltipDelay = 0) %>%
    visIgraphLayout()
}

makeDistPlot <- function(entry, type, isBlastp = FALSE) {
  if (!isBlastp) {
    msg("Making distribution plot for", entry)
    what <- substr(entry, 1, 1)
    scores <- switch(what,
      P = switch(type,
            JC = PFAMJaccardCoef(, entry),
            rJC = PFAMRunJaccardCoef(, entry),
            rHyperP = PFAMRunHyperP(, entry)
          ),
      T = switch(type,
            JC = TIGRFAMJaccardCoef(, entry),
            rJC = TIGRFAMRunJaccardCoef(, entry),
            rHyperP = TIGRFAMRunHyperP(, entry)
          ),
      K = switch(type,
            JC = KOJaccardCoef(, entry),
            rJC = KORunJaccardCoef(, entry),
            rHyperP = KORunHyperP(, entry)
          ),
    )
    scores <- scores[names(scores) != entry]
  } else {
    msg("Making distribution plot for blastp entry")
    scores <- entry[[type]]
    entry <- "BLASTP query"
  }
  switch(type,
    JC = plot_ly(
           x = scores,
           type = "histogram",
           xbins = list(start = 0, end = 1, size = 0.00625)
         ) %>%
           layout(
             title = paste0("JC scores for ", entry),
             xaxis = list(
                       title = "JC scores",
                       autotick = FALSE,
                       showline = TRUE,
                       dtick = 0.25,
                       ticklen = 5,
                       tickwidth = 1,
                       range = c(0, 1)
                     ),
             yaxis = list(
                       title = "Count",
                       type = "log"
                     )
           ),
    rJC = plot_ly(
           x = scores,
           type = "histogram",
           xbins = list(start = 0, end = 1, size = 0.00625)
         ) %>%
           layout(
             title = paste0("rJC scores for ", entry),
             xaxis = list(
                       title = "rJC scores",
                       autotick = FALSE,
                       showline = TRUE,
                       dtick = 0.25,
                       ticklen = 5,
                       tickwidth = 1,
                       range = c(0, 1)
                     ),
             yaxis = list(
                       title = "Count",
                       type = "log"
                     )
           ),
    rHyperP = plot_ly(
           x = as.integer(-log10(scores)),
           type = "histogram",
           xbins = list(start = 0, end = 320, size = 1)
         ) %>%
           layout(
             title = paste0("-log10(rHyperP) scores for ", entry),
             xaxis = list(
                       title = "-log10(rHyperP) scores",
                       autotick = FALSE,
                       showline = TRUE,
                       dtick = 50,
                       ticklen = 5,
                       tickwidth = 1,
                       range = c(0, 320)
                     ),
             yaxis = list(
                       title = "Count",
                       type = "log"
                     )
           )
  )
}

makeGOtable <- function(entry, globals = list(),
                        return.go = FALSE, isBlastp = FALSE) {
  if (!isBlastp) {
    msg("Running GO enrichment for", entry)
    what <- substr(entry, 1, 1)
  } else {
    msg("Running GO enrichment for blastp entry")
    what <- substr(rownames(entry)[1], 1, 1)
  }
  x <- getTableFiltered(entry, globals, isBlastp = isBlastp)
  if (nrow(x) <= 1) {
    msg("No results")
    return(data.frame(GO.ID = character(), Term = character(), Total = integer(),
                      Present = integer(), Expected = integer(), FP = numeric(),
                      BH = numeric()))
  }
  switch(globals$Ont,
    BP = {
      switch(what,
        P = {
          g <- factor(as.integer(names(PFAM2GObp) %in% rownames(x)))
          names(g) <- names(PFAM2GObp)
          g2g <- PFAM2GObp
        },
        T = {
          g <- factor(as.integer(names(TIGRFAM2GObp) %in% rownames(x)))
          names(g) <- names(TIGRFAM2GObp)
          g2g <- TIGRFAM2GObp
        }
      )
    },
    CC = {
      switch(what,
        P = {
          g <- factor(as.integer(names(PFAM2GOcc) %in% rownames(x)))
          names(g) <- names(PFAM2GOcc)
          g2g <- PFAM2GOcc
        },
        T = {
          g <- factor(as.integer(names(TIGRFAM2GOcc) %in% rownames(x)))
          names(g) <- names(TIGRFAM2GOcc)
          g2g <- TIGRFAM2GOcc
        }
      )
    },
    MF = {
      switch(what,
        P = {
          g <- factor(as.integer(names(PFAM2GOmf) %in% rownames(x)))
          names(g) <- names(PFAM2GOmf)
          g2g <- PFAM2GOmf
        },
        T = {
          g <- factor(as.integer(names(TIGRFAM2GOmf) %in% rownames(x)))
          names(g) <- names(TIGRFAM2GOmf)
          g2g <- TIGRFAM2GOmf
        }
      )
    }
  )
  if (length(levels(g)) < 2) {
    msg("No results")
    return(data.frame(GO.ID = character(), Term = character(), Total = integer(),
                      Present = integer(), Expected = integer(), FP = numeric(),
                      BH = numeric()))
  }
  tg <- new("topGOdata", allGenes = g, ontology = globals$Ont, gene2GO = g2g,
            annot = annFUN.gene2GO)
  tr <- suppressMessages(runTest(tg, "classic", "fisher"))
  gt <- GenTable(tg, FisherPval = tr, topNodes = length(score(tr)),
                 numChar = 100)
  gt$FisherPval <- as.numeric(gt$FisherPval)
  gt$FisherAdjPval <- p.adjust(gt$FisherPval, "fdr")
  gt <- gt[!is.na(gt$FisherAdjPval), ]
  gt <- gt[gt$Annotated >= globals$minTotal, ]
  gt <- gt[gt$Significant >= globals$minPresent, ]
  gt <- gt[gt$Expected >= globals$minExp, ]
  gt <- gt[gt$Significant > gt$Exp, ]
  gt <- gt[gt$FisherAdjPval <= globals$FDR, ]
  gt <- gt[gt$FisherPval <= globals$FP, ]
  colnames(gt) <- c("GO.ID", "Term", "Total", "Present", "Expected", "FP", "BH")
  rownames(gt) <- NULL
  gt <- gt[order(gt$FP), ]
  msg("GO enrichment table size:", nrow(gt))
  if (return.go) return(list(tr, gt, tg))
  gt
}

getKeggEnrich <- function(entry, globals = list(), isBlastp = FALSE) {
  if (!isBlastp) {
    msg("Running pathway enrichment for", entry)
  } else {
    msg("Running pathway enrichment for blastp entry")
  }
  x <- getTableFiltered(entry, globals, isBlastp = isBlastp)
  if (nrow(x) <= 1) {
    msg("No results")
    return(data.frame(Pathway.ID = character(), Description = character(), 
                      Total = integer(), Present = integer(),
                      Expected = integer(), FP = numeric(), BH = numeric()))
  }
  res <- enrichKegg(rownames(x), KO2Pathway, KEGGPathways)
  res <- res[res$Total >= globals$minTotal, ]
  res <- res[res$Present >= globals$minPresent, ]
  res <- res[res$Expected >= globals$minExp, ]
  res <- res[res$FP <= globals$FP, ]
  res <- res[res$BH <= globals$FDR, ]
  rownames(res) <- NULL
  msg("Pathway enrichment table size:", nrow(res))
  res
}

enrichKegg <- function(target, ko2pathway, pathways) {

  tt <- length(unique(names(ko2pathway)))
  target <- target[target %in% names(ko2pathway)]
  if (length(target) == 0) {
    return(data.frame(Pathway.ID = character(), Description = character(), 
                      Total = integer(), Present = integer(),
                      Expected = integer(), FP = numeric(), BH = numeric()))
  }

  t_in <- table(factor(ko2pathway[names(ko2pathway) %in% target],
                levels = unique(unname(ko2pathway))))
  t_total <- table(factor(ko2pathway, levels = unique(unname(ko2pathway))))

  t_df <- data.frame(Pathway.ID = unique(unname(ko2pathway)),
                     Present = as.integer(t_in),
                     Total = as.integer(t_total),
                     stringsAsFactors = FALSE)
  t_df$Expected <- t_df$Total / tt
  t_df$Expected <- round(t_df$Expected * length(target), 3)

  t_df$FP <- mapply(
    function(x, y) {
      if (x <= 1) return(1)
      x <- x - 1
      phyper(q = x, m = y,
             n = tt - y,
             k = length(target),
             lower.tail = FALSE)
    },
    t_df$Present, t_df$Total
  )

  t_df$BH <- p.adjust(t_df$FP, "BH")
  t_df$Description <- pathways[t_df$Pathway]
  t_df <- t_df[order(t_df$FP), ]

  t_df[, c("Pathway.ID", "Description", "Total", "Present", "Expected", "FP", "BH")]

}

makeGoNet <- function(entry, globals, isBlastp = FALSE) {
  if (!isBlastp) {
    msg("Making GO network for", entry)
    what <- substr(entry, 1, 1)
  } else {
    msg("Making GO network for blastp entry")
    what <- substr(rownames(entry)[1], 1, 1)
  }
  d <- getTableFiltered(entry, globals, isBlastp = isBlastp)
  go <- makeGOtable(entry, globals, return.go = TRUE, isBlastp = isBlastp)
  if (is.null(go) || is.data.frame(go) || nrow(go[[2]]) == 0) {
    msg("No results")
    return(NULL)
  }
  if (nrow(go[[2]]) > 20) go[[2]] <- go[[2]][1:20, ]
  msg("Found", nrow(go[[2]]), "GO hits")
  go2pf <- genesInTerm(go[[3]], go[[2]][["GO.ID"]])
  go2pf <- lapply(go2pf, function(x) x[x %in% rownames(d)])
  go2pf2 <- lapply(names(go2pf), function(x) data.frame(GO = x, PFAM = go2pf[[x]],
                                                        stringsAsFactors = FALSE))
  edgel <- as.matrix(do.call(rbind, go2pf2))
  edgel <- graph_from_edgelist(edgel)
  data <- toVisNetworkData(edgel)
  go2term <- go[[2]][["Term"]]
  names(go2term) <- go[[2]][["GO.ID"]]
  data$nodes$shape <- "circle"
  data$nodes$color.border <- "black"
  data$nodes$border.width <- 2
  data$nodes$title <- paste(
    data$nodes$id,
    switch(what, P = PFAMDesc[data$nodes$id], T = TIGRFAMDesc[data$nodes$id])
  )
  data$nodes$title[grepl("^GO:", data$nodes$id)] <- paste(
    data$nodes$id[grepl("^GO:", data$nodes$id)],
    go2term[data$nodes$id[grepl("^GO:", data$nodes$id)]]
  )
  data$nodes$color.background <- rep("gold", length(data$nodes$id))
  data$nodes$color.background[grepl("^GO:", data$nodes$id)] <- "tomato"
  visNetwork(nodes = data$nodes, edges = data$edges) %>%
    visNodes(physics = FALSE, borderWidthSelected = 2) %>%
    visInteraction(selectConnectedEdges = FALSE,
                   tooltipDelay = 1) %>%
    visIgraphLayout()
}

makePathwayNet <- function(entry, globals = list(), isBlastp = FALSE) {
  if (!isBlastp) {
    msg("Making pathway network for", entry)
  } else {
    msg("Making pathway network for blastp entry")
  }
  d <- getTableFiltered(entry, globals, isBlastp = isBlastp)
  pthwy <- getKeggEnrich(entry, globals, isBlastp = isBlastp)
  if (nrow(pthwy) == 0) {
    msg("No results")
    return(NULL)
  }
  if (nrow(pthwy) > 20) pthwy <- pthwy[1:20, ]
  msg("Found", nrow(pthwy), "pathway hits")
  pthwys <- pthwy$Pathway.ID
  ppp <- KO2Pathway[KO2Pathway %in% pthwys]
  ppp <- ppp[names(ppp) %in% rownames(d)]
  edgel <- data.frame(KO = names(ppp), Pathway = ppp, row.names = NULL)
  edgel <- as.matrix(edgel)
  edgel <- graph_from_edgelist(edgel)
  data <- toVisNetworkData(edgel)
  data$nodes$shape <- "circle"
  data$nodes$color.border <- "black"
  data$nodes$border.width <- 2
  data$nodes$title <- paste(
    data$nodes$id,
    KODesc[data$nodes$id]
  )
  data$nodes$title[grepl("^map", data$nodes$id)] <- paste(
    data$nodes$id[grepl("^map", data$nodes$id)],
    KEGGPathways[data$nodes$id[grepl("^map", data$nodes$id)]]
  )
  data$nodes$color.background <- rep("gold", length(data$nodes$id))
  data$nodes$color.background[grepl("^map", data$nodes$id)] <- "tomato"
  visNetwork(nodes = data$nodes, edges = data$edges) %>%
    visNodes(physics = FALSE, borderWidthSelected = 2) %>%
    visInteraction(selectConnectedEdges = FALSE,
                   tooltipDelay = 1) %>%
    visIgraphLayout()
}

getRuns <- function(x, y) {
  z <- integer(length(x))
  z[x & y] <- 1
  z[x & !y] <- 2
  z[!x & y] <- 3
  z <- S4Vectors::Rle(z)@values
  x2 <- logical(length(z))
  y2 <- logical(length(z))
  x2[z %in% 1:2] <- TRUE
  y2[z %in% c(1, 3)] <- TRUE
  matrix(c(x2, y2), ncol = 2)
}

JC <- function(x, y) {
  sum(x & y) / sum(x | y)
}

plotMatrixComp <- function(entry, comparison, runs = TRUE, isBlastp = FALSE, bId = "") {
  y <- comparison
  x <- entry
  what <- substr(y, 1, 1)
  if (isBlastp) {
    if (is.null(x)) return()
    x_col <- as.logical(x)
    x <- paste("BLASTP input", bId)
  } else {
    x_col <- as.logical(switch(what,
               P = PFAMTable(x),
               T = TIGRFAMTable(x),
               K = KOTable(x)
             ))
  }
  y_col <- as.logical(switch(what,
             P = PFAMTable(y),
             T = TIGRFAMTable(y),
             K = KOTable(y)
           ))
  msg("Making plotMatrixComp for", x, "and", y)
  out <- if (runs) getRuns(x_col, y_col) else matrix(c(x_col, y_col), ncol = 2)
  treesize <- nrow(out)
  main <- paste0(
    x, " (n = ", sum(out[, 1]), ")", " vs ", y, " (n = ", sum(out[, 2]), ")",
    if (runs) "\nRuns-adjusted tree size = " else "\nNaive tree size = ",
    treesize, if (runs) "\nrJC = " else "\nJC = ",
    round(JC(out[, 1], out[, 2]), 3),
    ", Ov = ", sum(out[, 1] & out[, 2])
  )
  ggplot(melt(out), aes(Var2, Var1)) +
    geom_raster(aes(fill = value)) +
    ggtitle(main) +
    scale_fill_manual(values = c("grey", "forestgreen")) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, vjust = -1, face = "bold"))
}

#------------------------------------------------------------------------------
# Server functions: blastp

checkQuery <- function(x) {
  if (grepl("^>", x))
    return(list(status = FALSE,
                reason = "Please make sure query is not fasta-formatted."))
  if (grepl("[^-ABCDEFGHIKLMNPQRSTUVWYZX*]", x))
    return(list(status = FALSE,
                reason = paste("An unknown character was found. Only the",
                               "following characters (as well as",
                               "whitespace) are allowed:",
                               "ABCDEFGHIKLMNPQRSTUVWYZX*-")))
  if (nchar(x) < 5 || nchar(x) > 10000)
    return(list(status = FALSE,
                reason = paste("Queries must be in between 5-10000 characters",
                               "long.")))
  return(list(status = TRUE))
}

cleanQuery <- function(query) {
  gsub("\\s+", "", toupper(query))
}

checkEmail <- function(x) {
  grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x),
        ignore.case = TRUE)
}

# Make sure to check for unacceptable characters! Also check length.
# According to the blast website, only * and - are accepted beyond the amino
# acid alphabet.
getBlastpId <- function(query, email, evalue, pident, qlenp, slenp, session) {
  msg("Blastp submission")
  query <- cleanQuery(query)
  isOK <- checkQuery(query)
  email <- gsub("\\s+", "", email)
  if (email != "" && !checkEmail(email)) {
    msg("Bad email.")
    return(list(Id = NULL, status = paste("Sorry, but the submitted email",
                                          "address is invalid. Please enter",
                                          "a valid email address or leave",
                                          "the field blank.")))
  }
  if (!isOK$status) {
    msg("Bad query.")
    return(list(Id = NULL, status = isOK$reason))
  }
  # Would also need to make sure all blastp params are exactly the same too..
  # if (file.exists("blastp/queries")) {
  #   queries <- suppressMessages(readr::read_tsv("blastp/queries",
  #                                               col_names = FALSE))
  #   colnames(queries) <- c("JOBID", "EMAIL", "SESSION", "EVALUE", "PIDENT",
  #                          "QLENP", "SLENP", "QUERY")
  #   if (query %in% queries$QUERY) {
  #     msg("Found existing query.")
  #     Id <- queries$JOBID[queries$QUERY == query]
  #     prev_email <- queries$EMAIL[queries$QUERY == query]
  #     # prev_evalue <- queries$EVALUE[queries$]
  #     if (Id %in% strtoi(list.files("blastp/app"))) {
  #       return(list(Id = Id,
  #                   status = paste("Note: a previously completed job contains an",
  #                                  "identical input. Use this job ID to view the",
  #                                  "results.")))
  #     } else if (email == prev_email) {
  #       return(list(Id = Id,
  #                   status = paste("Note: a job currently in progress contains an",
  #                                  "identical input. Use",
  #                                  "this job ID to check on its progress.")))
  #     }
  #   }
  # }
  Id <- as.integer(Sys.time())
  msg("Generated a new job ID:", Id)
  o <- paste0(c(Id, email, session, evalue, pident, qlenp, slenp, query),
              collapse = "\t")
  cat(o, file = "blastp/queries", append = TRUE, sep = "\n")
  list(Id = Id, status = paste("Your job has been submitted. Use this job ID",
                               "to monitor its progress and view results.",
                               "If you also submitted an email address then",
                               "you will be notified once your job is",
                               "complete."))
}

int2date <- function(x) format(as.POSIXct.numeric(x, origin = "1970-01-01"),
                               tz = "America/Toronto")

checkBlastpId <- function(Id) {
  msg("Checking blastp job ID:", Id)
  if (Id %in% strtoi(list.files("blastp/app"))) {
    msg("Found existing app/ file.")
    res <- readRDS(paste0("blastp/app/", Id))
    switch(res$status,
      "OK" = return(list(load = TRUE,
                         status = "Your job has been successfully completed.",
                         res = res)),
      "BLASTP_NOHITS" = return(list(load = FALSE,
                                    status = paste("Unfortunately no",
                                                   "hits were found with the",
                                                   "input blastp parameters.",
                                                   "Try again with more lenient",
                                                   "parameters."))),
      "BLASTP_TOOFEW" = return(list(load = FALSE,
                                    status = paste("Unfortunately too few",
                                                   "hits were found with the",
                                                   "input blastp parameters.",
                                                   "Try again with more lenient",
                                                   "parameters."))),
      "BLASTP_TOOMANY" = return(list(load = FALSE,
                                     status = paste("Unfortunately too many",
                                                    "hits were found for the",
                                                    "analysis to be performed.",
                                                    "Try again with more",
                                                    "stringent parameters."))),
      "BLASTP_ERROR" = return(list(load = FALSE,
                                   status = paste("An error occurred with",
                                                  "blastp. The admin has been",
                                                  "notified.")))
    )
  }
  if (Id %in% strtoi(list.files("blastp/results"))) {
    msg("Found existing results/ file.")
    return(list(load = FALSE,
                status = paste("Your job is currently in progress.")))
                # status = paste("Your job is currently in progress. BLASTP has",
                #                "been successfully run with your query, and",
                #                "phylocorrelation analyses are underway.")))
  }
  if (!file.exists("blastp/queries")) {
    msg("Couldn't find a queries file.")
    return(list(load = FALSE,
                status = "Sorry, no job with this ID exists."))
  }
  msg("Checking queries file.")
  queries <- suppressMessages(readr::read_tsv("blastp/queries",
                                              col_names = FALSE))
  colnames(queries) <- c("JOBID", "EMAIL", "SESSION", "EVALUE", "PIDENT",
                         "QLENP", "SLENP", "QUERY")
  queries <- queries[!queries$JOBID %in% strtoi(list.files("blastp/app")), ]
  if (nrow(queries) == 0 || !Id %in% queries$JOBID) {
    return(list(load = FALSE,
                status = "Sorry, no job with this ID exists."))
  }
  queue <- which(queries$JOBID == Id)
  if (queue == 1) {
    return(list(load = FALSE,
                status = paste("Your job is currently in progress.")))
  }
  if (queue == 2) {
    return(list(load = FALSE,
                status = paste("Your job is in the queue. There is 1 other",
                               "job ahead of it.")))
  }
  list(load = FALSE,
       status = paste("Your job is in the queue. There are", queue - 1,
                      "other jobs ahead of it."))
}

getBlastpInfo <- function(Id) {
  out <- structure(character(9),
                   names = c("Id", "Date1", "Date2", "Evalue", "Pident",
                             "Qlenp", "Slenp", "Query", "qOcc"))
  if (is.null(Id)) return(out)
  if (!file.exists("blastp/queries")) return(out)
  queryF <- suppressMessages(readr::read_tsv("blastp/queries",
                                             col_names = FALSE))
  queryF <- queryF[queryF[[1]] == Id, ]
  if (nrow(queryF) == 0) return(out)
  res <- readRDS(paste0("blastp/app/", Id))
  out["Id"] <- Id
  out["Date1"] <- int2date(Id)
  out["Date2"] <- res$FinishDate
  out["Evalue"] <- queryF[[4]]
  out["Pident"] <- queryF[[5]]
  out["Qlenp"] <- queryF[[6]]
  out["Slenp"] <- queryF[[7]]
  out["Query"] <- queryF[[8]]
  if (!is.null(res$Occ)) out["qOcc"] <- paste(res$Occ, "/ 27372")
  out
}

make_blastp_info_panel <- function(Id, jobInfo) {
  if (is.null(Id)) return("<B>No loaded jobs!</B>")
  paste0(
    "<div\n",
    "  <p><B>ID:</B> ", jobInfo["Id"], "</p>\n",
    "  <p><B>Date submitted:</B> ", jobInfo["Date1"], "</p>\n",
    "  <p><B>Date completed:</B> ", jobInfo["Date2"], "</p>\n",
    "  <p><B>Max E-value:</B> ", jobInfo["Evalue"], "</p>\n",
    "  <p><B>Min percent identity (%):</B> ", jobInfo["Pident"], "</p>\n",
    "  <p><B>Min query alignment (%):</B> ", jobInfo["Qlenp"], "</p>\n",
    "  <p><B>Min subject alignment (%):</B> ", jobInfo["Slenp"], "</p>\n",
    "  <p><B>Number of genomes with hits:</B> ", jobInfo["qOcc"], "</p>\n",
    "</div>"
  )
}

use_blastp_col <- function(x) if (is.null(x)) NA else x

cleanBlastpTable <- function(x) {
  data.frame(
    Links = use_blastp_col(x$Links),
    Description = use_blastp_col(x$Description),
    Compare = use_blastp_col(x$Compare),
    Ov = use_blastp_col(x$Ov),
    Occ = use_blastp_col(x$Occ),
    OccDiff = use_blastp_col(x$OccDiff),
    JC = round(use_blastp_col(x$JC), 3),
    rJC = round(use_blastp_col(x$rJC), 3),
    rHyperP = use_blastp_col(x$rHyperP),
    CS = use_blastp_col(x$CS),
    row.names = rownames(x)
  )
}

#------------------------------------------------------------------------------
# UI functions

p0 <- function(x, y) paste0(x, y)

make_tab_filter_gopath <- function(fam, blastp = FALSE, globals = list()) {

  ff <- function(x) paste0(x, if (blastp) "_BLASTP_" else "_", fam)

  tabPanel(br(), br(),
    if (fam != "KO") {
      radioButtons(
        ff("FILTER_GOPATH_ENRICH_ONTOLOGY"), "Ontology",
        choices = c("BP", "CC", "MF"),
        selected = globals[[ff("FILTER_GOPATH_ENRICH_ONTOLOGY_VALUE")]]
      )
    } else "",
    numericInput(
      ff("FILTER_GOPATH_ENRICH_MAX_FDR"), "Maximum BH",
      globals[[ff("FILTER_GOPATH_ENRICH_MAX_FDR_VALUE")]]
    ),
    numericInput(
      ff("FILTER_GOPATH_ENRICH_MAX_FP"), "Maximum FP",
      globals[[ff("FILTER_GOPATH_ENRICH_MAX_FP_VALUE")]]
    ),
    numericInput(
      ff("FILTER_GOPATH_ENRICH_MIN_TOTAL"), "Minimum total membership",
      globals[[ff("FILTER_GOPATH_ENRICH_MIN_TOTAL_VALUE")]]
    ),
    numericInput(
      ff("FILTER_GOPATH_ENRICH_MIN_PRESENT"), "Minimum present members",
      globals[[ff("FILTER_GOPATH_ENRICH_MIN_PRESENT_VALUE")]]
    ),
    numericInput(
      ff("FILTER_GOPATH_ENRICH_MIN_EXPECTED"), "Minimum expected members",
      globals[[ff("FILTER_GOPATH_ENRICH_MIN_EXPECTED_VALUE")]]
    )
  )

}

make_tab_filter_network <- function(fam, blastp = FALSE, globals = list()) {

  ff <- function(x) paste0(x, if (blastp) "_BLASTP_" else "_", fam)

  tabPanel(br(), br(),
    radioButtons(
      ff("FILTER_CORR_NETWORK_METRIC"), "Metric",
      choices = c("JC", "rJC", "CS"),
      selected = globals[[ff("FILTER_CORR_NETWORK_METRIC_VALUE")]]
    ),
    sliderInput(
      ff("FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE"),
      "Minimum of secondary edges (JC/rJC only)",
      min = 0, max = 1,
      value = globals[[ff("FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE_VALUE")]],
      step = 0.01, ticks = TRUE
    ),
    checkboxGroupInput(
      ff("FILTER_CORR_NETWORK_CS_SECONDARY_EDGE"),
      "Confidence Score of secondary edges (CS only)",
      choices = c("Very low", "Low", "High", "Very high"),
      selected = c("Low", "High", "Very high")
    )
  )

}

make_tab_filter_score_dist <- function(fam, blastp = FALSE, globals = list()) {

  ff <- function(x) paste0(x, if (blastp) "_BLASTP_" else "_", fam)

  tabPanel(br(), br(),
    radioButtons(
      ff("FILTER_SCORE_DIST_METRIC"), "Metric",
      choices = c("JC", "rJC", "rHyperP"),
      selected = globals[[ff("FILTER_SCORE_DIST_METRIC_VALUE")]]
    )
  )

}

make_tab_filter_matrix_comp <- function(fam, blastp = FALSE, globals = list()) {

  ff <- function(x) paste0(x, if (blastp) "_BLASTP_" else "_", fam)

  tabPanel(br(), br(),
    radioButtons(
      ff("FILTER_MATRIX_COMP_TYPE"), "Tree type",
      choices = c("Naive", "Runs-adjusted"),
      selected = globals[[ff("FILTER_MATRIX_COMP_TYPE_VALUE")]]
    )
  )

}

make_tab_filter_empty <- function() {
  tabPanel(br(), br(), "No options available for the current tab.")
}

make_side_panel_input_info <- function(fam, blastp = FALSE, globals = list()) {
  ff <- function(x) paste0(x, if (blastp) "_BLASTP_" else "_", fam)
  req(globals[[ff("INPUT")]])
  x <- switch(fam,
    PFAM = getPFAMinfo(globals[[ff("INPUT")]]),
    TIGRFAM = getTIGRFAMinfo(globals[[ff("INPUT")]]),
    KO = getKEGGinfo(globals[[ff("INPUT")]])
  )
  paste0(
    "<div\n",
    "  <p><B>ID:</B> ", x[[1]], "</p>\n",
    "  <p><B>Description:</B> ", x$description, "</p>\n",
    "  <p><B>Occurrences:</B> ", x$occurrences, " / 27372</p>\n",
    "  <p><B>Other links:</B> ", x$annotree, "</p>\n",
    "</div>"
  )
}

make_side_panel_table <- function(x) {
  if (is.null(x)) return(NULL)
  DT::datatable(x,
    options = list(dom = "pt", bSort = FALSE),
    selection = "none",
    colnames = "",
    class = "cell-border strip hover",
    escape = FALSE
  ) %>% formatStyle(0, cursor = "pointer")
}

check_thresh_msg <- function(THRESH_MSG_ID, keep) {
  if (length(THRESH_MSG_ID) > 0 && !keep) {
    removeNotification(THRESH_MSG_ID)
  }
  THRESH_MSG_ID
}

make_corr_table <- function(entry, fam, globals = list(), THRESH_MSG_ID, isBlastp = FALSE) {
  res <- getTableFiltered(entry, globals, isBlastp = isBlastp)
  if (nrow(res) == 0 && length(THRESH_MSG_ID) == 0) {
      THRESH_MSG_ID <- showNotification(
        "No hits. Try lowering the global thresholds.",
        type = "warning",
        closeButton = TRUE,
        duration = 0
      )
  }
  out <- DT::datatable(
    res,
    extensions = "Buttons",
    options = list(
      pageLength = 25,
      columnDefs = list(
        list(targets = 1:3, bSortable = FALSE),
        list(targets = c(1, 3), className = "dt-center")
      ),
      dom = "Bltip",
      buttons = list(
        list(
          extend = "collection",
          buttons = c("csv", "excel", "pdf"),
          text = "Download"
        )
      )
    ),
    selection = "none",
    container = switch(fam,
      PFAM = PFAMTableCols,
      TIGRFAM = TIGRFAMTableCols,
      KO = KEGGTableCols
    ),
    escape = FALSE
  ) %>% formatStyle(0, cursor = "pointer") %>%
    formatSignif("rHyperP", 3)
  list(table = out, THRESH_MSG_ID = THRESH_MSG_ID, keep = nrow(res) == 0)
}

make_main_go_enrich <- function(globals, fam, isBlastp = FALSE) {
  entry <- if (isBlastp) globals[[fam]] else globals$entry
  DT::datatable(
    switch(fam,
      PFAM = makeGOtable(entry, globals, isBlastp = isBlastp),
      TIGRFAM = makeGOtable(entry, globals, isBlastp = isBlastp),
      KO = getKeggEnrich(entry, globals, isBlastp = isBlastp)
    ),
    extensions = "Buttons",
    options = list(
      pageLength = 10,
      dom = "Bltip",
      buttons = list(
        list(
          extend = "collection",
          buttons = c("csv", "excel", "pdf"),
          text = "Download"
        )
      )
    ),
    container = if (fam == "KO") PathwayTableColumns else GOTableColumns,
    selection = "none"
  ) %>%
    formatSignif("FP", 3) %>%
    formatSignif("BH", 3)
}

make_list_globals <- function(fam, input = list(),
  SidePanelInputValues = list(),
  TabFilterCorrNetworkValues = list(),
  TabFilterGopathEnrichValues = list(),
  TabFilterScoreDistValues = list(),
  TabFilterMatrixCompValues = list(),
  BLASTP_LOADED = list(),
  blastp = FALSE) {

  ff <- function(x) paste0(x, if (blastp) "_BLASTP_" else "_", fam)

  out <- list(
    entry = SidePanelInputValues[[ff("INPUT")]],
    metric = TabFilterCorrNetworkValues[[ff("FILTER_CORR_NETWORK_METRIC_VALUE")]],
    min2nd = TabFilterCorrNetworkValues[[ff("FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE_VALUE")]],
    cswhich = TabFilterCorrNetworkValues[[ff("FILTER_CORR_NETWORK_CS_SECONDARY_EDGE")]],
    Ov = input[[ff("GLOBAL_FILTER_MIN_OV")]],
    Occ = input[[ff("GLOBAL_FILTER_MIN_OCC")]],
    OccDiff = input[[ff("GLOBAL_FILTER_MAX_OCCDIFF")]],
    JC = input[[ff("GLOBAL_FILTER_MIN_JC")]],
    rJC = input[[ff("GLOBAL_FILTER_MIN_RJC")]],
    rHyperP = input[[ff("GLOBAL_FILTER_MAX_RHYPERP")]],
    CS = input[[ff("GLOBAL_FILTER_MIN_PMF")]],
    Ont = TabFilterGopathEnrichValues[[ff("FILTER_GOPATH_ENRICH_ONTOLOGY_VALUE")]],
    FDR = TabFilterGopathEnrichValues[[ff("FILTER_GOPATH_ENRICH_MAX_FDR_VALUE")]],
    FP = TabFilterGopathEnrichValues[[ff("FILTER_GOPATH_ENRICH_MAX_FP_VALUE")]],
    minTotal = TabFilterGopathEnrichValues[[ff("FILTER_GOPATH_ENRICH_MIN_TOTAL_VALUE")]],
    minPresent = TabFilterGopathEnrichValues[[ff("FILTER_GOPATH_ENRICH_MIN_PRESENT_VALUE")]],
    minExp = TabFilterGopathEnrichValues[[ff("FILTER_GOPATH_ENRICH_MIN_EXPECTED_VALUE")]],
    type = TabFilterScoreDistValues[[ff("FILTER_SCORE_DIST_METRIC_VALUE")]],
    runs = TabFilterMatrixCompValues[[ff("FILTER_MATRIX_COMP_TYPE_VALUE")]],
    comparison = TabFilterMatrixCompValues[[ff("FILTER_MATRIX_COMP_COMPARISON_VALUE")]],
    PFAM = BLASTP_LOADED$Dat$PFAM,
    TIGRFAM = BLASTP_LOADED$Dat$TIGRFAM,
    KO = BLASTP_LOADED$Dat$KO,
    Vec = BLASTP_LOADED$Dat$Vec,
    Id = BLASTP_LOADED$Id
  )

  if (!is.null(out$runs)) out$runs <- out$runs == "Runs-adjusted"
  if (!is.null(out$rHyperP)) out$rHyperP <- 10^(-out$rHyperP)

  out

}

make_tab_main <- function(fam, blastp = FALSE) {

  ff <- function(x) paste0(x, if (blastp) "_BLASTP_" else "_", fam)

  tagList(

    br(),

    sidebarPanel(
      tabsetPanel(
        tabPanel("Input",
          br(),
          htmlOutput(ff("SIDE_PANEL_INPUT_INFO")),
          br(),
          if (!blastp) searchInput(
            ff("SIDE_PANEL_INPUT_SEARCH"),
            paste0(fam, " search"),
            btnSearch = icon("search"),
            btnReset = icon("remove")
          ),
          if (!blastp) textOutput(ff("SIDE_PANEL_INPUT_SEARCH_TEXT")),
          if (!blastp) DT::dataTableOutput(ff("SIDE_PANEL_INPUT_SEARCH_TABLE"))
        ),
        tabPanel("Global filters",
          br(),
          sliderInput(
            ff("GLOBAL_FILTER_MIN_OV"), "Minimum Ov",
            min = 0, max = 27372, value = 2, step = 1, ticks = TRUE
          ),
          sliderInput(
            ff("GLOBAL_FILTER_MIN_OCC"), "Minimum Occ",
            min = 0, max = 27372, value = 2, step = 1, ticks = TRUE
          ),
          sliderInput(
            ff("GLOBAL_FILTER_MAX_OCCDIFF"), "Maximum OccDiff",
            min = 0, max = 27372, value = 27372, step = 1, ticks = TRUE
          ),
          sliderInput(
            ff("GLOBAL_FILTER_MIN_JC"), "Minimum JC",
            min = 0, max = 1, value = 0, step = 0.01, ticks = TRUE
          ),
          sliderInput(
            ff("GLOBAL_FILTER_MIN_RJC"), "Minimum rJC",
            min = 0, max = 1, value = 0, step = 0.01, ticks = TRUE
          ),
          sliderInput(
            ff("GLOBAL_FILTER_MAX_RHYPERP"),
            "Maximum -log10(rHyperP)",
            min = 0, max = 320, value = 0, step = 1, ticks = TRUE
          ),
          checkboxGroupInput(
            ff("GLOBAL_FILTER_MIN_PMF"), "Confidence Score",
            choices = c("Very low", "Low", "High", "Very high"),
            selected = c("Low", "High", "Very high")
          )
        ),
        tabPanel("Tab Filters",
          uiOutput(ff("UI_SIDE_TAB_FILTERS"))
        )
      )
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Correlation Table",
          br(),
          if (!blastp) tags$em("Click on family IDs to use them as query."),
          br(), br(),
          DT::dataTableOutput(ff("MAIN_CORR_TABLE")),
          value = ff("TAB_MAIN_CORR_TABLE")
        ),
        tabPanel("Correlation Network",
          br(),
          tags$em("Max 200 nodes. Right-click to save image."),
          br(),
          visNetworkOutput(
            ff("MAIN_CORR_NETWORK"),
            height = "800px"
          ) %>% withSpinner(),
          value = ff("TAB_MAIN_CORR_NETWORK")
        ),
        tabPanel(switch(fam, KO = "Pathway Enrichment", "GO Enrichment"),
          br(),
          DT::dataTableOutput(ff("MAIN_GOPATH_ENRICH")) %>%
            withSpinner(),
          value = ff("TAB_MAIN_GOPATH_ENRICH")
        ),
        tabPanel("Score Distribution",
          br(),
          plotlyOutput(ff("MAIN_SCORE_DIST"), height = "auto") %>%
            withSpinner(),
          value = ff("TAB_MAIN_SCORE_DIST")
        ),
        tabPanel("Matrix Comparison",
          br(),
          tags$em("Right-click to save iamge."),
          br(),
          plotOutput(ff("MAIN_MATRIX_COMP"), height = "1000px") %>%
            withSpinner(),
          value = ff("TAB_MAIN_MATRIX_COMP")
        ),
        id = ff("MAIN_TABS")
      )
    )

  )

}

make_tab_blastp_input <- function() tagList(
  tabPanel("Input & jobs",
    tags$h3("Job submission"),
    br(),
    sidebarLayout(
      sidebarPanel(
        tags$h4("BLASTP parameters"),
        br(),
        numericInput("BLASTP_FILTER_EVAL", "Max E-value", 0.00001),
        numericInput("BLASTP_FILTER_PIDENT", "Min percent identity (%)", 30),
        numericInput("BLASTP_FILTER_QLENP", "Min query alignment (%)", 70),
        numericInput("BLASTP_FILTER_SLENP", "Min subject alignment (%)", 70)
      ),
      mainPanel(
        textAreaInput("BLASTP_INPUT_SEQUENCE",
          "Input a single amino acid sequence string:",
          value = "", placeholder = "", width = "600px", height = "200px"
        ),
        textInput("BLASTP_INPUT_EMAIL", "Email (recommended!):"),
        actionButton("BLASTP_INPUT_BUTTON", "Submit")
      )
    ),
    hr(),
    tags$h3("Job retrieval"),
    br(),
    sidebarLayout(
      sidebarPanel(
        tags$h4("Loaded job info"),
        br(),
        htmlOutput("BLASTP_JOB_INFO")
      ),
      mainPanel(
        textInput("BLASTP_CHECK_JOB_ID", "Enter your job ID:", ""),
        actionButton("BLASTP_CHECK_JOB_BUTTON", "Submit")
      )
    )
  , value = "BLASTP_TABS_INPUT")
)

make_tab_help <- function() tagList(
  tags$h2("FAQ"),

  tags$h3("1. I can't find my gene of interest."),
  tags$p("First, check whether your gene/protein is covered by an existing entry within the KEGG ,TIGRFAM, or PFAM database. In many cases, an entry may exist but may be labeled or described differently."),
  tags$p("If your gene does not exist, you may submit your (protein) sequence of interest and PhyloCorrelate will:"),
  tags$ol(
    tags$li("BLAST the sequence against the AnnoTree database"),
    tags$li("Identify all matches across the bacterial tree"),
    tags$li("Perform a co-occurrence analysis")
  ),

  tags$h3("2. PhyloCorrelate didn't find any correlations for my gene."),
  tags$p("Try exploring different parameters!"),
  tags$h4("a) Reduce the -logP threshold"),
  tags$p("Co-occurring genes may be present below the default threshold, especially if your query gene occurs infrequently across the tree (e.g., <100 occurrences). In these cases, it is statistically harder for these patterns to achieve a significant p-value."),
  tags$h4("b) Reduce the OccDiff parameter"),
  tags$p("This is the different in the number of occurrences of your query gene and potential correlating genes across the tree. If you increase this value, you may start detecting dependency-type relationships where your query gene depends on the existence of other genes for its function.")
)

#------------------------------------------------------------------------------
# Table column name hover tips

PFAMTableCols <- htmltools::withTags(table(
  class = "display",
  thead(
    tr(
      th("", title = "PFAM ID"),
      th("Links", title = "[P] pfam.xam.org [A] annotree.uwaterloo.ca"),
      th("Description", title = "PFAM Description"),
      th("Compare", title = "Load this comparison in the Matrix Comparison tab"),
      th("Ov", title = "Overlapping Occurrences"),
      th("Occ", title = "Total Occurrences"),
      th("OccDiff", title = "Difference of Occurrences"),
      th("JC", title = "Jaccard Coefficient"),
      th("rJC", title = "Runs-adjusted Jaccard Coefficient"),
      th("rHyperP", title = "Hypergeometric P-value"),
      th("CS",  title = "Confidence Score: Probability of a Matching GO:BP Term")
    )
  )
))

TIGRFAMTableCols <- htmltools::withTags(table(
  class = "display",
  thead(
    tr(
      th("", title = "TIGRFAM ID"),
      th("Links", title = "[T] tigrfams.jcvi.org [A] annotree.uwaterloo.ca"),
      th("Description", title = "TIGRFAM Description"),
      th("Compare", title = "Load this comparison in the Matrix Comparison tab"),
      th("Ov", title = "Overlapping Occurrences"),
      th("Occ", title = "Total Occurrences"),
      th("OccDiff", title = "Difference of Occurrences"),
      th("JC", title = "Jaccard Coefficient"),
      th("rJC", title = "Runs-adjusted Jaccard Coefficient"),
      th("rHyperP", title = "Hypergeometric P-value"),
      th("CS",  title = "Confidence Score: Probability of a Matching GO:BP Term")
    )
  )
))

KEGGTableCols <- htmltools::withTags(table(
  class = "display",
  thead(
    tr(
      th("", title = "KEGG Ortholog ID"),
      th("Links", title = "[K] genome.jp [A] annotree.uwaterloo.ca"),
      th("Description", title = "KEGG Ortholog Description"),
      th("Compare", title = "Load this comparison in the Matrix Comparison tab"),
      th("Ov", title = "Overlapping Occurrences"),
      th("Occ", title = "Total Occurrences"),
      th("OccDiff", title = "Difference of Occurrences"),
      th("JC", title = "Jaccard Coefficient"),
      th("rJC", title = "Runs-adjusted Jaccard Coefficient"),
      th("rHyperP", title = "Hypergeometric P-value"),
      th("CS",  title = "Confidence Score: Probability of a Matching Pathway")
    )
  )
))

GOTableColumns <- htmltools::withTags(table(
  class = "display",
  thead(
    tr(
      th("", title = ""),
      th("GO.ID", title = "Gene Ontology"),
      th("Term", title = "GO term description"),
      th("Total", title = "GO term group size"),
      th("Present", title = "Members with this term"),
      th("Expected", title = "Members expected to have this term"),
      th("FP", title = "Fisher's Exact Test P-value"),
      th("BH", title = "Benjamini & Hochberg Adjusted P-value")
    )
  )
))

PathwayTableColumns <- htmltools::withTags(table(
  class = "display",
  thead(
    tr(
      th("", title = ""),
      th("Pathway.ID", title = "Pathway ID"),
      th("Description", title = "Pathway description"),
      th("Total", title = "Pathway group size"),
      th("Present", title = "Members with this term"),
      th("Expected", title = "Members expected to have this term"),
      th("FP", title = "Fisher's Exact Test P-value"),
      th("BH", title = "Benjamini & Hochberg Adjusted P-value")
    )
  )
))

msg("Finished loading global.R")
