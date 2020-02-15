# PhyloCorrelate: blastp script.
# Benjamin Jean-Marie Tremblay
# 2019-09-17
# Last modified: 2020-02-14 22:42:09 (CET)
# This script is meant to run in the background, as a separate
# process from the actual app. Whenever a job is submitted, a new entry is
# added to the "queries" file. This is detected by the script and an analysis
# is launched. Only a single job is run at a time. The script quits whenever
# it runs out of jobs. The main app will run the script whenever it boots
# up, gets a new job, or checks on a job. If the script is already running
# the new instance quits, only ever leaving one instance running at a time.

# Note: do not ever be tempted to leave this script running continuously. R
# will always eventually crash. It may takes weeks, but it will happen.

## Script parameters ----------------------------------------------------------

BlastpVer <- 3L

AppDir <- "PhyloCorrApp3"
AppURL <- "https://phylocorrelate.uwaterloo.ca/app/"
Proteomes <- "~/PhyloCorrProteomes/AllProteomes.faa"

admin_email <- "benjmtremblay@gmail.com"

queries_columns <- c("JOBID", "EMAIL", "SESSION", "EVALUE",
                     "PIDENT", "QLENP", "SLENP", "QUERY")
BlastpSleepTime <- 5

setwd("~")

## ----------------------------------------------------------------------------

Pdir <- function(x = c()) paste0(c(AppDir, x), collapse = "/")
Ddir <- function(x = c()) Pdir(paste0(c("data", x), collapse = "/"))
Bdir <- function(x = c()) Pdir(paste0(c("blastp", x), collapse = "/"))
Mdir <- function(x = c()) Bdir(paste0(c("misc", x), collapse = "/"))
Adir <- function(x = c()) Bdir(paste0(c("app", x), collapse = "/"))
Rdir <- function(x = c()) Bdir(paste0(c("results", x), collapse = "/"))

checkDirs <- function() {
  if (!dir.exists(Pdir())) {
    print("Early exit: bad dir structure")
    q()
  }
}
checkDirs()

Date <- function() format(as.POSIXlt(Sys.time(), tz = "America/Toronto"))
msg <- function(...) {
  d <- Date()
  x <- paste0(c(d, ...), collapse = " ")
  cat(x, sep = "\n", file = Bdir("logs"), append  = TRUE)
}

pid <- Sys.getpid()

err_handler <- function() {
  msg("Error occurred!")
  cat(geterrmessage(), file = Bdir("logs"), append = TRUE)
  file.remove(Bdir("LIVE"))
  gmailr::gm_send_message(
    gmailr::gm_mime() %>%
      gmailr::gm_from("phylocorrelate@gmail.com") %>%
      gmailr::gm_to(admin_email) %>%
      gmailr::gm_subject("PhyloCorrelate BLASTP error") %>%
      gmailr::gm_text_body(geterrmessage())
  )
  msg("Quitting with pid", pid)
  q()
}

options(error = err_handler)

msg("Booting up with pid", pid)

dir.create(Bdir("results"), showWarnings = FALSE)
dir.create(Bdir("app"), showWarnings = FALSE)
if (file.exists(Bdir("KILL"))) {
  msg("Found KILL file, quitting", pid)
  q()
}

getOldPid <- function() {
  con <- file(Bdir("LIVE"))
  pid0 <- as.integer(suppressWarnings(readLines(con)))
  close(con)
  pid0
}

if (file.exists(Bdir("LIVE"))) {
  pid0 <- getOldPid()
  msg("Detected existing LIVE file, checking pid", pid0)
  status = try(ps::ps_status(ps::ps_handle(pid0)), silent = TRUE)
  if (!inherits(status, "try-error")) {
    msg("pid", pid0, "is running, quitting pid", pid)
    q()
  } else {
    msg("pid", pid0, "couldn't be found, running")
  }
}

cat(pid, file = Bdir("LIVE"))

suppressPackageStartupMessages(library(fst))

genomes <- readLines(con <- file(Ddir("genomes.txt"))) ; close(con)

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
  list(x = x2, y = y2)
}
JC <- function(x, y) {
  sum(x & y) / sum(x | y)
}
Ov <- function(x, y) {
  sum(x & y)
}
tdr_pfam <- function(x) {
  mod <- readRDS(Bdir("fdr/PFAM_TDR_approxfun.RDS"))
  mod(x)
}
tdr_pfam_rhp <- function(x) {
  mod <- readRDS(Bdir("fdr/PFAM_TDR_approxfun_rHyperP.RDS"))
  mod(x)
}
tdr_tigrfam <- function(x) {
  mod <- readRDS(Bdir("fdr/TIGRFAM_TDR_approxfun.RDS"))
  mod(x)
}
tdr_tigrfam_rhp <- function(x) {
  mod <- readRDS(Bdir("fdr/TIGRFAM_TDR_approxfun_rHyperP.RDS"))
  mod(x)
}
tdr_ko <- function(x) {
  mod <- readRDS(Bdir("fdr/KO_TDR_approxfun.RDS"))
  mod(x)
}
tdr_ko_rhp <- function(x) {
  mod <- readRDS(Bdir("fdr/KO_TDR_approxfun_rHyperP.RDS"))
  mod(x)
}
testHyper <- function(z, y) {
  fisher.test(
    matrix(c(
      sum(z & y), sum(y & !z),
      sum(!y & z), sum(!y & !z)
    ), ncol = 2),
    alternative = "g"
  )$p.value
}

msg("Running")

gmailr::gm_auth_configure(Bdir("credentials.json"))
gmailr::gm_auth(email = TRUE, cache = Bdir(".secret"))

send_mail <- function(Id, email) {
  gmailr::gm_send_message(
    gmailr::gm_mime() %>%
      gmailr::gm_from("phylocorrelate@gmail.com") %>%
      gmailr::gm_to(email) %>%
      gmailr::gm_subject(paste("PhyloCorrelate BLASTP job", Id))  %>%
      gmailr::gm_text_body(paste(
        "Your job is done! Visit", paste0(AppURL, "?blastp=", Id),
        "to view your results."
      ))
  )
}

send_mail_msg <- function(Id, email, msg) {
  gmailr::gm_send_message(
    gmailr::gm_mime() %>%
      gmailr::gm_from("phylocorrelate@gmail.com") %>%
      gmailr::gm_to(email) %>%
      gmailr::gm_subject(paste("PhyloCorrelate BLASTP job", Id)) %>%
      gmailr::gm_text_body(msg)
  )
}

int2date <- function(x) as.Date(as.POSIXct.numeric(x, origin = "1970-01-01"))

badBlast <- function(msg, JOBID, oA) {
    msg("Bad blast.")
    saveRDS(list(status = msg, Id = JOBID, version = BlastpVer), oA,
            compress = FALSE)
}

## ----------------------------------------------------------------------------

repeat {

  Sys.sleep(BlastpSleepTime)

  checkDirs()
  dir.create(Bdir(), FALSE)
  dir.create(Mdir(), FALSE)

  if (!file.exists(Bdir("queries"))) {
    msg("No new queries, BREAK")
    break
  }
  if (file.exists(Bdir("KILL"))) break
  if (!file.exists(Bdir("LIVE"))) {
    msg("LIVE file is missing, creating new file")
    cat(pid, file = Bdir("LIVE"))
  }
  pid0 <- getOldPid()
  if (pid0 != pid) {
    msg("Pid in LIVE file doesn't match current process, quitting")
    send_message(
      mime(
        To = admin_email,
        From = "phylocorrelate@gmail.com",
        Subject = "PhyloCorrApp/blastp",
        body = "malformed pid in LIVE file"
      )
    )
    break 
  }

  queries <- suppressMessages(readr::read_tsv(Bdir("queries"), col_names = FALSE))

  if (nrow(queries) == 0) {
    msg("No new queries, BREAK")
    break
  }
  if (ncol(queries) != length(queries_columns)) {
    msg("Detected badly formatted queries, moving to misc")
    send_message(
      mime(
        To = admin_email,
        From = "phylocorrelate@gmail.com",
        Subject = "PhyloCorrApp/blastp",
        body = "malformed queries file"
      )
    )
    file.copy(from = Bdir("queries"), to = Mdir(paste(Date(), "queries")))
    file.remove(Bdir("queries"))
  }

  colnames(queries) <- queries_columns
  queries <- queries[!queries$JOBID %in% strtoi(list.files(Adir())), ]
  if (nrow(queries) == 0) {
    msg("No new queries, BREAK")
    break
  }

  msg("Running job", queries$JOBID[1])

  f <- Bdir("query.faa")
  cat(c(">1", queries$QUERY[1]), file = f, sep = "\n", append = FALSE)

  o <- Rdir(queries$JOBID[1])
  oA <- Adir(queries$JOBID[1])
  qLen <- nchar(queries$QUERY[1])
  cmd <- paste("blastp -query", f, "-db ", Proteomes,
               "-out", o, "-outfmt '6 qseqid sseqid evalue pident length'",
               "-evalue", queries$EVALUE[1],
               "-max_target_seqs 100000000 -num_threads 2")

  if (queries$JOBID[1] %in% strtoi(list.files(Rdir()))) {
    msg("Detected previous blastp results")
  } else {
    msg("Running blastp")
    b <- system(cmd)
    if (b != 0) {
      send_message(
        mime(
          To = admin_email,
          From = "phylocorrelate@gmail.com",
          Subject = "PhyloCorrApp/blastp",
          body = paste0("failed blastp run (JOBID#", queries$JOBID, ")")
        )
      )
      badBlast("BLASTP_ERROR", queries$JOBID[1], oA)
      if (queries$EMAIL[1] != "")
        send_mail_msg(queries$JOBID[1], queries$EMAIL[1],
                      paste("Sorry, but an error occurred while running blastp.",
                            "The admin has been notified and a follow-up",
                            "email will be sent once the issue has been resolved."))
      next
    }
    msg("blastp run successful")
  }

  res <- suppressMessages(readr::read_tsv(o, col_names = FALSE))
  if (nrow(res) == 0) {
    badBlast("BLASTP_NOHITS", queries$JOBID[1], oA)
    next
  }

  colnames(res) <- c("qseqid", "sseqid", "evalue", "pident", "length")
  res <- res[res$pident >= queries$PIDENT[1], ]
  if (nrow(res) == 0) {
    badBlast("BLASTP_NOHITS", queries$JOBID[1], oA)
    next
  }

  res$qlenp <- res$length / qLen * 100
  rm(qLen)
  res <- res[res$qlenp >= queries$QLENP[1], ]
  if (nrow(res) == 0) {
    badBlast("BLASTP_NOHITS", queries$JOBID[1], oA)
    next
  }

  res$slenp <- vapply(strsplit(res$sseqid, "_L", fixed = TRUE),
                      function(x) rev(x)[1], character(1))
  res$slenp <- as.numeric(res$slenp)
  res$slenp <- res$length / res$slenp * 100
  res <- res[res$slenp >= queries$SLENP[1], ]
  if (nrow(res) == 0) {
    badBlast("BLASTP_NOHITS", queries$JOBID[1], oA)
    next
  }

  resCounts <- logical(length(genomes))
  names(resCounts) <- genomes
  for (i in seq_along(resCounts)) {
    if (any(grepl(names(resCounts)[i], res$sseqid))) resCounts[i] <- TRUE
  }

  if (sum(resCounts) <= 1) {
    badBlast("BLASTP_TOOFEW", queries$JOBID[1], oA)
    next
  }

  if (all(resCounts)) {
    badBlast("BLASTP_TOOMANY", queries$JOBID[1], oA)
    next
  }

  msg("Running analyses")

  msg("Comparing to PFAMs ..")

  PFAMTable <- read_fst(Ddir("PFAMTable.fst"))

  PFAM_Ov <- vapply(PFAMTable, function(x) Ov(resCounts, x),
                    integer(1))
  PFAM_rJC <- numeric(ncol(PFAMTable))
  PFAM_JC <- numeric(ncol(PFAMTable))
  PFAM_Hyper <- numeric(ncol(PFAMTable))
  for (i in seq_len(ncol(PFAMTable))) {
    PFAM_Hyper[i] <- testHyper(resCounts, PFAMTable[, i])
    PFAM_JC[i] <- JC(resCounts, PFAMTable[, i])
    PFAM_runs <- getRuns(resCounts, PFAMTable[, i])
    PFAM_rJC[i] <- JC(PFAM_runs$x, PFAM_runs$y)
  }
  PFAM_tdr <- tdr_pfam(PFAM_rJC)
  PFAM_tdr_rhp <- tdr_pfam_rhp(PFAM_Hyper)

  PFAMs <- colnames(PFAMTable)

  rm(PFAM_runs)
  rm(PFAMTable)

  msg("OK")
  msg("Comparing to TIGRFAMs ..")

  TIGRFAMTable <- read_fst(Ddir("TIGRFAMTable.fst"))

  TIGRFAM_Ov <- vapply(TIGRFAMTable, function(x) Ov(resCounts, x),
                       integer(1))
  TIGRFAM_rJC <- numeric(ncol(TIGRFAMTable))
  TIGRFAM_JC <- numeric(ncol(TIGRFAMTable))
  TIGRFAM_Hyper <- numeric(ncol(TIGRFAMTable))
  for (i in seq_len(ncol(TIGRFAMTable))) {
    TIGRFAM_Hyper[i] <- testHyper(resCounts, TIGRFAMTable[, i])
    TIGRFAM_JC[i] <- JC(resCounts, TIGRFAMTable[, i])
    TIGRFAM_runs <- getRuns(resCounts, TIGRFAMTable[, i])
    TIGRFAM_rJC[i] <- JC(TIGRFAM_runs$x, TIGRFAM_runs$y)
  }
  TIGRFAM_tdr <- tdr_tigrfam(TIGRFAM_rJC)
  TIGRFAM_tdr_rhp <- tdr_tigrfam_rhp(TIGRFAM_Hyper)

  TIGRFAMs <- colnames(TIGRFAMTable)

  rm(TIGRFAM_runs)
  rm(TIGRFAMTable)

  msg("OK")
  msg("Comparing to KOs ..")

  KOTable <- read_fst(Ddir("KOTable.fst"))

  KO_Ov <- vapply(KOTable, function(x) Ov(resCounts, x), integer(1))
  KO_rJC <- numeric(ncol(KOTable))
  KO_JC <- numeric(ncol(KOTable))
  KO_Hyper <- numeric(ncol(KOTable))
  for (i in seq_len(ncol(KOTable))) {
    KO_Hyper[i] <- testHyper(resCounts, KOTable[, i])
    KO_JC[i] <- JC(resCounts, KOTable[, i])
    KO_runs <- getRuns(resCounts, KOTable[, i])
    KO_rJC[i] <- JC(KO_runs$x, KO_runs$y)
  }
  KO_tdr <- tdr_ko(KO_rJC)
  KO_tdr_rhp <- tdr_ko_rhp(KO_Hyper)

  KOs <- colnames(KOTable)

  rm(KO_runs)
  rm(KOTable)

  msg("OK")
  msg("Prepping output ..")

  PFAMCounts <- readRDS(Ddir("PFAMCounts.RDS"))
  PFAMLinks <- readRDS(Ddir("PFAMLinks2.RDS"))
  PFAMDesc <- readRDS(Ddir("PFAMDescriptions.RDS"))
  PFAMButtons <- readRDS(Ddir("PFAMButtons2.RDS"))

  msg("  Making PFAM results")
  
  msg(class(PFAMLinks))
  msg(class(PFAMDesc))
  msg(class(PFAMButtons))
  msg(class(PFAM_Ov))
  msg(class(PFAMCounts))
  msg(class(PFAMCounts - sum(resCounts)))
  msg(class(PFAM_JC))
  msg(class(PFAM_rJC))
  msg(class(PFAM_Hyper))
  msg(class(PFAM_tdr))
  msg(class(PFAM_tdr_rhp))
  msg(class(PFAMs))

  PFAM_all <- data.frame(
    Links = PFAMLinks[PFAMs],
    Description = PFAMDesc[PFAMs],
    Compare = PFAMButtons[PFAMs],
    Ov = PFAM_Ov,
    Occ = PFAMCounts,
    OccDiff = abs(PFAMCounts - sum(resCounts)),
    JC = PFAM_JC,
    rJC = PFAM_rJC,
    rHyperP = PFAM_Hyper,
    PMF = PFAM_tdr,
    PMF2 = PFAM_tdr_rhp,
    row.names = PFAMs,
    stringsAsFactors = FALSE
  )
  PFAM_all <- PFAM_all[order(PFAM_all$rJC, decreasing = TRUE), ]

  rm(PFAM_Ov)
  rm(PFAMButtons)
  rm(PFAM_rJC)
  rm(PFAM_JC)
  rm(PFAM_Hyper)
  rm(PFAM_tdr)
  rm(PFAM_tdr_rhp)
  rm(PFAMCounts)
  rm(PFAMLinks)
  rm(PFAMDesc)
  rm(PFAMs)

  msg("  Making TIGRFAM results")

  TIGRFAMCounts <- readRDS(Ddir("TIGRFAMCounts.RDS"))
  TIGRFAMLinks <- readRDS(Ddir("TIGRFAMLinks2.RDS"))
  TIGRFAMDesc <- readRDS(Ddir("TIGRFAMDescriptions.RDS"))
  TIGRFAMButtons <- readRDS(Ddir("TIGRFAMButtons2.RDS"))

  TIGRFAM_all <- data.frame(
    Links = TIGRFAMLinks[TIGRFAMs],
    Description = TIGRFAMDesc[TIGRFAMs],
    Compare = TIGRFAMButtons[TIGRFAMs],
    Ov = TIGRFAM_Ov,
    Occ = TIGRFAMCounts,
    OccDiff = abs(TIGRFAMCounts - sum(resCounts)),
    JC = TIGRFAM_JC,
    rJC = TIGRFAM_rJC,
    rHyperP = TIGRFAM_Hyper,
    PMF = TIGRFAM_tdr,
    PMF2 = TIGRFAM_tdr_rhp,
    row.names = TIGRFAMs,
    stringsAsFactors = FALSE
  )
  TIGRFAM_all <- TIGRFAM_all[order(TIGRFAM_all$rJC, decreasing = TRUE), ]

  rm(TIGRFAM_Ov)
  rm(TIGRFAMButtons)
  rm(TIGRFAM_rJC)
  rm(TIGRFAM_JC)
  rm(TIGRFAM_tdr)
  rm(TIGRFAM_tdr_rhp)
  rm(TIGRFAM_Hyper)
  rm(TIGRFAMCounts)
  rm(TIGRFAMLinks)
  rm(TIGRFAMDesc)
  rm(TIGRFAMs)

  msg("  Making KO results")

  KOCounts <- readRDS(Ddir("KOCounts.RDS"))
  KOLinks <- readRDS(Ddir("KOLinks2.RDS"))
  KODesc <- readRDS(Ddir("KODescriptions.RDS"))
  KOButtons <- readRDS(Ddir("KOButtons2.RDS"))

  KO_all <- data.frame(
    Links = KOLinks[KOs],
    Description = KODesc[KOs],
    Compare = KOButtons[KOs],
    Ov = KO_Ov,
    Occ = KOCounts,
    OccDiff = abs(KOCounts - sum(resCounts)),
    JC = KO_JC,
    rJC = KO_rJC,
    rHyperP = KO_Hyper,
    PMF = KO_tdr,
    PMF2 = KO_tdr_rhp,
    row.names = KOs,
    stringsAsFactors = FALSE
  )
  KO_all <- KO_all[order(KO_all$rJC, decreasing = TRUE), ]

  rm(KO_Ov)
  rm(KOButtons)
  rm(KO_rJC)
  rm(KO_JC)
  rm(KO_Hyper)
  rm(KO_tdr)
  rm(KO_tdr_rhp)
  rm(KOCounts)
  rm(KOLinks)
  rm(KODesc)
  rm(KOs)

  msg("  Combining results")

  finalOut <- list(
    status = "OK",
    version = BlastpVer,
    FinishDate = Date(),
    Vec = resCounts,
    Occ = sum(resCounts),
    Id = queries$JOBID[1],
    PFAM = PFAM_all,
    TIGRFAM = TIGRFAM_all,
    KO = KO_all
  )

  msg("  Saving..")

  saveRDS(finalOut, oA, compress = FALSE)

  rm(PFAM_all)
  rm(TIGRFAM_all)
  rm(KO_all)

  rm(res)
  rm(resCounts)
  rm(finalOut)
  rm(o)
  rm(oA)
  rm(cmd)
  rm(f)

  if (queries$EMAIL[1] != "") send_mail(queries$JOBID[1], queries$EMAIL[1])

  msg("Job done for ID", queries$JOBID[1])

  rm(queries)

}

## ----------------------------------------------------------------------------

msg("Detected KILL order, quitting with pid", pid)
if (file.exists(Bdir("LIVE"))) file.remove(Bdir("LIVE"))
q()
