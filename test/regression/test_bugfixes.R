#!/usr/bin/env Rscript
# Regression tests for the three bugs reported 2026-11-17 (Sara Trigueros).
# Run:  Rscript test/regression/test_bugfixes.R
# Each test must FAIL on the unfixed code and PASS after the fix.

suppressPackageStartupMessages({library(DESeq2); library(topGO)})

# the fixes under test live here
args <- commandArgs(trailingOnly = FALSE)
f    <- sub("^--file=", "", args[grep("^--file=", args)])
repo <- if (length(f)) dirname(dirname(dirname(normalizePath(f)))) else "."
source(file.path(repo, "scripts", "NyaRNAthoDEG_helpers.R"))
fails <- 0
ok <- function(name, cond, info="") {
  cat(sprintf("%-58s %s%s\n", name, if (isTRUE(cond)) "PASS" else "FAIL",
              if (nzchar(info)) paste0("  (", info, ")") else ""))
  if (!isTRUE(cond)) fails <<- fails + 1
}

# ── 1. VST must work on a small dataset (< 1000 genes) ───────────────────────
set.seed(1)
n <- 500
cts <- matrix(rnbinom(n * 6, mu = 100, size = 1/0.3), ncol = 6,
              dimnames = list(paste0("g", 1:n), paste0("s", 1:6)))
dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(
         cts, data.frame(cond = factor(rep(c("c","t"), each = 3))), ~cond), quiet = TRUE))

r <- try(.vst_safe(dds), silent = TRUE)
ok("vst on 500 genes does not error", !inherits(r, "try-error"))

# ── 2. GenTable's statistic column must be usable at topNodes = 1 ────────────
genes <- paste0("G", sprintf("%04d", 1:300))
gos   <- c("GO:0006355","GO:0005524","GO:0003723","GO:0009414","GO:0005737")
set.seed(42)
geneID2GO <- setNames(lapply(seq_along(genes),
                   function(i) sample(gos, sample(1:3, 1))), genes)
gv <- factor(as.integer(genes %in% sample(genes, 30))); names(gv) <- genes
GOdata <- suppressMessages(new("topGOdata", ontology = "BP", allGenes = gv,
                               annot = annFUN.gene2GO, gene2GO = geneID2GO))
res <- suppressMessages(runTest(GOdata, statistic = "fisher"))

tbl1 <- .gentable_fisher(GOdata, res, 1)
ok("GenTable(topNodes=1) exposes a usable 'Fishers' column",
   !is.null(tbl1$Fishers) && length(suppressWarnings(as.numeric(tbl1$Fishers))) == nrow(tbl1),
   sprintf("nrow=%d", nrow(tbl1)))
r2 <- try({tbl1$log10F <- -log10(suppressWarnings(as.numeric(tbl1$Fishers))); TRUE}, silent = TRUE)
ok("log10F assignment on a 1-row table does not error", !inherits(r2, "try-error"))

# ── 3. The AHRD fixture + tolerant reader ───────────────────────────────────
ahrd_path <- file.path(repo, "test", "test_ahrd.txt")
fan_path  <- file.path(repo, "test", "test_fantasia.txt")
first <- readLines(ahrd_path, n = 1, warn = FALSE)
ok("test_ahrd.txt starts with the '# AHRD-Version' comment", grepl("^#", first), first)
a <- .read_ahrd(ahrd_path)
f <- read.delim(fan_path, header = FALSE)
ok("AHRD and FANTASIA fixtures expose the same gene count",
   nrow(a) == nrow(f), sprintf("ahrd=%d fantasia=%d", nrow(a), nrow(f)))
ok("AHRD header row is a header, not a gene",
   !grepl("^SMELv", colnames(a)[1]), colnames(a)[1])

# the tolerant reader must also cope with a file that has no preamble
tmp <- tempfile(fileext = ".txt")
writeLines(readLines(ahrd_path)[-(1:2)], tmp)
b <- suppressWarnings(.read_ahrd(tmp))
ok("reader handles a hand-built file with no preamble",
   nrow(b) == nrow(a), sprintf("no-preamble=%d with-preamble=%d", nrow(b), nrow(a)))

cat(sprintf("\n%s (%d failing)\n", if (fails == 0) "ALL PASS" else "FAILURES", fails))
quit(status = if (fails == 0) 0 else 1)
