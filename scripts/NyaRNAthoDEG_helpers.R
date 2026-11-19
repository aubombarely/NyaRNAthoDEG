# Helpers shared by NyaRNAthoDEG_DESeq2.R, kept separate so they can be
# exercised by test/regression/test_bugfixes.R without running the pipeline.

# Variance-stabilising transform that also works on small gene sets.
#
# DESeq2's vst() samples `nsub` genes (default 1000) to fit the dispersion
# trend and stops with "less than 'nsub' rows" when the dataset has fewer.
# DESeq2's own recommendation in that case is varianceStabilizingTransformation().
.vst_safe <- function(dds, blind = FALSE, nsub = 1000) {
  if (nrow(dds) < nsub) {
    cat(sprintf("      %d genes (< %d): using varianceStabilizingTransformation()\n",
                nrow(dds), nsub))
    varianceStabilizingTransformation(dds, blind = blind)
  } else {
    vst(dds, blind = blind, nsub = nsub)
  }
}

# GenTable() with a reliably named statistic column.
#
# topGO's GenTable() loses the name of the statistic column when the table has
# a single row: it comes back as the deparsed expression
# "apply(l, 2, format.FUN, dig = 2, eps = 1e-30)" instead of the argument name.
# Any later `tbl$Fishers` is then NULL, and assigning a derived column fails
# with "replacement has 0 rows, data has 1". The statistic is always the last
# column, so name it positionally.
.gentable_fisher <- function(GOdata, gsea_res, topnodes) {
  tbl <- topGO::GenTable(GOdata, Fishers = gsea_res, useLevels = TRUE,
                         topNodes = topnodes)
  if (ncol(tbl) > 0) colnames(tbl)[ncol(tbl)] <- "Fishers"
  tbl
}

# Read an AHRD-style annotation whether or not it carries the two-line preamble.
#
# Real AHRD output starts with "# AHRD-Version ..." and a blank line before the
# header, so the pipeline historically read it with skip=2. A hand-built file
# without that preamble silently lost its header row and the first gene. Decide
# the skip from the file itself instead of assuming.
.read_ahrd <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  skip  <- if (length(first) && grepl("^\\s*#", first)) 2L else 0L
  df    <- read.delim(path, skip = skip, stringsAsFactors = FALSE)
  if (skip == 0L) {
    cat("      note: no '# AHRD-Version' preamble found, reading from line 1\n")
  }
  df
}
