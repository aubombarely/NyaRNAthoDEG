#!/usr/bin/env Rscript
# NyaRNAthoDEG_DESeq2.R — Module 7: differential expression + GO enrichment
# Called by NyaRNAthoDEG.py via Rscript
#
# Key arguments
#   --sample_info     TSV: SampleName, Filename, [condition cols], [Replicate]
#   --htseq_dir       directory containing {SampleName}.htseq.counts files
#   --annot_ahrd      AHRD TSV (col1=GeneID, col4=Description, col6=GO terms)
#   --annot_fantasia  FANTASIA file in topGO mapping format (GeneID TAB GO1,GO2,...)
#   --output_dir      directory for all outputs
#   --prefix          filename prefix
#   --deseq_design    R formula string, e.g. "~ DevelopmentalStage + Temperature"
#   --deseq_stratify  [optional] column to stratify sub-analyses by
#   --deseq_contrast  [optional] column to contrast within each stratum
#                     (auto-inferred from --deseq_design when omitted)
#   --padj            adjusted p-value threshold (default 0.05)
#   --min_counts      minimum total read count per gene (default 10)
#   --go_categories   comma-separated GO namespaces (default BP)
#   --go_pval         topGO p-value cut-off for table output (default 0.01)
#
# DESeq2 design examples
# ----------------------
# Single factor:
#   --deseq_design "~ Temperature"
#
# Two factors, additive:
#   --deseq_design "~ DevelopmentalStage + Temperature"
#
# Two factors with interaction:
#   --deseq_design "~ DevelopmentalStage + Temperature + DevelopmentalStage:Temperature"
#
# Stratified analysis — full model + per-level contrasts:
#   --deseq_design   "~ DevelopmentalStage + Temperature"
#   --deseq_stratify DevelopmentalStage
#   # → full model on all samples, PLUS one run per DevelopmentalStage level
#   #   each contrasting Temperature (L vs H) within that level
#
# Explicit contrast factor:
#   --deseq_design   "~ DevelopmentalStage + Temperature"
#   --deseq_stratify DevelopmentalStage
#   --deseq_contrast Temperature   # auto-inferred when omitted

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(topGO)
  library(ggplot2)
  library(RColorBrewer)
  library(pheatmap)
  library(UpSetR)
})

# ── CLI arguments ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--sample_info",    type="character", default=NULL),
  make_option("--htseq_dir",      type="character", default=NULL),
  make_option("--annot_ahrd",     type="character", default=NULL),
  make_option("--annot_fantasia", type="character", default=NULL),
  make_option("--output_dir",     type="character", default=NULL),
  make_option("--prefix",         type="character", default=NULL),
  make_option("--deseq_design",   type="character", default=NULL,
              help="R formula string, e.g. '~ Stage + Temperature'"),
  make_option("--deseq_stratify", type="character", default="",
              help="Column to stratify sub-analyses by"),
  make_option("--deseq_contrast", type="character", default="",
              help="Column to contrast within each stratum (auto-inferred if omitted)"),
  make_option("--padj",           type="double",    default=0.05),
  make_option("--min_counts",     type="integer",   default=10L),
  make_option("--go_categories",  type="character", default="BP"),
  make_option("--go_pval",        type="double",    default=0.01)
)
opt <- parse_args(OptionParser(option_list=option_list))

for (arg in c("sample_info","htseq_dir","annot_ahrd","annot_fantasia",
              "output_dir","prefix","deseq_design")) {
  if (is.null(opt[[arg]])) stop(sprintf("Missing required argument: --%s", arg))
}

dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)
go_cats <- trimws(strsplit(opt$go_categories, ",")[[1]])

cat("=================================================================\n")
cat("  NyaRNAthoDEG — DESeq2 + topGO module\n")
cat("=================================================================\n")
cat(sprintf("  Sample info    : %s\n", opt$sample_info))
cat(sprintf("  HTSeq dir      : %s\n", opt$htseq_dir))
cat(sprintf("  Output dir     : %s\n", opt$output_dir))
cat(sprintf("  Prefix         : %s\n", opt$prefix))
cat(sprintf("  DESeq2 design  : %s\n", opt$deseq_design))
if (nchar(opt$deseq_stratify) > 0)
  cat(sprintf("  Stratify by    : %s\n", opt$deseq_stratify))
if (nchar(opt$deseq_contrast) > 0)
  cat(sprintf("  Contrast       : %s\n", opt$deseq_contrast))
cat(sprintf("  padj           : %.3f\n", opt$padj))
cat(sprintf("  min_counts     : %d\n",   opt$min_counts))
cat(sprintf("  GO categories  : %s\n",   opt$go_categories))
cat(sprintf("  GO p-val       : %.3f\n", opt$go_pval))
cat("\n")


# ── 1. Load sample table and annotation files ─────────────────────────────────
cat("[1] Loading sample table and annotations\n")

sample_info <- read.delim(opt$sample_info, stringsAsFactors=FALSE, check.names=FALSE)

# Identify condition columns (everything except SampleName, Filename, Replicate)
reserved_cols    <- c("SampleName", "Filename", "Replicate")
condition_cols   <- setdiff(colnames(sample_info), reserved_cols)
has_replicate    <- "Replicate" %in% colnames(sample_info)

# Build the sampleTable expected by DESeqDataSetFromHTSeqCount
sampleTable <- data.frame(
  sampleName = sample_info$SampleName,
  fileName   = sample_info$Filename,
  stringsAsFactors = FALSE
)
for (col in condition_cols) {
  sampleTable[[col]] <- factor(sample_info[[col]])
}
if (has_replicate) {
  sampleTable$Replicate <- factor(sample_info$Replicate)
}

cat(sprintf("  Samples          : %d\n", nrow(sampleTable)))
cat(sprintf("  Condition cols   : %s\n", paste(condition_cols, collapse=", ")))

# Parse the DESeq2 design formula
deseq_formula <- tryCatch(
  as.formula(opt$deseq_design),
  error=function(e) stop(sprintf("Invalid --deseq_design formula '%s': %s",
                                 opt$deseq_design, conditionMessage(e)))
)
design_terms <- all.vars(deseq_formula)
cat(sprintf("  Design formula   : %s\n", opt$deseq_design))
cat(sprintf("  Design terms     : %s\n", paste(design_terms, collapse=", ")))

# Validate that design terms exist in the sample table
missing_terms <- setdiff(design_terms, colnames(sampleTable))
if (length(missing_terms) > 0) {
  stop(sprintf("Design term(s) not found in sample sheet: %s",
               paste(missing_terms, collapse=", ")))
}

# Resolve contrast column for stratified analysis
do_stratify   <- nchar(opt$deseq_stratify) > 0
contrast_col  <- NULL
if (do_stratify) {
  if (!(opt$deseq_stratify %in% colnames(sampleTable))) {
    stop(sprintf("--deseq_stratify column '%s' not found in sample sheet",
                 opt$deseq_stratify))
  }
  contrast_col <- if (nchar(opt$deseq_contrast) > 0) {
    opt$deseq_contrast
  } else {
    setdiff(design_terms, opt$deseq_stratify)[1]
  }
  cat(sprintf("  Stratify by      : %s\n", opt$deseq_stratify))
  cat(sprintf("  Contrast within  : %s\n", contrast_col))
  strata <- levels(factor(sampleTable[[opt$deseq_stratify]]))
  cat(sprintf("  Strata           : %s\n", paste(strata, collapse=", ")))
}

# AHRD annotation: columns 1 (GeneID), 4 (Description), 6 (GO)
annot_ahrd_raw     <- read.delim(opt$annot_ahrd, stringsAsFactors=FALSE)
annot_by_hom_table <- annot_ahrd_raw[, c(1, 4, 6)]
colnames(annot_by_hom_table) <- c("GeneID", "Description", "GO_by_Homol")

hom_tgo_path <- file.path(opt$output_dir,
                           paste0(opt$prefix, "_HomologyTopGO.txt"))
write.table(annot_ahrd_raw[, c(1, 6)], hom_tgo_path,
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# FANTASIA annotation (topGO format + two-column table)
annot_fantasia_raw <- read.delim(opt$annot_fantasia,
                                 header=FALSE, stringsAsFactors=FALSE)
colnames(annot_fantasia_raw) <- c("GeneID", "GO_by_FANTASIA")

cat(sprintf("  AHRD genes       : %d\n", nrow(annot_ahrd_raw)))
cat(sprintf("  FANTASIA genes   : %d\n", nrow(annot_fantasia_raw)))


# ── 2. Load HTSeq count data ──────────────────────────────────────────────────
cat("\n[2] Loading HTSeq count data\n")
dds_full <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  directory   = opt$htseq_dir,
  design      = deseq_formula
)
total_genes <- nrow(dds_full)
total_count <- sum(counts(dds_full))
cat(sprintf("  Total genes  : %d\n", total_genes))
cat(sprintf("  Total counts : %d\n", total_count))


# ── 3. Filter low-count genes ─────────────────────────────────────────────────
cat(sprintf("\n[3] Filtering genes with total counts < %d\n", opt$min_counts))
dds_full <- dds_full[rowSums(counts(dds_full)) >= opt$min_counts, ]
cat(sprintf("  Genes after filter : %d  (%.1f%% retained)\n",
            nrow(dds_full), 100 * nrow(dds_full) / total_genes))


# ── 4. Full-model DESeq2 ──────────────────────────────────────────────────────
cat(sprintf("\n[4] Full-model DESeq2 (%s)\n", opt$deseq_design))
dds_full  <- DESeq(dds_full, quiet=TRUE)
res_full  <- results(dds_full, alpha=opt$padj)
n_sig_all <- sum(res_full$padj < opt$padj, na.rm=TRUE)
cat(sprintf("  DEGs (padj < %.2f) : %d\n", opt$padj, n_sig_all))


# ── 5. Stratified sub-analyses ────────────────────────────────────────────────
strat_res_list <- list()

if (do_stratify) {
  cat(sprintf("\n[5] Stratified analysis — per level of '%s', contrast '%s'\n",
              opt$deseq_stratify, contrast_col))
  strat_formula <- as.formula(paste("~", contrast_col))

  for (stratum in strata) {
    cat(sprintf("  Stratum: %s\n", stratum))
    sel <- colData(dds_full)[[opt$deseq_stratify]] == stratum
    dds_s        <- dds_full[, sel]
    design(dds_s) <- strat_formula
    dds_s         <- DESeq(dds_s, quiet=TRUE)
    res_s         <- results(dds_s, alpha=opt$padj)
    n_sig_s       <- sum(res_s$padj < opt$padj, na.rm=TRUE)
    cat(sprintf("    DEGs (padj < %.2f) : %d\n", opt$padj, n_sig_s))
    strat_res_list[[stratum]] <- res_s
  }
} else {
  cat("\n[5] No stratification requested — skipping\n")
}


# ── 6. Build merged output table ─────────────────────────────────────────────
cat("\n[6] Building merged annotated table\n")
count_tbl        <- as.data.frame(counts(dds_full, normalized=TRUE))
count_tbl$GeneID <- rownames(count_tbl)
merge_tbl        <- count_tbl

# Full model results
de_full_df <- as.data.frame(res_full)
colnames(de_full_df) <- paste0("ALL_", colnames(de_full_df))
de_full_df$GeneID    <- rownames(de_full_df)
merge_tbl <- merge(merge_tbl, de_full_df, by="GeneID")

# Per-stratum results
for (stratum in names(strat_res_list)) {
  stag_tag   <- gsub("[^A-Za-z0-9]", "", as.character(stratum))
  de_s_df    <- as.data.frame(strat_res_list[[stratum]])
  colnames(de_s_df) <- paste0(stag_tag, "_", colnames(de_s_df))
  de_s_df$GeneID    <- rownames(de_s_df)
  merge_tbl <- merge(merge_tbl, de_s_df, by="GeneID")
}

merge_tbl <- merge(merge_tbl, annot_by_hom_table, by="GeneID", all.x=TRUE)
merge_tbl <- merge(merge_tbl, annot_fantasia_raw,  by="GeneID", all.x=TRUE)

out_csv <- file.path(opt$output_dir, paste0(opt$prefix, "_DEG_annotated.csv"))
write.csv(merge_tbl, out_csv, row.names=FALSE)
cat(sprintf("  Written : %s\n", out_csv))


# ── 7. PCA plot ───────────────────────────────────────────────────────────────
cat("\n[7] PCA plot (VST)\n")
vsdata   <- vst(dds_full, blind=FALSE)
pca_data <- plotPCA(vsdata, intgroup=condition_cols, returnData=TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

pca_fig <- ggplot(pca_data, aes(x=PC1, y=PC2, color=group)) +
  geom_point(size=3) +
  xlab(sprintf("PC1: %d%% variance", pct_var[1])) +
  ylab(sprintf("PC2: %d%% variance", pct_var[2])) +
  theme_bw(base_size=12) +
  ggtitle(sprintf("PCA — VST-normalized counts (%s)", paste(condition_cols, collapse=" + ")))

pca_base <- file.path(opt$output_dir, paste0(opt$prefix, "_PCA"))
ggsave(paste0(pca_base, ".pdf"), pca_fig, width=7, height=5)
ggsave(paste0(pca_base, ".png"), pca_fig, width=7, height=5, dpi=150)
cat(sprintf("  Written : %s.pdf / .png\n", pca_base))


# ── 8. Sample distance heatmap ────────────────────────────────────────────────
cat("\n[8] Sample distance heatmap\n")
sample_dists <- dist(t(assay(vsdata)))
dist_matrix  <- as.matrix(sample_dists)
row_labels   <- apply(as.data.frame(colData(vsdata)[, condition_cols, drop=FALSE]),
                      1, paste, collapse="-")
rownames(dist_matrix) <- row_labels
colnames(dist_matrix) <- NULL
heat_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

heatmap_path <- file.path(opt$output_dir,
                           paste0(opt$prefix, "_SampleDist_heatmap.pdf"))
pdf(heatmap_path, width=8, height=7)
pheatmap(dist_matrix,
         clustering_distance_rows=sample_dists,
         clustering_distance_cols=sample_dists,
         col=heat_colors, show_rownames=TRUE, show_colnames=TRUE)
invisible(dev.off())
cat(sprintf("  Written : %s\n", heatmap_path))


# ── 9. UpSet plot ─────────────────────────────────────────────────────────────
cat("\n[9] UpSet plot\n")
sig_gene_lists <- list()

if (length(strat_res_list) > 0) {
  for (stratum in names(strat_res_list)) {
    stag_tag <- gsub("[^A-Za-z0-9]", "", as.character(stratum))
    padj_col <- paste0(stag_tag, "_padj")
    if (padj_col %in% colnames(merge_tbl)) {
      sig_ids <- merge_tbl[
        !is.na(merge_tbl[[padj_col]]) & merge_tbl[[padj_col]] < opt$padj,
        "GeneID"
      ]
      sig_gene_lists[[as.character(stratum)]] <- sig_ids
      cat(sprintf("  %s : %d significant genes\n", stratum, length(sig_ids)))
    }
  }
} else {
  # No stratification: one entry for the full model
  sig_ids <- merge_tbl[
    !is.na(merge_tbl$ALL_padj) & merge_tbl$ALL_padj < opt$padj,
    "GeneID"
  ]
  sig_gene_lists[["ALL"]] <- sig_ids
  cat(sprintf("  ALL : %d significant genes\n", length(sig_ids)))
}

if (length(sig_gene_lists) >= 2) {
  upset_path <- file.path(opt$output_dir, paste0(opt$prefix, "_UpSet.pdf"))
  pdf(upset_path, width=10, height=6)
  print(upset(fromList(sig_gene_lists),
              nsets=length(sig_gene_lists), order.by="freq"))
  invisible(dev.off())
  cat(sprintf("  Written : %s\n", upset_path))
} else {
  cat("  Skipping UpSet (fewer than 2 comparison groups)\n")
}


# ── 10. topGO enrichment ──────────────────────────────────────────────────────

.bubble_plot <- function(tbl, fig_title) {
  tbl$log10F <- suppressWarnings(as.numeric(tbl$log10F))
  tbl <- tbl[!is.na(tbl$log10F) & is.finite(tbl$log10F), ]
  if (nrow(tbl) == 0) return(NULL)
  ggplot(tbl, aes(x=Term, y=log10F, size=Significant, fill=log10F)) +
    geom_point(shape=21) +
    scale_fill_continuous(low="royalblue", high="tomato") +
    xlab("") + ylab("-log10(p value)") +
    labs(title=fig_title) +
    geom_hline(yintercept=c(-log10(0.05), -log10(0.01), -log10(0.001)),
               linetype="dotted", colour="black", linewidth=0.8) +
    theme_bw(base_size=12) +
    theme(axis.text.y=element_text(vjust=0.5), legend.position="right") +
    coord_flip()
}

.run_topgo <- function(sig_lists, geneID2GO, go_cat,
                       annot_label, output_dir, prefix, go_pval) {
  for (grp in names(sig_lists)) {
    grp_tag    <- gsub("[^A-Za-z0-9]", "", grp)
    gene_names <- names(geneID2GO)
    gene_vec   <- factor(as.integer(gene_names %in% sig_lists[[grp]]))
    names(gene_vec) <- gene_names
    fig_title  <- sprintf("%s %s-based GSEA — %s DE genes", go_cat, annot_label, grp)
    cat(sprintf("    %s\n", fig_title))

    tryCatch({
      GOdata   <- new("topGOdata", ontology=go_cat, allGenes=gene_vec,
                      annot=annFUN.gene2GO, gene2GO=geneID2GO)
      gsea_res <- runTest(GOdata, statistic="fisher")
      pval_df  <- as.data.frame(score(gsea_res))
      topnodes <- max(1L, sum(pval_df < go_pval))
      gsea_tbl <- GenTable(GOdata, Fishers=gsea_res,
                           useLevels=TRUE, topNodes=topnodes)
      gsea_tbl$log10F <- -log10(suppressWarnings(as.numeric(gsea_tbl$Fishers)))

      out_base <- file.path(output_dir,
        sprintf("%s_GOSEA_%s_%s_%s", prefix, annot_label, go_cat, grp_tag))
      write.table(gsea_tbl, paste0(out_base, "_table.txt"),
                  sep="\t", quote=FALSE, row.names=FALSE)

      fig <- .bubble_plot(gsea_tbl, fig_title)
      if (!is.null(fig)) {
        h <- max(4, nrow(gsea_tbl) * 0.3 + 2)
        pdf(paste0(out_base, ".pdf"), width=10, height=h)
        print(fig); invisible(dev.off())
        jpeg(paste0(out_base, ".jpeg"), width=700, height=round(h * 80))
        print(fig); invisible(dev.off())
      }
      cat(sprintf("      Written : %s\n", out_base))
    }, error=function(e) {
      cat(sprintf("      WARNING: topGO failed for %s / %s: %s\n",
                  grp, go_cat, conditionMessage(e)))
    })
  }
}

cat("\n[10] topGO enrichment analysis\n")

if (length(sig_gene_lists) > 0) {
  cat("  >> Homology-based GO annotations\n")
  geneID2GO_hom <- tryCatch(readMappings(hom_tgo_path), error=function(e) {
    cat(sprintf("  WARNING: could not load Homology topGO file: %s\n", conditionMessage(e)))
    NULL
  })
  if (!is.null(geneID2GO_hom)) {
    for (go_cat in go_cats)
      .run_topgo(sig_gene_lists, geneID2GO_hom, go_cat,
                 "Homology", opt$output_dir, opt$prefix, opt$go_pval)
  }

  cat("  >> FANTASIA-based GO annotations\n")
  geneID2GO_fan <- tryCatch(readMappings(opt$annot_fantasia), error=function(e) {
    cat(sprintf("  WARNING: could not load FANTASIA topGO file: %s\n", conditionMessage(e)))
    NULL
  })
  if (!is.null(geneID2GO_fan)) {
    for (go_cat in go_cats)
      .run_topgo(sig_gene_lists, geneID2GO_fan, go_cat,
                 "FANTASIA", opt$output_dir, opt$prefix, opt$go_pval)
  }
} else {
  cat("  No significant gene lists — skipping topGO\n")
}


# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n=================================================================\n")
cat("  DESeq2 module complete\n")
cat(sprintf("  Total genes (input)  : %d\n", total_genes))
cat(sprintf("  Genes after filter   : %d\n", nrow(dds_full)))
cat(sprintf("  Full-model DEGs      : %d\n", n_sig_all))
if (length(strat_res_list) > 0) {
  for (stratum in names(strat_res_list)) {
    stag_tag <- gsub("[^A-Za-z0-9]", "", as.character(stratum))
    padj_col <- paste0(stag_tag, "_padj")
    if (padj_col %in% colnames(merge_tbl)) {
      n <- sum(!is.na(merge_tbl[[padj_col]]) & merge_tbl[[padj_col]] < opt$padj)
      cat(sprintf("  DEGs %-22s: %d\n", stratum, n))
    }
  }
}
cat(sprintf("  Output directory     : %s\n", opt$output_dir))
cat("=================================================================\n")
