# Changelog

## [v0.1.1] — 2026-07-07

### Fixed
- **Module 7 (DESeq2)** — AHRD annotation files were not parsed correctly
  when the file begins with the standard `# AHRD-Version X.Y` comment line
  and a blank line before the tab-separated header. `read.delim` was treating
  the comment as the header (1 column) and failing with "more columns than
  column names". Fixed by adding `skip=2` to jump over the two preamble lines.

### Added
- **Module 7 (DESeq2)** — two new summary TSVs written to the output directory:
  - `{prefix}_sample_stats.tsv` — per-sample library metrics: raw total
    counts, genes detected, genes above `--min_counts`, DESeq2 size factor,
    and size-factor-normalised total counts
  - `{prefix}_DEG_summary.tsv` — per-comparison DE metrics: n tested, n DEG,
    n up-regulated, n down-regulated, n genes with NA padj, and % DEG;
    one row for the full model plus one row per stratum when `--deseq_stratify`
    is used

---

## [v0.1.0] — 2026-06-25

### Added

#### Pipeline modules
- **Module 1** — FastQC pre-trim quality reports (per sample, R1 + R2)
- **Module 2** — fastp adapter trimming and quality filtering (paired-end)
- **Module 3** — FastQC post-trim quality reports
- **Module 4** — MultiQC aggregated report (pre-trim + post-trim + fastp)
- **Module 5** — STAR splice-aware alignment (paired-end; sorted, indexed BAM)
- **Module 6** — HTSeq gene-level read counting
- **Module 7** — DESeq2 differential expression + topGO GO enrichment (R)

#### Flexible sample sheet
- Per-sample `StartingModule` (1–7) lets FASTQs, BAMs, and count files coexist
  in the same run without global skip flags
- `FilePath` column accepts two comma-separated FASTQ paths (R1,R2) for
  modules 1–5, a single BAM path for module 6, or a count file path for module 7
- Any column beyond `SampleName`, `StartingModule`, `FilePath`, `Replicate` is
  automatically treated as a DESeq2 condition/covariate factor

#### DESeq2 design
- `--deseq_design` — arbitrary R formula; column names map directly to sample
  sheet condition columns (e.g. `"~ DevelopmentalStage + Temperature"`)
- `--deseq_stratify` — runs a separate DESeq2 analysis per unique level of the
  specified column (e.g. one run per developmental stage)
- `--deseq_contrast` — factor to contrast within each stratum; auto-inferred
  from `--deseq_design` when omitted

#### Module 7 outputs
- Full-model DEG table merged with per-stratum DE results and AHRD + FANTASIA
  annotations (`{prefix}_DEG_annotated.csv`)
- PCA plot (VST-normalized, coloured by condition combination)
- Sample distance heatmap (hierarchical clustering)
- UpSet plot showing DEG overlap across strata
- topGO Fisher enrichment with Homology-based GO annotations (derived from AHRD)
- topGO Fisher enrichment with FANTASIA GO annotations
- Bubble plots (PDF + JPEG) per GO namespace per stratum

#### Pipeline infrastructure
- Checkpoint / resume for every module and every sample — re-running the same
  command skips completed steps automatically
- `--force` to rerun all steps from scratch regardless of existing outputs
- `--dry_run` to validate inputs and print planned steps without executing
- `--skip_module1` – `--skip_module7` as global module overrides
- `mod01_` – `mod07_` directory prefixes in `workdir/` and `results/` so
  outputs self-sort in chronological order
- Run log written to `logs/Run_NyaRNAthoDEG.log` with timestamp, user, host,
  OS, working directory, and full command
- Run summary JSON (`results/{prefix}.run_summary.json`) with all parameters
  and resource usage (wall time, peak RAM, CO₂ emissions)
- Optional carbon footprint tracking via codecarbon
  (`--disable_co2_tracking` to opt out)

#### Project files
- `envs/NyaRNAthoDEG.yaml` — conda environment
- `assets/NyaRNAthoDEG_logo.svg` — project logo
- `LICENSE` — MIT
- `CITATION.cff` — machine-readable citation metadata
- `README.md` — full documentation (modules, sample sheet, options, workflow,
  quicktest, citations)
- `test/make_test_data.py` — reproducible synthetic data generator
  (`random.seed(42)`; 8 samples × 500 genes with baked-in DE signal)
- `test/test_samples.tsv`, `test/htseq_counts/`, `test/test_ahrd.txt`,
  `test/test_fantasia.txt` — committed test inputs
- `test/README.md` — quicktest command and expected output table
