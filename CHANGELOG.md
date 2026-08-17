# Changelog

## [v0.2.1] — 2026-08-17

### Fixed

- **Subprocess diagnostic output (including the new "Contrast direction"
  line from v0.2.0) never reached `Run_NyaRNAthoDEG.log`.** `_run()`
  invoked non-data-capturing commands (Module 7's R script, FastQC,
  fastp, STAR, `samtools index`) with `stdout=None`, which inherits
  straight to the terminal -- only the invoked command line itself was
  ever written to the log, not anything the subprocess printed. Found
  immediately after v0.2.0: a real rerun's `Run_NyaRNAthoDEG.log` had no
  trace of the R script's own `[1]`..`[10]` step banners or its
  "Contrast direction" confirmation, defeating the point of always
  printing that line. Fixed by streaming stdout (with stderr merged in,
  so R's own warnings/messages are captured too) line-by-line through
  `_log()` as the subprocess runs -- still visible live in the terminal,
  now also permanently recorded. The data-capturing path
  (`capture_stdout=True`, used by `htseq-count` to write counts straight
  to a file) is untouched -- verified directly that it still returns
  stdout byte-exact with no log formatting or merged stderr, since that
  path relies on `result.stdout` being written verbatim to
  `{sample}.htseq.counts`. Also verified the streaming path captures
  both stdout and stderr lines correctly, and that the on-error
  diagnostic still surfaces the captured output before exiting.

## [v0.2.0] — 2026-08-17

### Fixed

- **log2FoldChange sign was reversed relative to what every report assumed.**
  `NyaRNAthoDEG_DESeq2.R` builds condition-column factors with plain
  `factor()` — no explicit `levels=`, no `relevel()` — and calls `results()`
  with no explicit `contrast=` anywhere (full model or stratified). R's
  `factor()` defaults to alphabetical level order, and DESeq2's `results()`
  with no `contrast=` picks the coefficient for the LAST (non-reference)
  level vs. the first. For a `Mutant` column with values `MPO`/`WT`, that
  silently made `MPO` (alphabetically first) the reference, so every
  `log2FoldChange` was actually `log2(WT/MPO)` — the opposite of what every
  written report assumed ("positive = higher in MPO"). Confirmed against
  `resultsNames(dds)` (literally `Mutant_WT_vs_MPO`) and independently
  against raw counts on two genes with unambiguous, large expression
  differences: both had their reported direction backwards. This affected
  every up/down call in the full model, and every stratum, in every run to
  date -- not a stratum-specific or edge-case bug.

### Added

- **`--deseq_ref_level`**: pins which level of the contrast column is the
  reference/baseline, so `log2FoldChange = log2(other_level / this_level)`
  is deterministic and matches what you tell it to mean, instead of
  silently depending on alphabetical factor ordering. Implemented via an
  explicit `contrast=c(col, other, ref)` in every `results()` call (full
  model and each stratum) rather than relying on DESeq2's default
  coefficient selection.
- **The resolved contrast direction is now always printed to the log**
  (e.g. `Contrast direction : log2( MPO / WT )`), whether or not
  `--deseq_ref_level` is given, specifically so a reversed convention is
  visible in the log instead of only discoverable by cross-checking raw
  counts against a written report, as happened here. Threaded through the
  Python wrapper (`NyaRNAthoDEG.py`) as `--deseq_ref_level`, including
  `--dry_run` reporting and `run_summary.json`.
- Verified directly against real HTSeq count data (not synthetic): with
  `--deseq_ref_level WT`, `NtK326IBMCPv1C09-g573` (MPO counts 21/33/14 vs.
  WT counts 0/2/4 -- MPO unambiguously higher) now correctly returns a
  positive `log2FoldChange`, and p-values/padj are consistent with the
  pre-fix run (only the sign changed, as expected -- the Wald test itself
  is symmetric with respect to which level is called reference).

### Note for prior runs

Every past `DEG_annotated.csv` / report produced before this fix has
`log2FoldChange` sign, and therefore every up/down label, reversed. No
DESeq2 rerun is required to correct existing results -- negate
`log2FoldChange` and swap "up"/"down" labels in place; p-values and padj
are unaffected.

## [v0.1.2] — 2026-07-15

### Fixed
- **Environment installation** — `conda env create` could fail with a
  `LinkError` on `bioconductor-go.db` (`yq: Error starting jq:
  FileNotFoundError`). `jq` is required by that package's post-link script
  but isn't a declared dependency upstream. Added `jq` to
  `envs/NyaRNAthoDEG.yaml` and documented a two-step install
  (`conda create` with `jq` first, then `conda env update`) in the README to
  guarantee link order.

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
