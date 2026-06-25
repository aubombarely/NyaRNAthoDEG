# Changelog

## [v0.1.0] — 2026-06-25

### Added
- Module 1: FastQC pre-trim quality reports (per sample, R1 + R2)
- Module 2: fastp adapter trimming and quality filtering (paired-end)
- Module 3: FastQC post-trim quality reports
- Module 4: MultiQC aggregated QC report (pre + post)
- Module 5: STAR splice-aware alignment (paired-end, sorted BAM)
- Module 6: HTSeq gene-level read counting
- Module 7: DESeq2 full-model + per-stage differential expression
- Module 7: PCA and sample-distance heatmap (VST-normalized)
- Module 7: UpSet plot for DEG overlap across developmental stages
- Module 7: topGO enrichment with Homology-based and FANTASIA annotations
- Checkpoint/resume support for all modules
- `--force`, `--dry_run`, `--skip_module1`–`--skip_module7` flags
- Carbon footprint tracking via codecarbon (opt-out with `--disable_co2_tracking`)
- Run summary JSON with resource usage
- Conda environment file `envs/NyaRNAthoDEG.yaml`
