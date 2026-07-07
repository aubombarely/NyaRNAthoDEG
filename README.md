<p align="center">
<img src="assets/NyaRNAthoDEG_logo.svg" width="260"/>
</p>

<p align="center">
<img src="https://img.shields.io/badge/version-v0.1.1-teal"/>
<img src="https://img.shields.io/badge/python-3.10%2B-blue"/>
<img src="https://img.shields.io/badge/R-4.3%2B-276DC3"/>
<img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey"/>
<a href="CHANGELOG.md"><img src="https://img.shields.io/badge/changelog-v0.1.1-orange"/></a>
</p>

---

# NyaRNAthoDEG

**NyaRNAthoDEG** is an end-to-end RNA-seq pipeline covering read quality
control, trimming, splice-aware alignment, gene counting, differential
expression analysis, and GO enrichment.

The name blends **Nyarlathotep** — the *Crawling Chaos*, messenger of the
Outer Gods in Lovecraft's mythos — with **RNA** and **DEG** (differentially
expressed genes), reflecting a pipeline that crawls through raw reads and
delivers expression results.

---

## Overview

| Module | Tool | Input | Output |
|---|---|---|---|
| 1 | FastQC | raw FASTQs | per-sample QC reports (pre-trim) |
| 2 | fastp | raw FASTQs | clean FASTQs + HTML/JSON reports |
| 3 | FastQC | clean FASTQs | per-sample QC reports (post-trim) |
| 4 | MultiQC | FastQC + fastp reports | aggregated QC report |
| 5 | STAR | clean FASTQs + genome index | sorted, indexed BAM files |
| 6 | HTSeq | BAM + GTF | per-sample gene count tables |
| 7 | DESeq2 + topGO (R) | count tables + annotations | DEG tables, PCA, heatmap, UpSet, GO bubble plots |

Each module is independently checkpointed. Individual samples can enter the
pipeline at any module via the `StartingModule` column in the sample sheet,
allowing FASTQs, BAMs, and count files to coexist in the same run.

---

## Requirements

### Conda installation

```bash
conda env create -f envs/NyaRNAthoDEG.yaml
conda activate NyaRNAthoDEG
```

### Dependencies

| Tool | Version | Purpose |
|---|---|---|
| Python | ≥ 3.10 | pipeline orchestration |
| FastQC | any | read quality reports |
| fastp | any | adapter trimming and QC |
| MultiQC | any | aggregated QC report |
| STAR | ≥ 2.7 | splice-aware alignment |
| samtools | ≥ 1.15 | BAM indexing |
| HTSeq | ≥ 2.0 | gene-level read counting |
| R | ≥ 4.3 | DESeq2 + topGO analysis |
| DESeq2 | ≥ 1.40 | differential expression |
| topGO | ≥ 2.52 | GO enrichment |
| ggplot2 | any | PCA and bubble plots |
| pheatmap | any | sample distance heatmap |
| UpSetR | any | UpSet plot |
| optparse | any | R CLI argument parsing |
| codecarbon | any | carbon footprint tracking *(optional)* |

### STAR genome index

Build the genome index once before running the pipeline:

```bash
STAR --runMode genomeGenerate \
     --genomeDir /path/to/star_index \
     --genomeFastaFiles genome.fa \
     --sjdbGTFfile annotation.gtf \
     --runThreadN 8
```

---

## Sample sheet

The sample sheet is a **tab-separated file** with the following columns:

| Column | Required | Description |
|---|---|---|
| `SampleName` | yes | unique sample identifier |
| `StartingModule` | yes | integer 1–7; module at which this sample enters the pipeline |
| `FilePath` | yes | StartingModule 1–5: two comma-separated FASTQ paths (R1,R2) · StartingModule 6: one sorted BAM path · StartingModule 7: one HTSeq count file path |
| `Replicate` | no | replicate ID; stored as metadata, not used in the DESeq2 design |
| *any other columns* | no | treated as condition / covariate factors for DESeq2 (column names must match `--deseq_design`) |

### StartingModule values

| Value | Entry point | FilePath content |
|---|---|---|
| 1 | raw FASTQs → all modules | `R1.fastq.gz,R2.fastq.gz` |
| 2 | raw FASTQs → skip FastQC pre | `R1.fastq.gz,R2.fastq.gz` |
| 3 | clean FASTQs → skip FastQC pre + fastp | `R1.clean.fastq.gz,R2.clean.fastq.gz` |
| 5 | clean FASTQs → skip QC; start at STAR | `R1.clean.fastq.gz,R2.clean.fastq.gz` |
| 6 | sorted BAM → start at HTSeq | `sample.bam` |
| 7 | HTSeq count file → start at DESeq2 | `sample.counts` |

### Example — mixed entry points

```
SampleName      StartingModule  FilePath                                DevelopmentalStage  Temperature  Replicate
Stage1_L_rep1   1               /data/S1L1_R1.fq.gz,/data/S1L1_R2.fq.gz  Stage1           L            1
Stage1_L_rep2   1               /data/S1L2_R1.fq.gz,/data/S1L2_R2.fq.gz  Stage1           L            2
Stage1_H_rep1   6               /data/S1H1.bam                            Stage1            H            1
Stage1_H_rep2   7               /data/S1H2.counts                         Stage1            H            2
Stage2_L_rep1   1               /data/S2L1_R1.fq.gz,/data/S2L1_R2.fq.gz  Stage2           L            1
Stage2_H_rep1   1               /data/S2H1_R1.fq.gz,/data/S2H1_R2.fq.gz  Stage2           H            1
```

---

## Usage

```
python3 scripts/NyaRNAthoDEG.py \
    --sample_sheet   samples.tsv \
    --genome_dir     /path/to/star_index \
    --gtf            /path/to/annotation.gtf \
    --annot_ahrd     /path/to/ahrd_annotation.txt \
    --annot_fantasia /path/to/fantasia_annotation.txt \
    --output         my_run \
    --deseq_design   "~ DevelopmentalStage + Temperature" \
    --deseq_stratify DevelopmentalStage \
    --threads        8
```

---

## Options

### Required inputs

| Option | Description |
|---|---|
| `--sample_sheet` | TSV sample sheet (see format above) |
| `--genome_dir` | STAR genome index directory |
| `--gtf` | GTF annotation file for HTSeq counting |
| `--annot_ahrd` | AHRD annotation TSV (col1=GeneID, col4=Description, col6=GO terms) |
| `--annot_fantasia` | FANTASIA GO annotation in topGO format (`GeneID TAB GO1,GO2,...`) |
| `--output` | Output root directory; its **name** becomes the run prefix |

### DESeq2 design

| Option | Default | Description |
|---|---|---|
| `--deseq_design` | — | R formula string; column names must match condition columns in the sample sheet |
| `--deseq_stratify` | *(none)* | Condition column to stratify by; runs one DESeq2 analysis per level |
| `--deseq_contrast` | *(auto-inferred)* | Factor to contrast within each stratum; defaults to the other term in `--deseq_design` |

#### Design examples

```bash
# Single factor
--deseq_design "~ Temperature"

# Two factors, additive
--deseq_design "~ DevelopmentalStage + Temperature"

# Two factors with interaction
--deseq_design "~ DevelopmentalStage + Temperature + DevelopmentalStage:Temperature"

# Full model + stratified sub-analyses
# Produces: one full-model run AND one per DevelopmentalStage level (~Temperature each)
--deseq_design   "~ DevelopmentalStage + Temperature"
--deseq_stratify DevelopmentalStage

# Equivalent with explicit contrast
--deseq_design   "~ DevelopmentalStage + Temperature"
--deseq_stratify DevelopmentalStage
--deseq_contrast Temperature
```

### Alignment (STAR)

| Option | Default | Description |
|---|---|---|
| `--star_min_intron` | 20 | `--alignIntronMin` passed to STAR |
| `--star_max_intron` | 1000000 | `--alignIntronMax` passed to STAR |

### Trimming (fastp)

| Option | Default | Description |
|---|---|---|
| `--fastp_min_length` | 36 | minimum read length after trimming |
| `--fastp_quality` | 20 | Phred quality threshold |

### Counting (HTSeq)

| Option | Default | Description |
|---|---|---|
| `--htseq_stranded` | `no` | library strandedness: `yes`, `no`, or `reverse` |

### Differential expression / GO

| Option | Default | Description |
|---|---|---|
| `--padj` | 0.05 | adjusted p-value threshold for DEG calling |
| `--min_counts` | 10 | minimum total read count per gene |
| `--go_pval` | 0.01 | p-value cut-off for topGO table output |
| `--go_categories` | `BP` | comma-separated GO namespaces: `BP`, `CC`, `MF` |

### General

| Option | Default | Description |
|---|---|---|
| `--threads` | 4 | threads for FastQC, fastp, and STAR |
| `--format` | `pdf` | plot format(s): `pdf`, `png`, `svg` — comma-separated |
| `--force` | off | rerun all steps even if outputs already exist |
| `--dry_run` | off | validate inputs and print steps, then exit without running |
| `--disable_co2_tracking` | off | disable codecarbon carbon footprint tracking |
| `--version` | — | print version and exit |

### Module skip flags

Each flag is a global override that skips the specified module for **all** samples.
Per-sample entry points are controlled instead via `StartingModule` in the sample sheet.

| Flag | Skips |
|---|---|
| `--skip_module1` | FastQC pre-trim |
| `--skip_module2` | fastp trimming |
| `--skip_module3` | FastQC post-trim |
| `--skip_module4` | MultiQC |
| `--skip_module5` | STAR alignment |
| `--skip_module6` | HTSeq counting |
| `--skip_module7` | DESeq2 + topGO |

---

## Output directory layout

```
{output}/
├── results/
│   ├── mod04_multiqc/
│   │   └── multiqc_report.html
│   ├── mod07_deseq2/
│   │   ├── {prefix}_DEG_annotated.csv       normalized counts + DE results + annotations
│   │   ├── {prefix}_sample_stats.tsv        per-sample library metrics (see below)
│   │   ├── {prefix}_DEG_summary.tsv         per-comparison DEG counts (see below)
│   │   ├── {prefix}_PCA.pdf/.png            PCA plot (VST-normalized)
│   │   ├── {prefix}_SampleDist_heatmap.pdf  sample distance heatmap
│   │   ├── {prefix}_UpSet.pdf               DEG overlap across strata
│   │   ├── {prefix}_HomologyTopGO.txt        topGO mapping (derived from AHRD)
│   │   ├── {prefix}_GOSEA_Homology_BP_{stratum}_table.txt
│   │   ├── {prefix}_GOSEA_Homology_BP_{stratum}.pdf/.jpeg
│   │   ├── {prefix}_GOSEA_FANTASIA_BP_{stratum}_table.txt
│   │   └── {prefix}_GOSEA_FANTASIA_BP_{stratum}.pdf/.jpeg
│   └── {prefix}.run_summary.json
├── workdir/
│   ├── mod01_fastqc_pre/        FastQC HTML + zip (raw reads)
│   ├── mod02_fastp/             clean FASTQs + fastp HTML/JSON reports
│   ├── mod03_fastqc_post/       FastQC HTML + zip (clean reads)
│   ├── mod05_star/{sample}/     STAR BAM files + logs
│   ├── mod06_htseq/             {sample}.htseq.counts files
│   └── {prefix}_sample_info.tsv generated sample table for R
└── logs/
    ├── Run_NyaRNAthoDEG.log
    └── {prefix}.emissions.csv   (if codecarbon is installed)
```

### DEG annotated table columns

The main output `{prefix}_DEG_annotated.csv` contains:

- **Normalized counts** — one column per sample (DESeq2 size-factor normalized)
- **Full-model DE results** — `ALL_baseMean`, `ALL_log2FoldChange`, `ALL_lfcSE`, `ALL_stat`, `ALL_pvalue`, `ALL_padj`
- **Per-stratum DE results** — `{stratum}_baseMean`, `{stratum}_log2FoldChange`, ..., `{stratum}_padj` (one block per level of `--deseq_stratify`)
- **AHRD annotation** — `Description`, `GO_by_Homol`
- **FANTASIA annotation** — `GO_by_FANTASIA`

### Sample statistics table

`{prefix}_sample_stats.tsv` contains one row per sample. Counts are captured
**before** the `--min_counts` filter so the table reflects the original library.

| Column | Description |
|---|---|
| `SampleName` | sample identifier |
| *condition columns* | one column per condition factor from the sample sheet |
| `raw_total_counts` | sum of all raw HTSeq counts |
| `genes_detected` | number of genes with at least one count |
| `genes_above_mincount` | number of genes at or above `--min_counts` |
| `size_factor` | DESeq2 size-factor normalisation coefficient |
| `norm_total_counts` | raw total counts divided by the size factor |

### DEG summary table

`{prefix}_DEG_summary.tsv` contains one row per comparison: the full model
plus one row per stratum when `--deseq_stratify` is used.

| Column | Description |
|---|---|
| `Comparison` | `Full_model` or the stratum name |
| `n_tested` | genes tested after `--min_counts` filtering |
| `n_DEG` | significant DEGs (padj < `--padj`) |
| `n_up` | up-regulated DEGs (log2FoldChange > 0) |
| `n_down` | down-regulated DEGs (log2FoldChange < 0) |
| `n_NA_padj` | genes excluded by DESeq2 (outliers or all-zero counts) |
| `pct_DEG` | percentage of tested genes that are DEG |

---

## Run log

Every run writes `logs/Run_NyaRNAthoDEG.log` containing:

- Date, user, server, OS, working directory, full command
- Timestamped log of every tool call executed
- Checkpoint messages when intermediate outputs are reused

---

## Run summary

`results/{prefix}.run_summary.json` is written at the end of every run:

```json
{
  "date": "2026-06-25 14:30:00",
  "version": "v0.1.0",
  "n_samples": 8,
  "samples": ["Stage1_L_rep1", "..."],
  "condition_cols": ["DevelopmentalStage", "Temperature"],
  "parameters": {
    "deseq_design": "~ DevelopmentalStage + Temperature",
    "deseq_stratify": "DevelopmentalStage",
    "padj": 0.05,
    "..."
  },
  "resource_usage": {
    "wall_clock_s": 1842.3,
    "peak_mem_mb": 4210.1,
    "emissions_kg_CO2eq": 0.000412
  }
}
```

---

## Carbon footprint

When [codecarbon](https://codecarbon.io) is installed, NyaRNAthoDEG tracks
CPU energy consumption and reports CO₂-equivalent emissions in the run summary
and in `logs/{prefix}.emissions.csv`.

```bash
conda install -c conda-forge codecarbon   # optional
```

Disable tracking with `--disable_co2_tracking`.

---

## Recommended workflow

### 1 — Build STAR genome index (once per genome)

```bash
STAR --runMode genomeGenerate \
     --genomeDir /path/to/star_index \
     --genomeFastaFiles genome.fa \
     --sjdbGTFfile annotation.gtf \
     --runThreadN 8
```

### 2 — Prepare the sample sheet

Create a TSV with one row per sample. Set `StartingModule` to the appropriate
value for each sample (1 for raw FASTQs, 6 for BAMs, 7 for count files).

### 3 — Dry run to validate

```bash
python3 scripts/NyaRNAthoDEG.py \
    --sample_sheet   samples.tsv \
    --genome_dir     /path/to/star_index \
    --gtf            annotation.gtf \
    --annot_ahrd     ahrd.txt \
    --annot_fantasia fantasia.txt \
    --output         my_run \
    --deseq_design   "~ Stage + Treatment" \
    --deseq_stratify Stage \
    --dry_run
```

### 4 — Full run

```bash
python3 scripts/NyaRNAthoDEG.py \
    --sample_sheet   samples.tsv \
    --genome_dir     /path/to/star_index \
    --gtf            annotation.gtf \
    --annot_ahrd     ahrd.txt \
    --annot_fantasia fantasia.txt \
    --output         my_run \
    --deseq_design   "~ Stage + Treatment" \
    --deseq_stratify Stage \
    --threads        8 \
    --htseq_stranded reverse \
    --format         pdf,png
```

### 5 — Resume after interruption

Re-run the exact same command. Completed steps are detected via checkpoints
and skipped automatically. Use `--force` to rerun everything from scratch.

### 6 — Run only DESeq2 (skip QC and alignment)

```bash
python3 scripts/NyaRNAthoDEG.py \
    --sample_sheet   samples.tsv \
    --genome_dir     /path/to/star_index \
    --gtf            annotation.gtf \
    --annot_ahrd     ahrd.txt \
    --annot_fantasia fantasia.txt \
    --output         my_run \
    --deseq_design   "~ Stage + Treatment" \
    --skip_module1 --skip_module2 --skip_module3 \
    --skip_module4 --skip_module5 --skip_module6
```

---

## Quicktest

```bash
conda activate NyaRNAthoDEG

# Generate synthetic test data (run once)
python3 test/make_test_data.py

# Run Module 7 only on 8 synthetic samples (< 2 min on a laptop)
python3 scripts/NyaRNAthoDEG.py \
    --sample_sheet   test/test_samples.tsv \
    --genome_dir     test/test_genome_dir \
    --gtf            test/test_annotation.gtf \
    --annot_ahrd     test/test_ahrd.txt \
    --annot_fantasia test/test_fantasia.txt \
    --output         test/test_run \
    --deseq_design   "~ DevelopmentalStage + Temperature" \
    --deseq_stratify DevelopmentalStage \
    --skip_module1 --skip_module2 --skip_module3 \
    --skip_module4 --skip_module5 --skip_module6
```

See [`test/README.md`](test/README.md) for expected outputs and DEG counts.

---

## Third-party tools and citations

| Tool | Reference |
|---|---|
| **FastQC** | Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| **fastp** | Chen S. et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560 |
| **MultiQC** | Ewels P. et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354 |
| **STAR** | Dobin A. et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1), 15–21. https://doi.org/10.1093/bioinformatics/bts635 |
| **samtools** | Danecek P. et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008 |
| **HTSeq** | Putri G. et al. (2022). Analysing high-throughput sequencing data in Python with HTSeq 2.0. *Bioinformatics*, 38(10), 2943–2945. https://doi.org/10.1093/bioinformatics/btac166 |
| **DESeq2** | Love M.I. et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. https://doi.org/10.1186/s13059-014-0550-8 |
| **topGO** | Alexa A. & Rahnenführer J. (2023). topGO: Enrichment Analysis for Gene Ontology. R package. https://bioconductor.org/packages/topGO |
| **codecarbon** | Courty P.-L. et al. (2023). Codecarbon: estimate and track carbon emissions from machine learning computing. https://codecarbon.io |
