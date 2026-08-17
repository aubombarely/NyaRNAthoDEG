#!/usr/bin/env python3
"""
NyaRNAthoDEG — Nyarlathotep-inspired RNA-seq pipeline
Nyarlathotep (the Crawling Chaos) + RNA + DEG

Modules
  1  FastQC    raw-read quality reports
  2  fastp     adapter trimming and quality filtering
  3  FastQC    post-trim quality reports
  4  MultiQC   aggregated QC report
  5  STAR      splice-aware alignment
  6  HTSeq     gene-level read counting
  7  DESeq2/R  differential expression + GO enrichment (topGO)

Sample sheet (TSV) — required columns:
  SampleName      unique sample identifier
  StartingModule  integer 1–7; module at which to enter the pipeline
  FilePath        StartingModule 1–5 → two comma-separated FASTQ paths (R1,R2)
                  StartingModule 6   → one BAM path (sorted, indexed)
                  StartingModule 7   → one HTSeq count file path
  Replicate       replicate ID (optional; metadata only)
  <any other columns> → treated as condition/covariate factors for DESeq2

Example sample sheet:
  SampleName  StartingModule  FilePath                             DevelopmentalStage  Temperature  Replicate
  sampleA     1               /data/A_R1.fq.gz,/data/A_R2.fq.gz  Stage1              L            1
  sampleB     6               /data/B.bam                          Stage1              H            1
  sampleC     7               /data/C.counts                       Stage2              L            2
"""

VERSION = "v0.2.1"

import argparse, csv, getpass, json, os, platform, warnings
import resource, shutil, subprocess, sys, time
from collections import Counter
from datetime import datetime
from pathlib import Path

_LOG_FH = None

# columns consumed by the pipeline; everything else becomes a DESeq2 factor
_RESERVED_COLS = {"SampleName", "StartingModule", "FilePath", "Replicate"}


# ── logging ───────────────────────────────────────────────────────────────────

def _log(msg: str) -> None:
    ts   = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, file=sys.stderr)
    if _LOG_FH is not None:
        print(line, file=_LOG_FH, flush=True)

def _banner(title: str) -> None:
    bar = "─" * (len(title) + 4)
    _log(f"┌{bar}┐")
    _log(f"│  {title}  │")
    _log(f"└{bar}┘")


_QUOTE_LINES = [
    "\"And through this revolting graveyard of the universe the muffled,",
    " maddening beating of drums, and thin, monotonous whine of blasphemous",
    " flutes from inconceivable, unlighted chambers beyond Time; the blind,",
    " voiceless, mindless gargoyles whose soul and messenger is Nyarlathotep.\"",
    "                                  — H.P. Lovecraft, Nyarlathotep (1920)",
]

def _print_quote() -> None:
    width = max(len(l) for l in _QUOTE_LINES) + 4
    border = "─" * width
    _log(f"┌{border}┐")
    for line in _QUOTE_LINES:
        padding = width - len(line) - 1
        _log(f"│ {line}{' ' * padding}│")
    _log(f"└{border}┘")


# ── external tool helpers ─────────────────────────────────────────────────────

def _require_tool(name: str) -> str:
    tool = shutil.which(name)
    if tool is None:
        print(f"ERROR: '{name}' not found in PATH.\n"
              f"       Install with:  conda install -c bioconda {name}",
              file=sys.stderr)
        sys.exit(1)
    return tool

def _run(cmd: list, capture_stdout: bool = False,
         env: dict = None, cwd: Path = None) -> subprocess.CompletedProcess:
    _log(f"  $ {' '.join(str(c) for c in cmd)}")
    if capture_stdout:
        # stdout is DATA the caller writes verbatim to a file (e.g. htseq-count
        # counts) -- must stay byte-exact, so it's captured as before with no
        # streaming/logging and no merging with stderr.
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True, env=env, cwd=cwd)
        if result.returncode != 0:
            print(f"ERROR: command failed (exit {result.returncode}):\n"
                  f"{result.stderr[-3000:]}", file=sys.stderr)
            sys.exit(1)
        return result

    # stdout here is diagnostic/progress text for humans (e.g. the R script's
    # own [1]..[10] step banners) -- previously inherited straight to the
    # terminal via stdout=None and never reached Run_NyaRNAthoDEG.log, so
    # anything a subprocess printed (including NyaRNAthoDEG_DESeq2.R's
    # "Contrast direction" line -- added specifically so a reversed sign
    # convention would be visible in the permanent log, not just live in a
    # terminal someone may not have been watching) was lost from the log.
    # Stream it line-by-line through _log() instead, so it's both still
    # visible live and now actually recorded. stderr is merged into the same
    # stream (R's own warnings/messages go there) so the on-error diagnostic
    # below captures everything, not just half of it.
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, env=env, cwd=cwd, bufsize=1)
    captured_lines = []
    for line in proc.stdout:
        line = line.rstrip("\n")
        captured_lines.append(line)
        _log(f"    {line}")
    proc.wait()
    if proc.returncode != 0:
        print(f"ERROR: command failed (exit {proc.returncode}):\n"
              + "\n".join(captured_lines)[-3000:], file=sys.stderr)
        sys.exit(1)
    return subprocess.CompletedProcess(cmd, proc.returncode)


# ── checkpoint ────────────────────────────────────────────────────────────────

def _checkpoint(path: Path, label: str, force: bool) -> bool:
    if not force and path.exists() and path.stat().st_size > 0:
        _log(f"  [checkpoint] {label} — {path.name} already exists, skipping")
        return True
    return False


# ── input validation ──────────────────────────────────────────────────────────

def _validate_inputs(pairs: list) -> None:
    ok = True
    for flag, path in pairs:
        if not path.exists():
            print(f"ERROR: {flag} not found: {path}", file=sys.stderr)
            ok = False
    if not ok:
        sys.exit(1)


# ── sample sheet ──────────────────────────────────────────────────────────────

def parse_sample_sheet(path: Path) -> tuple:
    """
    Returns (samples, condition_cols).

    Each sample dict contains:
      sample_name, starting_module, file_paths, replicate, conditions
      + r1/r2 (sm <= 5), bam (sm == 6), or count_file (sm == 7)
    """
    required = {"SampleName", "StartingModule", "FilePath"}
    samples  = []

    with open(path) as fh:
        reader     = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing    = required - set(fieldnames)
        if missing:
            print(f"ERROR: sample sheet missing columns: {', '.join(sorted(missing))}",
                  file=sys.stderr)
            sys.exit(1)

        condition_cols = [c for c in fieldnames if c not in _RESERVED_COLS]
        has_replicate  = "Replicate" in fieldnames

        for row in reader:
            name = row["SampleName"].strip()
            try:
                sm = int(row["StartingModule"].strip())
            except ValueError:
                print(f"ERROR: StartingModule must be an integer for sample '{name}'",
                      file=sys.stderr)
                sys.exit(1)

            if sm not in range(1, 8):
                print(f"ERROR: StartingModule must be 1–7 for sample '{name}', got {sm}",
                      file=sys.stderr)
                sys.exit(1)

            file_paths = [Path(p.strip()) for p in row["FilePath"].split(",")]

            if sm <= 5 and len(file_paths) != 2:
                print(f"ERROR: '{name}' StartingModule={sm} requires two comma-separated "
                      f"FASTQ paths (R1,R2) in FilePath", file=sys.stderr)
                sys.exit(1)
            if sm in (6, 7) and len(file_paths) != 1:
                print(f"ERROR: '{name}' StartingModule={sm} requires exactly one file path "
                      f"in FilePath", file=sys.stderr)
                sys.exit(1)

            s = {
                "sample_name":     name,
                "starting_module": sm,
                "file_paths":      file_paths,
                "replicate":       row.get("Replicate", "").strip() if has_replicate else "",
                "conditions":      {c: row[c].strip() for c in condition_cols},
            }

            if sm <= 5:
                s["r1"] = file_paths[0]
                s["r2"] = file_paths[1]
            elif sm == 6:
                s["bam"] = file_paths[0]
            else:
                s["count_file"] = file_paths[0]

            samples.append(s)

    if not samples:
        print("ERROR: sample sheet is empty", file=sys.stderr)
        sys.exit(1)

    _log(f"  Loaded {len(samples)} sample(s) from {path.name}")
    _log(f"  Condition columns : {', '.join(condition_cols) if condition_cols else 'none'}")
    for sm, n in sorted(Counter(s["starting_module"] for s in samples).items()):
        _log(f"  StartingModule={sm} : {n} sample(s)")

    return samples, condition_cols


# ── module 1 & 3: FastQC ─────────────────────────────────────────────────────

def run_fastqc(samples: list, outdir: Path, threads: int,
               force: bool, label: str = "pre") -> None:
    """label='pre' uses r1/r2; label='post' uses r1_clean/r2_clean."""
    if not samples:
        return
    fastqc = _require_tool("fastqc")
    outdir.mkdir(parents=True, exist_ok=True)
    for s in samples:
        pairs = ([(s["r1"], "R1"), (s["r2"], "R2")] if label == "pre"
                 else [(s["r1_clean"], "R1"), (s["r2_clean"], "R2")])
        for read_path, rn in pairs:
            stem = read_path.name
            for ext in (".fastq.gz", ".fastq", ".fq.gz", ".fq"):
                stem = stem.replace(ext, "")
            chk = outdir / f"{stem}_fastqc.html"
            if _checkpoint(chk, f"fastqc-{label} {s['sample_name']} {rn}", force):
                continue
            _run([fastqc, str(read_path.resolve()),
                  "--outdir", str(outdir),
                  "--threads", str(threads),
                  "--quiet"])


# ── module 2: fastp ───────────────────────────────────────────────────────────

def run_fastp(samples: list, workdir: Path, threads: int,
              min_length: int, quality: int, force: bool) -> None:
    """
    Samples with StartingModule <= 2 are trimmed by fastp.
    Samples with StartingModule 3–5 already have clean FASTQs — r1_clean/r2_clean
    is set to their provided paths without running fastp.
    """
    fastp_bin = _require_tool("fastp")
    fastp_dir = workdir / "mod02_fastp"
    fastp_dir.mkdir(parents=True, exist_ok=True)
    for s in samples:
        sm   = s["starting_module"]
        name = s["sample_name"]
        if sm <= 2:
            r1_out = fastp_dir / f"{name}_R1.clean.fastq.gz"
            if not _checkpoint(r1_out, f"fastp {name}", force):
                r2_out = fastp_dir / f"{name}_R2.clean.fastq.gz"
                _run([fastp_bin,
                      "--in1",  str(s["r1"].resolve()),
                      "--in2",  str(s["r2"].resolve()),
                      "--out1", str(r1_out),
                      "--out2", str(r2_out),
                      "--html", str(fastp_dir / f"{name}_fastp.html"),
                      "--json", str(fastp_dir / f"{name}_fastp.json"),
                      "--thread",                  str(threads),
                      "--length_required",         str(min_length),
                      "--qualified_quality_phred", str(quality),
                      "--detect_adapter_for_pe"])
            s["r1_clean"] = r1_out
            s["r2_clean"] = fastp_dir / f"{name}_R2.clean.fastq.gz"
        elif sm <= 5:
            s["r1_clean"] = s["r1"]
            s["r2_clean"] = s["r2"]


# ── module 4: MultiQC ─────────────────────────────────────────────────────────

def run_multiqc(input_dirs: list, outdir: Path, force: bool) -> None:
    multiqc = _require_tool("multiqc")
    outdir.mkdir(parents=True, exist_ok=True)
    chk = outdir / "multiqc_report.html"
    if _checkpoint(chk, "multiqc", force):
        return
    existing = [str(d) for d in input_dirs if d.exists()]
    if not existing:
        _log("  WARNING: no FastQC/fastp directories found for MultiQC — skipping")
        return
    _run([multiqc, "--force", "--outdir", str(outdir)] + existing)


# ── module 5: STAR ────────────────────────────────────────────────────────────

def run_star(samples: list, genome_dir: Path, workdir: Path,
             threads: int, min_intron: int, max_intron: int, force: bool) -> None:
    star     = _require_tool("STAR")
    samtools = _require_tool("samtools")
    star_dir = workdir / "mod05_star"
    for s in samples:
        if s["starting_module"] > 5:
            continue
        name       = s["sample_name"]
        sample_dir = star_dir / name
        sample_dir.mkdir(parents=True, exist_ok=True)
        bam = sample_dir / "Aligned.sortedByCoord.out.bam"
        if not _checkpoint(bam, f"STAR {name}", force):
            _run([star,
                  "--runThreadN",            str(threads),
                  "--genomeDir",             str(genome_dir),
                  "--readFilesIn",           str(s["r1_clean"]), str(s["r2_clean"]),
                  "--readFilesCommand",      "zcat",
                  "--outSAMtype",            "BAM", "SortedByCoordinate",
                  "--outFileNamePrefix",     str(sample_dir) + "/",
                  "--outSAMattributes",      "NH", "HI", "AS", "NM",
                  "--alignIntronMin",        str(min_intron),
                  "--alignIntronMax",        str(max_intron),
                  "--outFilterMultimapNmax", "10"],
                 cwd=sample_dir)
            _run([samtools, "index", str(bam)])
        s["bam"] = bam


# ── module 6: HTSeq ───────────────────────────────────────────────────────────

def run_htseq(samples: list, gtf: Path, workdir: Path,
              stranded: str, force: bool) -> None:
    """
    Skips samples with StartingModule == 7 (count file provided directly).
    Samples with StartingModule == 6 use their provided BAM from file_paths.
    """
    htseq     = _require_tool("htseq-count")
    htseq_dir = workdir / "mod06_htseq"
    htseq_dir.mkdir(parents=True, exist_ok=True)
    for s in samples:
        if s["starting_module"] == 7:
            continue
        name       = s["sample_name"]
        count_file = htseq_dir / f"{name}.htseq.counts"
        if not _checkpoint(count_file, f"htseq {name}", force):
            result = _run([htseq,
                           "--format",   "bam",
                           "--order",    "pos",
                           "--stranded", stranded,
                           str(s["bam"]),
                           str(gtf.resolve())],
                          capture_stdout=True)
            count_file.write_text(result.stdout)
        s["count_file"] = count_file


# ── module 7: DESeq2 + topGO (R) ─────────────────────────────────────────────

def run_deseq2(samples: list, condition_cols: list,
               workdir: Path, results_dir: Path,
               annot_ahrd: Path, annot_fantasia: Path,
               prefix: str, padj: float, min_counts: int,
               deseq_design: str, deseq_stratify: str, deseq_contrast: str,
               deseq_ref_level: str,
               go_categories: str, go_pval: float,
               r_script: Path, force: bool) -> None:
    rscript   = _require_tool("Rscript")
    htseq_dir = workdir / "mod06_htseq"
    htseq_dir.mkdir(parents=True, exist_ok=True)

    # Bring StartingModule==7 count files into htseq_dir so DESeq2 sees one directory
    for s in samples:
        if s["starting_module"] == 7:
            dest = htseq_dir / f"{s['sample_name']}.htseq.counts"
            if not dest.exists():
                shutil.copy2(s["count_file"], dest)
            s["count_file"] = dest

    # Write sample_info TSV for the R script
    has_replicate    = any(s["replicate"] for s in samples)
    sample_info_path = workdir / f"{prefix}_sample_info.tsv"
    with open(sample_info_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        header = ["SampleName", "Filename"] + condition_cols
        if has_replicate:
            header.append("Replicate")
        writer.writerow(header)
        for s in samples:
            row = ([s["sample_name"], f"{s['sample_name']}.htseq.counts"]
                   + [s["conditions"].get(c, "") for c in condition_cols])
            if has_replicate:
                row.append(s["replicate"])
            writer.writerow(row)

    deseq2_dir = results_dir / "mod07_deseq2"
    deseq2_dir.mkdir(parents=True, exist_ok=True)
    chk = deseq2_dir / f"{prefix}_DEG_annotated.csv"
    if _checkpoint(chk, "DESeq2+topGO", force):
        return

    cmd = [rscript, str(r_script),
           "--sample_info",    str(sample_info_path),
           "--htseq_dir",      str(htseq_dir),
           "--annot_ahrd",     str(annot_ahrd),
           "--annot_fantasia", str(annot_fantasia),
           "--output_dir",     str(deseq2_dir),
           "--prefix",         prefix,
           "--deseq_design",   deseq_design,
           "--padj",           str(padj),
           "--min_counts",     str(min_counts),
           "--go_categories",  go_categories,
           "--go_pval",        str(go_pval)]
    if deseq_stratify:
        cmd += ["--deseq_stratify", deseq_stratify]
    if deseq_contrast:
        cmd += ["--deseq_contrast", deseq_contrast]
    if deseq_ref_level:
        cmd += ["--deseq_ref_level", deseq_ref_level]

    _run(cmd)


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    warnings.filterwarnings("ignore", category=FutureWarning, message=".*pynvml.*")
    _print_quote()
    ap = argparse.ArgumentParser(
        prog="NyaRNAthoDEG",
        description="RNA-seq pipeline: QC → fastp → STAR → HTSeq → DESeq2/topGO",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
DESeq2 design examples
----------------------
Single factor:
  --deseq_design "~ Temperature"

Two factors, additive:
  --deseq_design "~ DevelopmentalStage + Temperature"

Two factors with interaction:
  --deseq_design "~ DevelopmentalStage + Temperature + DevelopmentalStage:Temperature"

Stratified analysis (full model + per-level contrasts):
  --deseq_design    "~ DevelopmentalStage + Temperature"
  --deseq_stratify  DevelopmentalStage
  # → full model on all samples, plus one DESeq2 run per DevelopmentalStage level
  #   each contrasting Temperature (L vs H) within that level

Specify the contrast factor explicitly:
  --deseq_design    "~ DevelopmentalStage + Temperature"
  --deseq_stratify  DevelopmentalStage
  --deseq_contrast  Temperature
  # equivalent to the above; --deseq_contrast is auto-inferred when omitted

StartingModule values
---------------------
  1  raw FASTQs (R1,R2) → all modules
  2  raw FASTQs → skip FastQC pre; start at fastp
  3  clean FASTQs (already trimmed) → skip FastQC pre + fastp
  5  clean FASTQs → skip QC; start at STAR
  6  sorted BAM → skip QC + STAR; start at HTSeq
  7  HTSeq count file → skip QC + STAR + HTSeq; start at DESeq2
""",
    )

    # inputs
    ap.add_argument("--sample_sheet",   required=True, type=Path,
                    help="TSV: SampleName, StartingModule, FilePath, [conditions...], [Replicate]")
    ap.add_argument("--genome_dir",     required=True, type=Path,
                    help="STAR genome index directory")
    ap.add_argument("--gtf",            required=True, type=Path,
                    help="GTF annotation for HTSeq counting")
    ap.add_argument("--annot_ahrd",     required=True, type=Path,
                    help="AHRD functional annotation TSV (col1=GeneID, col4=Desc, col6=GO)")
    ap.add_argument("--annot_fantasia", required=True, type=Path,
                    help="FANTASIA GO annotation in topGO format (GeneID TAB GO1,GO2,...)")
    ap.add_argument("--output",         required=True, type=Path,
                    help="Output directory (name becomes the run prefix)")

    # DESeq2 design
    ap.add_argument("--deseq_design",   default=None,
                    help=("DESeq2 formula. Column names must match condition columns "
                          "in the sample sheet. See examples below."))
    ap.add_argument("--deseq_stratify", default="",
                    help=("Condition column to stratify by. For each unique level, a "
                          "separate DESeq2 run contrasts --deseq_contrast within that level. "
                          "Example: --deseq_stratify DevelopmentalStage"))
    ap.add_argument("--deseq_contrast", default="",
                    help=("Factor to contrast within each stratum. Defaults to the other "
                          "factor in --deseq_design when omitted. "
                          "Example: --deseq_contrast Temperature"))
    ap.add_argument("--deseq_ref_level", default="",
                    help=("Reference/baseline level of the contrast factor (e.g. 'WT'), "
                          "so log2FoldChange = log2(other_level / this_level). Without "
                          "this, R's factor() defaults to alphabetical order, which "
                          "silently determines the sign of every log2FoldChange -- set "
                          "this explicitly rather than relying on that default. The "
                          "resolved direction is always printed to the log either way. "
                          "Example: --deseq_ref_level WT"))

    # STAR
    ap.add_argument("--star_min_intron", type=int, default=20,
                    help="STAR --alignIntronMin (default: 20)")
    ap.add_argument("--star_max_intron", type=int, default=1_000_000,
                    help="STAR --alignIntronMax (default: 1000000)")

    # fastp
    ap.add_argument("--fastp_min_length", type=int, default=36,
                    help="Minimum read length after trimming (default: 36)")
    ap.add_argument("--fastp_quality",    type=int, default=20,
                    help="Phred quality threshold for fastp (default: 20)")

    # HTSeq
    ap.add_argument("--htseq_stranded", default="no",
                    choices=["yes", "no", "reverse"],
                    help="HTSeq strandedness (default: no)")

    # DESeq2 thresholds / GO
    ap.add_argument("--padj",          type=float, default=0.05,
                    help="Adjusted p-value threshold for DEG calling (default: 0.05)")
    ap.add_argument("--min_counts",    type=int,   default=10,
                    help="Minimum total read count per gene (default: 10)")
    ap.add_argument("--go_pval",       type=float, default=0.01,
                    help="p-value cut-off for topGO reporting (default: 0.01)")
    ap.add_argument("--go_categories", default="BP",
                    help="Comma-separated GO categories: BP, CC, MF (default: BP)")

    # pipeline control
    ap.add_argument("--threads",      type=int, default=4,
                    help="Threads for FastQC, fastp, and STAR (default: 4)")
    ap.add_argument("--skip_module1", action="store_true",
                    help="Skip Module 1 — FastQC pre-trim (global override)")
    ap.add_argument("--skip_module2", action="store_true",
                    help="Skip Module 2 — fastp trimming (global override)")
    ap.add_argument("--skip_module3", action="store_true",
                    help="Skip Module 3 — FastQC post-trim (global override)")
    ap.add_argument("--skip_module4", action="store_true",
                    help="Skip Module 4 — MultiQC")
    ap.add_argument("--skip_module5", action="store_true",
                    help="Skip Module 5 — STAR alignment (global override)")
    ap.add_argument("--skip_module6", action="store_true",
                    help="Skip Module 6 — HTSeq counting (global override)")
    ap.add_argument("--skip_module7", action="store_true",
                    help="Skip Module 7 — DESeq2 + topGO")
    ap.add_argument("--force", action="store_true",
                    help="Rerun all steps even if outputs already exist")
    ap.add_argument("--dry_run", action="store_true",
                    help="Validate inputs and print steps, then exit")
    ap.add_argument("--disable_co2_tracking", action="store_true",
                    help="Disable carbon footprint tracking")
    ap.add_argument("--format", default="pdf",
                    help="Plot format(s): pdf, png, svg — comma-separated (default: pdf)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")

    args = ap.parse_args()

    # resolve absolute paths
    args.sample_sheet   = args.sample_sheet.resolve()
    args.genome_dir     = args.genome_dir.resolve()
    args.gtf            = args.gtf.resolve()
    args.annot_ahrd     = args.annot_ahrd.resolve()
    args.annot_fantasia = args.annot_fantasia.resolve()
    args.output         = args.output.resolve()

    # validate static inputs
    _validate_inputs([
        ("--sample_sheet",   args.sample_sheet),
        ("--genome_dir",     args.genome_dir),
        ("--gtf",            args.gtf),
        ("--annot_ahrd",     args.annot_ahrd),
        ("--annot_fantasia", args.annot_fantasia),
    ])

    samples, condition_cols = parse_sample_sheet(args.sample_sheet)

    # validate per-sample input files
    file_checks = []
    for s in samples:
        sm = s["starting_module"]
        if sm <= 5:
            file_checks += [(f"r1 ({s['sample_name']})", s["r1"]),
                            (f"r2 ({s['sample_name']})", s["r2"])]
        elif sm == 6:
            file_checks.append((f"bam ({s['sample_name']})", s["bam"]))
        else:
            file_checks.append((f"counts ({s['sample_name']})", s["count_file"]))
    _validate_inputs(file_checks)

    # R script must sit alongside this script
    r_script = Path(__file__).parent / "NyaRNAthoDEG_DESeq2.R"
    if not args.skip_module7:
        _validate_inputs([("NyaRNAthoDEG_DESeq2.R", r_script)])
        if not args.deseq_design:
            print("ERROR: --deseq_design is required when running Module 7. "
                  "Use --skip_module7 to skip DESeq2.", file=sys.stderr)
            sys.exit(1)

    # output directories
    run_dir  = args.output
    results  = run_dir / "results"
    workdir  = run_dir / "workdir"
    logs_dir = run_dir / "logs"
    for d in (results, workdir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    prefix = run_dir.name

    # convenience: pre-compute filtered sample lists
    pre_samples  = [s for s in samples if s["starting_module"] == 1]
    post_samples = [s for s in samples if s["starting_module"] <= 3]

    fastqc_pre_dir  = workdir / "mod01_fastqc_pre"
    fastqc_post_dir = workdir / "mod03_fastqc_post"

    # dry run
    if args.dry_run:
        _banner("Dry run — no steps will be executed")
        _log(f"  Sample sheet      : {args.sample_sheet}")
        _log(f"  Samples           : {len(samples)}")
        _log(f"  Condition columns : {', '.join(condition_cols) if condition_cols else 'none'}")
        _log(f"  Genome dir        : {args.genome_dir}")
        _log(f"  GTF               : {args.gtf}")
        _log(f"  Output            : {run_dir}/")
        if args.deseq_design:
            _log(f"  DESeq2 design     : {args.deseq_design}")
        if args.deseq_stratify:
            _log(f"  Stratify by       : {args.deseq_stratify}")
            _log(f"  Contrast          : {args.deseq_contrast or '(inferred)'}")
        _log(f"  Contrast ref level : {args.deseq_ref_level or '(alphabetically first level -- set --deseq_ref_level to control this explicitly)'}")
        _log("  Steps that would run:")
        if not args.skip_module1 and pre_samples:
            _log(f"    [1] FastQC pre-trim  →  workdir/mod01_fastqc_pre/  ({len(pre_samples)} sample(s))")
        if not args.skip_module2:
            _log(f"    [2] fastp            →  workdir/mod02_fastp/")
        if not args.skip_module3 and post_samples:
            _log(f"    [3] FastQC post-trim →  workdir/mod03_fastqc_post/ ({len(post_samples)} sample(s))")
        if not args.skip_module4:
            _log( "    [4] MultiQC          →  results/mod04_multiqc/")
        if not args.skip_module5:
            _log( "    [5] STAR             →  workdir/mod05_star/")
        if not args.skip_module6:
            _log( "    [6] HTSeq            →  workdir/mod06_htseq/")
        if not args.skip_module7:
            _log( "    [7] DESeq2 + topGO   →  results/mod07_deseq2/")
        _log("  Exiting (--dry_run).")
        sys.exit(0)

    # open log file
    global _LOG_FH
    log_path = logs_dir / "Run_NyaRNAthoDEG.log"
    _LOG_FH  = open(log_path, "w")
    sep = "=" * 62
    _LOG_FH.write(f"{sep}\n  NyaRNAthoDEG {VERSION}  —  Run Log\n{sep}\n")
    _LOG_FH.write(f"Date      : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    _LOG_FH.write(f"User      : {getpass.getuser()}\n")
    _LOG_FH.write(f"Server    : {platform.node()}\n")
    _LOG_FH.write(f"OS        : {platform.system()} {platform.release()} "
                  f"({platform.machine()})\n")
    _LOG_FH.write(f"Directory : {os.getcwd()}\n")
    _LOG_FH.write(f"Command   : {' '.join(sys.argv)}\n")
    _LOG_FH.write(f"{sep}\n\n")
    _LOG_FH.flush()

    if args.force:
        _log("--force set: all steps will rerun regardless of existing outputs")
    elif workdir.exists() and any(workdir.iterdir()):
        _log("Existing workdir found — resuming from checkpoints "
             "(use --force to rerun all steps from scratch)")

    # carbon tracker
    _tracker = None
    if args.disable_co2_tracking:
        _log("  Carbon footprint tracking disabled (--disable_co2_tracking)")
    else:
        try:
            try:
                import pkg_resources  # noqa: F401
            except ModuleNotFoundError:
                import types as _t, importlib.metadata as _m, importlib as _il
                from pathlib import Path as _P
                _shim = _t.ModuleType("pkg_resources")
                def _get_dist(name):
                    try:
                        d = _m.distribution(name)
                        d.version = d.metadata["Version"]
                        return d
                    except Exception:
                        return None
                def _resource_filename(pkg, resource):
                    try:
                        mod = _il.import_module(pkg)
                        return str(_P(mod.__file__).parent / resource)
                    except Exception:
                        return resource
                _shim.get_distribution    = _get_dist
                _shim.resource_filename   = _resource_filename
                _shim.DistributionNotFound = Exception
                sys.modules["pkg_resources"] = _shim
            from codecarbon import EmissionsTracker
            _tracker = EmissionsTracker(
                output_dir=str(logs_dir),
                output_file=f"{prefix}.emissions.csv",
                project_name="NyaRNAthoDEG",
                log_level="warning",
            )
            _tracker.start()
            _log("  codecarbon tracker started")
        except ImportError as e:
            _log(f"  codecarbon not installed — carbon tracking skipped ({e})")
        except Exception as e:
            _log(f"  codecarbon failed to start — carbon tracking skipped ({e})")

    t_start = time.monotonic()

    # ── Module 1: FastQC pre-trim ─────────────────────────────────────────────
    if not args.skip_module1 and pre_samples:
        _banner("Module 1 — FastQC (pre-trim)")
        run_fastqc(pre_samples, fastqc_pre_dir, args.threads, args.force, label="pre")

    # ── Module 2: fastp ───────────────────────────────────────────────────────
    if not args.skip_module2:
        _banner("Module 2 — fastp trimming")
        run_fastp(samples, workdir, args.threads,
                  args.fastp_min_length, args.fastp_quality, args.force)
    else:
        fastp_dir = workdir / "mod02_fastp"
        for s in samples:
            sm = s["starting_module"]
            if sm <= 2:
                s["r1_clean"] = fastp_dir / f"{s['sample_name']}_R1.clean.fastq.gz"
                s["r2_clean"] = fastp_dir / f"{s['sample_name']}_R2.clean.fastq.gz"
            elif sm <= 5:
                s["r1_clean"] = s["r1"]
                s["r2_clean"] = s["r2"]

    # ── Module 3: FastQC post-trim ────────────────────────────────────────────
    if not args.skip_module3 and post_samples:
        _banner("Module 3 — FastQC (post-trim)")
        run_fastqc(post_samples, fastqc_post_dir, args.threads, args.force, label="post")

    # ── Module 4: MultiQC ─────────────────────────────────────────────────────
    if not args.skip_module4:
        _banner("Module 4 — MultiQC")
        run_multiqc([fastqc_pre_dir, fastqc_post_dir, workdir / "mod02_fastp"],
                    results / "mod04_multiqc", args.force)

    # ── Module 5: STAR ────────────────────────────────────────────────────────
    if not args.skip_module5:
        star_samples = [s for s in samples if s["starting_module"] <= 5]
        if star_samples:
            _banner("Module 5 — STAR alignment")
            run_star(samples, args.genome_dir, workdir, args.threads,
                     args.star_min_intron, args.star_max_intron, args.force)
    else:
        star_dir = workdir / "mod05_star"
        for s in samples:
            if s["starting_module"] <= 5:
                s["bam"] = star_dir / s["sample_name"] / "Aligned.sortedByCoord.out.bam"

    # ── Module 6: HTSeq ───────────────────────────────────────────────────────
    if not args.skip_module6:
        htseq_samples = [s for s in samples if s["starting_module"] <= 6]
        if htseq_samples:
            _banner("Module 6 — HTSeq counting")
            run_htseq(samples, args.gtf, workdir, args.htseq_stranded, args.force)
    else:
        htseq_dir = workdir / "mod06_htseq"
        for s in samples:
            if s["starting_module"] <= 6:
                s["count_file"] = htseq_dir / f"{s['sample_name']}.htseq.counts"

    # ── Module 7: DESeq2 + topGO ──────────────────────────────────────────────
    if not args.skip_module7:
        _banner("Module 7 — DESeq2 + topGO (R)")
        run_deseq2(samples, condition_cols, workdir, results,
                   args.annot_ahrd, args.annot_fantasia,
                   prefix, args.padj, args.min_counts,
                   args.deseq_design, args.deseq_stratify, args.deseq_contrast,
                   args.deseq_ref_level,
                   args.go_categories, args.go_pval,
                   r_script, args.force)

    # ── resource usage & run summary ──────────────────────────────────────────
    elapsed_s   = time.monotonic() - t_start
    ru          = resource.getrusage(resource.RUSAGE_SELF)
    peak_mem_mb = (ru.ru_maxrss / (1024 * 1024)
                   if platform.system() == "Darwin"
                   else ru.ru_maxrss / 1024)

    emissions_kg = None
    if _tracker is not None:
        try:
            emissions_kg = _tracker.stop()
        except Exception:
            pass

    summary = {
        "date":    datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "version": VERSION,
        "input_sample_sheet":   str(args.sample_sheet),
        "input_genome_dir":     str(args.genome_dir),
        "input_gtf":            str(args.gtf),
        "input_annot_ahrd":     str(args.annot_ahrd),
        "input_annot_fantasia": str(args.annot_fantasia),
        "n_samples":        len(samples),
        "samples":          [s["sample_name"] for s in samples],
        "condition_cols":   condition_cols,
        "parameters": {
            "threads":           args.threads,
            "fastp_min_length":  args.fastp_min_length,
            "fastp_quality":     args.fastp_quality,
            "htseq_stranded":    args.htseq_stranded,
            "star_min_intron":   args.star_min_intron,
            "star_max_intron":   args.star_max_intron,
            "deseq_design":      args.deseq_design,
            "deseq_stratify":    args.deseq_stratify,
            "deseq_contrast":    args.deseq_contrast,
            "deseq_ref_level":   args.deseq_ref_level,
            "padj":              args.padj,
            "min_counts":        args.min_counts,
            "go_pval":           args.go_pval,
            "go_categories":     args.go_categories,
        },
        "resource_usage": {
            "wall_clock_s":       round(elapsed_s, 1),
            "peak_mem_mb":        round(peak_mem_mb, 1),
            "emissions_kg_CO2eq": emissions_kg,
        },
    }

    summary_path = results / f"{prefix}.run_summary.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")

    _banner("NyaRNAthoDEG complete")
    _log(f"  Samples processed : {len(samples)}")
    _log(f"  Wall time         : {elapsed_s:.1f} s")
    _log(f"  Peak memory       : {peak_mem_mb:.1f} MB")
    if emissions_kg is not None:
        _log(f"  CO2 emissions     : {emissions_kg:.6f} kg CO2eq")
    _log(f"  Results           : {results}/")
    _log(f"  Run log           : {log_path}")

    if _LOG_FH is not None:
        _LOG_FH.close()


if __name__ == "__main__":
    main()
