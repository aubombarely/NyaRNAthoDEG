#!/usr/bin/env python3
"""
make_test_data.py — generate synthetic test data for NyaRNAthoDEG
Run once before the quicktest command.

Outputs (all inside this script's directory):
  test_samples.tsv          sample sheet  (StartingModule=7, 8 samples)
  htseq_counts/             synthetic HTSeq count files
  test_ahrd.txt             AHRD-format annotation  (500 genes)
  test_fantasia.txt         FANTASIA topGO-format annotation
  test_annotation.gtf       minimal dummy GTF (CLI placeholder only)
  test_genome_dir/          empty directory  (CLI placeholder only)

Experiment design
  DevelopmentalStage : Stage1, Stage2
  Temperature        : L (low), H (high)
  Replicates         : 2 per condition  →  8 samples total

Baked-in signal (detectable by DESeq2):
  genes 1–100    up in Temperature H  (~4× fold change)
  genes 101–200  down in Temperature H  (~4× fold change)
  genes 201–300  up in Stage2  (~3× fold change)
  genes 301–400  down in Stage2  (~3× fold change)
  genes 401–500  no effect (background)
"""

import csv, math, random
from pathlib import Path

random.seed(42)

HERE     = Path(__file__).parent.resolve()
N_GENES  = 500
GENE_IDS = [f"SMELv4.1_{i:05d}" for i in range(1, N_GENES + 1)]

# Baseline mean expression per gene (Poisson-ish, log-normal)
BASELINE = [max(5, int(random.lognormvariate(4.0, 1.5))) for _ in range(N_GENES)]

# Small set of GO terms used for annotation
GO_TERMS = [
    "GO:0006096", "GO:0006412", "GO:0055114", "GO:0006260",
    "GO:0007049", "GO:0006351", "GO:0006355", "GO:0005524",
    "GO:0003723", "GO:0016740", "GO:0003824", "GO:0005634",
    "GO:0005737", "GO:0005739",
]


# ── helpers ───────────────────────────────────────────────────────────────────

def _sim_count(gene_idx: int, stage: str, temperature: str) -> int:
    """Simulate a single read count with added Gaussian noise."""
    base = BASELINE[gene_idx]

    # Temperature effect
    if gene_idx < 100:
        t_eff = 4.0 if temperature == "H" else 1.0
    elif gene_idx < 200:
        t_eff = 0.25 if temperature == "H" else 1.0
    else:
        t_eff = 1.0

    # Stage effect
    if 200 <= gene_idx < 300:
        s_eff = 3.0 if stage == "Stage2" else 1.0
    elif 300 <= gene_idx < 400:
        s_eff = 0.33 if stage == "Stage2" else 1.0
    else:
        s_eff = 1.0

    mu    = base * t_eff * s_eff
    noise = random.gauss(0, math.sqrt(mu) + 1)
    return max(0, round(mu + noise))


def _go_terms_for_gene(gene_idx: int) -> list:
    """Assign 1–3 GO terms semi-randomly per gene."""
    rng = random.Random(gene_idx * 7919)  # deterministic per gene
    n   = rng.randint(1, 3)
    return rng.sample(GO_TERMS, n)


# ── 1. HTSeq count files ──────────────────────────────────────────────────────

SAMPLES = [
    ("Stage1_L_rep1", "Stage1", "L", "1"),
    ("Stage1_L_rep2", "Stage1", "L", "2"),
    ("Stage1_H_rep1", "Stage1", "H", "1"),
    ("Stage1_H_rep2", "Stage1", "H", "2"),
    ("Stage2_L_rep1", "Stage2", "L", "1"),
    ("Stage2_L_rep2", "Stage2", "L", "2"),
    ("Stage2_H_rep1", "Stage2", "H", "1"),
    ("Stage2_H_rep2", "Stage2", "H", "2"),
]

counts_dir = HERE / "htseq_counts"
counts_dir.mkdir(exist_ok=True)

for name, stage, temp, rep in SAMPLES:
    rng_rep = random.Random(hash((name, stage, temp, rep)) & 0xFFFFFFFF)
    out = counts_dir / f"{name}.htseq.counts"
    with open(out, "w") as fh:
        for idx, gid in enumerate(GENE_IDS):
            # Add per-replicate noise
            c = _sim_count(idx, stage, temp)
            c = max(0, round(c * rng_rep.uniform(0.85, 1.15)))
            fh.write(f"{gid}\t{c}\n")
        # HTSeq summary lines (required; DESeq2 ignores them)
        fh.write("__no_feature\t312\n")
        fh.write("__ambiguous\t87\n")
        fh.write("__too_low_aQual\t0\n")
        fh.write("__not_aligned\t1024\n")
        fh.write("__alignment_not_unique\t256\n")
    print(f"  Written: {out.name}")

print(f"HTSeq count files → {counts_dir}/\n")


# ── 2. Sample sheet ───────────────────────────────────────────────────────────

sample_sheet = HERE / "test_samples.tsv"
with open(sample_sheet, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["SampleName", "StartingModule", "FilePath",
                     "DevelopmentalStage", "Temperature", "Replicate"])
    for name, stage, temp, rep in SAMPLES:
        fp = counts_dir / f"{name}.htseq.counts"
        writer.writerow([name, 7, str(fp), stage, temp, rep])

print(f"Sample sheet      → {sample_sheet.name}\n")


# ── 3. AHRD annotation ────────────────────────────────────────────────────────
# Columns: GeneID, Blast-Hit-Accession, AHRD-Quality-Code,
#          Human-Readable-Description, InterPro-ID, Gene-Ontology-ID

DESCRIPTIONS = [
    "photosystem II protein", "ribosomal protein L", "heat shock protein",
    "WRKY transcription factor", "cytochrome P450", "F-box protein",
    "serine/threonine kinase", "ubiquitin ligase", "peroxidase",
    "alcohol dehydrogenase", "hypothetical protein", "unknown function",
]

ahrd_file = HERE / "test_ahrd.txt"
with open(ahrd_file, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["Protein-Accession", "Blast-Hit-Accession",
                     "AHRD-Quality-Code", "Human-Readable-Description",
                     "InterPro-ID", "Gene-Ontology-ID"])
    for idx, gid in enumerate(GENE_IDS):
        rng_g  = random.Random(idx * 1013)
        desc   = rng_g.choice(DESCRIPTIONS) + f" {idx + 1}"
        go_str = ", ".join(_go_terms_for_gene(idx))
        writer.writerow([gid, f"sp|P{idx:05d}|TEST", "*", desc, "", go_str])

print(f"AHRD annotation   → {ahrd_file.name}\n")


# ── 4. FANTASIA annotation (topGO format) ─────────────────────────────────────
# Format: GeneID<TAB>GO1,GO2,...   (no header)

fantasia_file = HERE / "test_fantasia.txt"
with open(fantasia_file, "w") as fh:
    for idx, gid in enumerate(GENE_IDS):
        go_list = _go_terms_for_gene(idx)
        fh.write(f"{gid}\t{', '.join(go_list)}\n")

print(f"FANTASIA annot.   → {fantasia_file.name}\n")


# ── 5. Dummy GTF (CLI placeholder) ───────────────────────────────────────────

gtf_file = HERE / "test_annotation.gtf"
with open(gtf_file, "w") as fh:
    fh.write('##description: minimal dummy GTF for NyaRNAthoDEG test\n')
    for idx, gid in enumerate(GENE_IDS[:10]):
        start = idx * 2000 + 1
        end   = start + 1800
        fh.write(f'Chr1\t.\tgene\t{start}\t{end}\t.\t+\t.\t'
                 f'gene_id "{gid}"; gene_name "{gid}";\n')

print(f"Dummy GTF         → {gtf_file.name}  (placeholder, not used at StartingModule=7)\n")


# ── 6. Dummy genome dir (CLI placeholder) ────────────────────────────────────

genome_dir = HERE / "test_genome_dir"
genome_dir.mkdir(exist_ok=True)
(genome_dir / "genomeParameters.txt").write_text(
    "# dummy STAR genome dir for NyaRNAthoDEG test\n"
)

print(f"Dummy genome dir  → {genome_dir.name}/  (placeholder, not used at StartingModule=7)\n")


# ── Summary ───────────────────────────────────────────────────────────────────

print("=" * 55)
print("  Test data generated successfully.")
print(f"  Samples    : {len(SAMPLES)}")
print(f"  Genes      : {N_GENES}")
print(f"  Signal     : genes 1–200 Temperature DE")
print(f"               genes 201–400 Stage DE")
print()
print("  Quicktest command:")
print()
print("    python3 scripts/NyaRNAthoDEG.py \\")
print("        --sample_sheet   test/test_samples.tsv \\")
print("        --genome_dir     test/test_genome_dir \\")
print("        --gtf            test/test_annotation.gtf \\")
print("        --annot_ahrd     test/test_ahrd.txt \\")
print("        --annot_fantasia test/test_fantasia.txt \\")
print("        --output         test/test_run \\")
print('        --deseq_design   "~ DevelopmentalStage + Temperature" \\')
print("        --deseq_stratify DevelopmentalStage \\")
print("        --skip_module1 --skip_module2 --skip_module3 \\")
print("        --skip_module4 --skip_module5 --skip_module6")
print("=" * 55)
