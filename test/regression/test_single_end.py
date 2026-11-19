#!/usr/bin/env python3
"""Single-end support (added 2026-11-19).

Many public abiotic-stress RNA-seq experiments are single-end; before this the
sample sheet loader rejected them outright.

Run:  python3 test/regression/test_single_end.py
"""
import sys, tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import NyaRNAthoDEG as N  # noqa: E402

FAILS = 0
def ok(name, cond, info=""):
    global FAILS
    print(f"{name:<58} {'PASS' if cond else 'FAIL'}{'  (' + info + ')' if info else ''}")
    if not cond: FAILS += 1

def sheet(tmp, rows):
    p = Path(tmp) / "s.tsv"
    p.write_text("SampleName\tStartingModule\tFilePath\tCondition\n" +
                 "".join(f"{n}\t{sm}\t{fp}\t{c}\n" for n, sm, fp, c in rows))
    return p

with tempfile.TemporaryDirectory() as tmp:
    for f in ("a_R1.fastq.gz", "a_R2.fastq.gz", "b.fastq.gz"):
        (Path(tmp) / f).write_bytes(b"")

    # single-end: one FASTQ path
    s_se, _ = N.parse_sample_sheet(sheet(tmp, [("b", 1, f"{tmp}/b.fastq.gz", "stress")]))
    ok("single-end sheet loads", len(s_se) == 1)
    ok("single-end sample is marked SE", s_se[0].get("layout") == "SE", str(s_se[0].get("layout")))
    ok("single-end sample has no r2", s_se[0].get("r2") is None)

    # paired-end still works
    s_pe, _ = N.parse_sample_sheet(sheet(tmp, [("a", 1, f"{tmp}/a_R1.fastq.gz,{tmp}/a_R2.fastq.gz", "ctrl")]))
    ok("paired-end sample is marked PE", s_pe[0].get("layout") == "PE", str(s_pe[0].get("layout")))
    ok("paired-end keeps r2", s_pe[0].get("r2") is not None)

    # command builders
    se, pe = dict(s_se[0]), dict(s_pe[0])
    se["r1_clean"], se["r2_clean"] = se["r1"], None
    pe["r1_clean"], pe["r2_clean"] = pe["r1"], pe["r2"]

    c = N.fastp_cmd("fastp", se, Path(tmp), 4, 36, 20)
    ok("fastp SE: no --in2/--out2", "--in2" not in c and "--out2" not in c)
    ok("fastp SE: no --detect_adapter_for_pe", "--detect_adapter_for_pe" not in c)
    c = N.fastp_cmd("fastp", pe, Path(tmp), 4, 36, 20)
    ok("fastp PE: keeps --in2 and the PE adapter flag",
       "--in2" in c and "--detect_adapter_for_pe" in c)

    c = N.star_cmd("STAR", se, Path(tmp), Path(tmp), 4, 20, 100000)
    i = c.index("--readFilesIn")
    ok("STAR SE: one reads file", c[i+1].endswith("b.fastq.gz") and c[i+2].startswith("--"))
    c = N.star_cmd("STAR", pe, Path(tmp), Path(tmp), 4, 20, 100000)
    i = c.index("--readFilesIn")
    ok("STAR PE: two reads files", c[i+1].endswith("_R1.fastq.gz") and c[i+2].endswith("_R2.fastq.gz"))

print(f"\n{'ALL PASS' if FAILS == 0 else 'FAILURES'} ({FAILS} failing)")
sys.exit(0 if FAILS == 0 else 1)
