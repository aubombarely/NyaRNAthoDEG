#!/usr/bin/env python3
"""Regression tests for the HTSeq assignment-rate check (added 2026-11-19).

The check exists to catch a silent failure: htseq-count run with the wrong
--stranded value still produces a counts file and DESeq2 still returns
"significant" genes, but most reads were never assigned to a feature.

Run:  python3 test/regression/test_count_qc.py
"""
import sys, tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from NyaRNAthoDEG import count_assignment_stats  # noqa: E402

FAILS = 0

def ok(name, cond, info=""):
    global FAILS
    print(f"{name:<58} {'PASS' if cond else 'FAIL'}{'  (' + info + ')' if info else ''}")
    if not cond:
        FAILS += 1

def write_counts(tmp, genes, specials):
    p = Path(tmp) / "x.htseq.counts"
    lines = [f"g{i}\t{c}" for i, c in enumerate(genes)]
    lines += [f"{k}\t{v}" for k, v in specials.items()]
    p.write_text("\n".join(lines) + "\n")
    return p

with tempfile.TemporaryDirectory() as tmp:
    # a healthy sample: 9,000 assigned of 10,000
    p = write_counts(tmp, [1000] * 9, {"__no_feature": 600, "__ambiguous": 100,
                                       "__too_low_aQual": 0, "__not_aligned": 200,
                                       "__alignment_not_unique": 100})
    st = count_assignment_stats(p)
    ok("assigned reads counted correctly", st["assigned"] == 9000, str(st["assigned"]))
    ok("total reads counted correctly", st["total"] == 10000, str(st["total"]))
    ok("assignment rate computed", abs(st["rate"] - 0.9) < 1e-9, f"{st['rate']:.3f}")

    # the strandedness-mismatch signature: almost everything in __no_feature
    p = write_counts(tmp, [10] * 9, {"__no_feature": 9000, "__ambiguous": 10,
                                     "__too_low_aQual": 0, "__not_aligned": 800,
                                     "__alignment_not_unique": 100})
    st = count_assignment_stats(p)
    ok("mismatch signature gives a low rate", st["rate"] < 0.05, f"{st['rate']:.3f}")
    ok("dominant loss category identified",
       st["dominant_loss"] == "__no_feature", st["dominant_loss"])

    # a counts file with no summary lines at all must not crash
    p = write_counts(tmp, [5] * 4, {})
    st = count_assignment_stats(p)
    ok("file without __ lines is handled", st["assigned"] == 20 and st["rate"] is None,
       f"assigned={st['assigned']} rate={st['rate']}")

# the bundled fixtures must look healthy
fx = Path(__file__).resolve().parents[1] / "htseq_counts" / "Stage1_L_rep1.htseq.counts"
if fx.exists():
    st = count_assignment_stats(fx)
    ok("bundled test fixture has a sane assignment rate",
       st["rate"] is not None and st["rate"] > 0.5, f"{st['rate']:.3f}")

print(f"\n{'ALL PASS' if FAILS == 0 else 'FAILURES'} ({FAILS} failing)")
sys.exit(0 if FAILS == 0 else 1)
