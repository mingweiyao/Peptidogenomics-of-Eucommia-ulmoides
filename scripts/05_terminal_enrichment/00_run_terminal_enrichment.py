# -*- coding: utf-8 -*-
"""Terminal-enrichment suite one-stop: runs the baseline permutation -> four constrained null models -> plot.

Usage: py 00_run_terminal_enrichment.py [--out output root] [--rej archive root] [--gen13 figure middleware root]
Outputs (<--out>/SupFig5_terminal_enrichment/):
  Fig3C_terminal_enrichment_permutation.csv              baseline (17/17 significant)
  Fig3C_terminal_enrichment_permutation_NCrestricted.csv non-CDS-restricted
  Fig3C_terminal_enrichment_GCexact.csv                  GC-exact
  Fig3C_terminal_enrichment_GCmatched.csv                GC-stratified-matched
  Fig3C_terminal_enrichment_ExprRestricted.csv           expression-restricted (1,940 vs 1,652.8)
  SuppFig5_terminal_enrichment.png/pdf                   figure (A baseline / B expression-restricted)

Note: this is a runner that invokes the five model scripts and the plot script in order;
      see the header comments of the 01-05 null-model scripts for implementation details.
"""
import argparse
import os
import runpy
import sys

_ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--rej", default=r"./data/rejected 20260204 horticultural plant journal",
                 help="upstream archive root")
_ap.add_argument("--gen13", default=r"./data/figure_middleware",
                 help="figure middleware root")
args = _ap.parse_args()

OUT = os.path.abspath(os.path.expanduser(args.out))
HERE = os.path.dirname(os.path.abspath(__file__))

STEPS = [
    "01_run_baseline_null.py",
    "02_run_noncoding_restricted_null.py",
    "03_run_gc_exact_null.py",
    "04_run_gc_matched_null.py",
    "05_run_expression_restricted_null.py",
    "06_plot_terminal_enrichment.py",
]
for s in STEPS:
    print(f"\n========== Running {s} ==========", flush=True)
    sys.argv = [s, "--out", OUT, "--rej", os.path.abspath(os.path.expanduser(args.rej)),
                "--gen13", os.path.abspath(os.path.expanduser(args.gen13))]
    runpy.run_path(os.path.join(HERE, s), run_name="__main__")
print("\nAll done ->", os.path.join(OUT, "SupFig5_terminal_enrichment"))
