# -*- coding: utf-8 -*-
"""Cache every sheet of peptide_analysis_results.xlsx (aa pipeline output) as pkl files.

To use your own re-exported peptide_analysis_results.xlsx:
    py 01_prepare_source_tables.py --peptide_xlsx "your file" --out "output root"
    (also delete all old pkls in <out>/_source_cache/ - the script skips caches that already exist,
      so stale files would silently keep being used! Then rerun scripts 02-05 and downstream in order.)
   The script verifies the file SHA-256: identical to the archived version = same file, downstream
   numbers reproduce exactly; different = a re-export, row counts may deviate from the manuscript targets (warning printed).

Note: Eu_sp_finally.xlsx (the hand-curated original) is not part of this chain - it is an archived
    input, not a reproduction product. Its three final sheets (CP 266,830 / NCP 343,601 / GFF 336,380) are
    exactly what 02_build_cp_ncp_tables.py and 01_clean_gff_coordinates.py re-derive, so it is neither read nor written here.

Inputs (read-only):
  peptide_analysis_results.xlsx (set via --peptide_xlsx; archived authoritative version =
      SHA256 = 16a365f9bf02056e392bbc8d8089d10058c39cb8d29cbb10a1ecbce34ebc781d，
      97,764,606 bytes; 4 sheets: Annotated_Peptides/Annotated_Proteins/Genomic_Features/removed_n)
Outputs:
  <out>/_source_cache/par_<sheet>.pkl (set via --out, default ./output)
Usage: py 01_prepare_source_tables.py --peptide_xlsx ./data/aa_pipeline/peptide_analysis_results.xlsx --out ./output
"""
import hashlib
import pandas as pd
from pathlib import Path
import os as _os
import argparse
import os

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--peptide_xlsx", default=r"./data/aa_pipeline/peptide_analysis_results.xlsx", help="aa pipeline result xlsx (input)")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))


AUTH_SHA256 = "16a365f9bf02056e392bbc8d8089d10058c39cb8d29cbb10a1ecbce34ebc781d"
AUTH_BYTES = 97764606

OUT = Path(OUTROOT + "/_source_cache")
OUT.mkdir(parents=True, exist_ok=True)

peptide_xlsx = args.peptide_xlsx

# ---- File provenance check (SHA-256 and byte size) ----
h = hashlib.sha256()
with open(peptide_xlsx, "rb") as f:
    for chunk in iter(lambda: f.read(1 << 20), b""):
        h.update(chunk)
digest = h.hexdigest()
size = Path(peptide_xlsx).stat().st_size
print(f"Using: {peptide_xlsx}")
print(f"  Bytes: {size:,} (archived authoritative version {AUTH_BYTES:,})")
if digest == AUTH_SHA256:
    print("  SHA-256 matches the archived version -> same file, downstream numbers reproduce exactly [OK]")
else:
    print(f"  SHA256 = {digest}")
    print("  [WARN] Does not match the archived version: this is a re-export. Downstream row counts/numbers")
    print("    Confirm the export used the same aa pipeline version and settings.")

xl = pd.ExcelFile(peptide_xlsx)
for sheet in xl.sheet_names:
    out = OUT / f"par_{sheet}.pkl"
    if out.exists():
        print(f"Skipping existing cache: {out.name}", flush=True)
        continue
    df = pd.read_excel(xl, sheet_name=sheet)
    df.to_pickle(out)
    print(f"{Path(peptide_xlsx).name} [{sheet}]: {df.shape[0]} rows x {df.shape[1]} cols", flush=True)
print(f"DONE → {OUT}")
