# -*- coding: utf-8 -*-
"""SHA-256 checksums of key input files (R2-Minor12): genome / GFF / six-frame peptide DB / aa pipeline output.

Verification: every file must match the archived reference_file_checksums.csv (response-letter reference).
Inputs: four files under REJ/00raw_data/ (~1.1 GB, 1-2 min to read)
Output: output/18_revision_supplements/reference_file_checksums_recomputed.csv + on-screen comparison
Usage: py 02_verify_input_checksums.py --out ./output --rej ./data/rejected_20260204
"""
import hashlib
import pandas as pd
import os
import os as _os
import argparse

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--rej", default=r"./data/rejected 20260204 horticultural plant journal", help="upstream archive root")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))

REJ = os.path.abspath(os.path.expanduser(args.rej))
GLM = "./data/response_letter_reference"
OUT = OUTROOT + "/18_revision_supplements"
os.makedirs(OUT, exist_ok=True)

FILES = [
    os.path.join(args.rej, '00raw_data/Eu_genome.fasta'),
    os.path.join(args.rej, '00raw_data/GWHBISF00000000.gff'),
    REJ + "/00raw_data/Eu_peptide_database_customized_5.fa",
    os.path.join(args.rej, '00raw_data/peptide_analysis_results.xlsx'),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


rows = []
for p in FILES:
    digest = sha256(p)
    rows.append({"File": os.path.basename(p), "Bytes": os.path.getsize(p), "SHA256": digest})
    print(f"{os.path.basename(p)}  {digest}")
res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/reference_file_checksums_recomputed.csv", index=False)

ref = pd.read_csv(GLM + "/reference_file_checksums.csv")
m = res.merge(ref, on="File", suffixes=("_now", "_ref"))
m["match"] = m["SHA256_now"] == m["SHA256_ref"]
print("\nComparison with the archived reference:")
print(m[["File", "match"]].to_string(index=False))
print("PASS" if m["match"].all() else "MISMATCH found! Check the data version")
