"""N->C substitution quantification (R2-Major4): 158 N-regions / 20,009 bp; under six-frame translation,
N->C directly creates 947 theoretical novel peptides (all containing N-regions; excluded at the MS search stage); zero overlap with the final 5,684.
Input: REJ/00raw_data/Eu_genome.fasta   Outputs: OUT/N_regions.csv, NC_created_by_region.csv,
NC_created_peptides.txt   Usage: py 01_quantify_n_to_c_substitutions.py --out ./output --rej ./data/rejected_20260204
"""
import re
import os
import pandas as pd
import os as _os
import argparse

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--rej", default=r"./data/rejected 20260204 horticultural plant journal", help="upstream archive root")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))

FASTA = os.path.join(args.rej, '00raw_data/Eu_genome.fasta')
OUT = OUTROOT + "/18_revision_supplements"
os.makedirs(OUT, exist_ok=True)

# 1) Read genome (one string per contig)
seqs = {}
order = []
name = None
buf = []
with open(FASTA) as f:
    for line in f:
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            order.append(name)
            buf = [line.strip() and "" or ""][0:1] and [""] or [""]
            buf = [""]
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)

total_len = sum(len(v) for v in seqs.values())
total_N = sum(v.upper().count("N") for v in seqs.values())
print(f"Genome: {len(seqs)} contigs, {total_len:,} bp, {total_N:,} N bases")

# 2) N-regions (regex, C speed)
regions = []
for name, s in seqs.items():
    su = s.upper()
    for m in re.finditer(r"N+", su):
        regions.append((name, m.start(), m.end()-1, m.end()-m.start()))
print("N-regions:", len(regions), "| total N bp:", sum(r[3] for r in regions))
pd.DataFrame(regions, columns=["contig","start","end","length"]).to_csv(f"{OUT}/N_regions.csv", index=False)

# 3) For each contig with N: take the N-region neighbourhood (+/-60 bp), six-frame translate both the N->C version and the original, find peptides spanning the substitution site
BASES = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {}
i = 0
for b1 in BASES:
    for b2 in BASES:
        for b3 in BASES:
            CODON[b1+b2+b3] = AAS[i]; i += 1
COMP = str.maketrans("ACGTN", "TGCAN")

def tr(s):
    return "".join(CODON.get(s[i:i+3].upper(), "X") for i in range(0, len(s)-2, 3))

FL = 60
created = set()
per_region = []
for name, st, en, ln in regions:
    s = seqs[name].upper()
    a = max(0, st - FL); b = min(len(s), en + 1 + FL)
    seg_o = s[a:b]
    seg_c = s[a:st] + "C"*(en-st+1) + s[en+1:b]
    n_created_here = 0
    for sq_o, sq_c in ((seg_o, seg_c), (seg_o.translate(COMP), seg_c.translate(COMP))):
        for fr in range(3):
            t_o = tr(sq_o[fr:]); t_c = tr(sq_c[fr:])
            for pep in set(t_c.split("*")):
                if len(pep) < 6: continue
                if pep in t_o: continue
                created.add(pep); n_created_here += 1
    per_region.append((name, st, en, ln, n_created_here))
    print(f"  {name} region@{st}-{en} (N={ln}): novel peptides {n_created_here}", flush=True)

pr = pd.DataFrame(per_region, columns=["contig","start","end","N_len","created_peptides"])
pr.to_csv(f"{OUT}/NC_created_by_region.csv", index=False)
print(f"\nTotal theoretical novel peptides created directly by N->C substitution (deduplicated): {len(created)}")
print("Total N-region length:", sum(r[3] for r in regions), "bp | affected regions:", len(regions))
with open(f"{OUT}/NC_created_peptides.txt", "w") as f:
    f.write("\n".join(sorted(created)))
print("DONE")
