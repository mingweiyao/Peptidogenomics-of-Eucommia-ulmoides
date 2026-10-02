"""R1-Minor2: NCP sequence conservation screen (exact k-mer matching, k=8).
References: Arabidopsis/poplar/grape proteomes + the annotated Eucommia proteome (Eu_protein.faa)."""
import pandas as pd, numpy as np, gzip, os
import os as _os
import argparse
import os

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--eu_faa", default=r"./data/figure_middleware/Eu_protein.faa", help="annotated Eucommia proteome fasta")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))

def load_fa(p):
    seqs = []
    cur = None
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None: seqs.append(cur)
                cur = ""
            else:
                cur += line.strip()
    if cur is not None: seqs.append(cur)
    return seqs

def kmers(prots, k=8):
    km = set()
    for pr in prots:
        if len(pr) < k: continue
        for i in range(len(pr) - k + 1):
            km.add(pr[i:i+k])
    return km

HERE = os.path.dirname(os.path.abspath(__file__))
s4 = pd.read_csv(OUTROOT + "/_base_data_chain/supplementary_table4_with_recomputed_props.csv", low_memory=False)
s4["Pep_ID"] = s4["Chr"].astype(str).str.replace("GWHBISF","",regex=False).map(lambda v: str(int(v)))+"_"+s4["Start"].astype(str)+"_"+s4["End"].astype(str)+"_"+s4["Strand"]
u = s4.drop_duplicates(subset=["Pep_ID"]).copy()
ncp_seqs = u["Sequence"].astype(str).tolist()
print("NCP unique:", len(u))

species = {
 "Arabidopsis_thaliana": "Arabidopsis_thaliana.fa.gz",
 "Populus_trichocarpa": "Populus_trichocarpa.fa.gz",
 "Vitis_vinifera": "Vitis_vinifera.fa.gz",
}
eu_faa = args.eu_faa
OUT = OUTROOT + "/18_revision_supplements"
os.makedirs(OUT, exist_ok=True)
if os.path.exists(eu_faa):
    species["Eucommia_annotated"] = eu_faa

res = {}
for name, f in species.items():
    p = f if os.path.isabs(f) else os.path.join(HERE, f)
    prots = load_fa(p)
    km = kmers(prots, 8)
    hits = sum(1 for s in ncp_seqs if s in km or any(s[i:i+8] in km for i in range(len(s)-7)))
    res[name] = (len(prots), hits)
    print(f"{name}: proteins {len(prots)} | 8-mer hits in NCPs {hits}/{len(ncp_seqs)} ({100*hits/len(ncp_seqs):.2f}%)")

out = pd.DataFrame([(k, v[0], v[1], round(100*v[1]/len(u),2)) for k, v in res.items()],
                   columns=["Species","Proteins","NCPs_with_8mer_match","Percent"])
out.to_csv(os.path.join(OUT, "conservation_kmer_screen.csv"), index=False)
print("Wrote conservation_kmer_screen.csv")
