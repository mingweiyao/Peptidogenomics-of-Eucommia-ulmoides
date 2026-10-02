# -*- coding: utf-8 -*-
"""Stress WGCNA archive reconstruction and rule back-inference (no WGCNA rerun; bit-exact reconciliation against the rejected-archive products only).

Five steps:
 1) dt_st_count_fpkm.txt -> t() -> MAD>0 universe (expect 21,608, incl. 267 NCPs)
 2) hub_genes/genes_<color>_trait_CK.txt (9 files = full module member lists, original pipeline Step 8 products)
    + grey (universe remainder) -> rebuild gene_module_info.txt
 3) Back-infer the veen.xlsx list-generation rules:
      st_module(181) = all NCPs of modules {brown,green,magenta,pink,red,turquoise,yellow}
      dt_module(231) = all NCPs of modules {black,blue,brown,pink,red,turquoise,yellow}
                       + EuNCP_trans_2492 (the only magenta member in the DT list; historical quirk of the source table)
    cross-checked against the source table module_trait_associate_pep.xlsx (first-row IDs swallowed by Excel as headers in the DT/ST sheets)
 4) Cross-check: novel union of the 9 high-TOM pair files in module_specific_results = 245
 5) Candidate chain: dt 231 intersect nominal DE 51 -> 35; st 181 intersect nominal DE 35 -> 27

Inputs (--rej root; must mirror the submitted archive layout):
  02figure/figure6/WGCNA/dt_st_count_fpkm.txt
  02figure/figure6/WGCNA/hub_genes/genes_<color>_trait_CK.txt                  ×9
  02figure/figure6/WGCNA/module_specific_results/<color>_high_conf_pairs.txt   ×9
  02figure/figure6/veen.xlsx
  02figure/figure6/module_trait_associate_pep.xlsx
Outputs (<out>/S5_stress_candidates/):
  gene_module_info_reconstructed.txt     21,608 rows of Gene/Module/Type
  stress_wgcna_reconstruction_report.csv per-module reconciliation table
"""
import argparse
import os

import pandas as pd

MODS = ['black', 'blue', 'brown', 'green', 'magenta', 'pink', 'red', 'turquoise', 'yellow']
DT_MODS = ['black', 'blue', 'brown', 'pink', 'red', 'turquoise', 'yellow']
ST_MODS = ['brown', 'green', 'magenta', 'pink', 'red', 'turquoise', 'yellow']
DT_EXTRA = 'EuNCP_trans_2492'

_ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--rej", default=r"./data/rejected 20260204 horticultural plant journal",
                 help="rejected archive root (same layout suffices)")
args = _ap.parse_args()
REJ = os.path.abspath(os.path.expanduser(args.rej))
OUT = os.path.join(os.path.abspath(os.path.expanduser(args.out)), "S5_stress_candidates")
os.makedirs(OUT, exist_ok=True)
W = os.path.join(REJ, "02figure", "figure6", "WGCNA")
F6 = os.path.join(REJ, "02figure", "figure6")

def ok(flag):
    return "✅" if flag else "❌"

# ---------- ① universe ----------
d = pd.read_csv(os.path.join(W, "dt_st_count_fpkm.txt"), sep="\t", index_col=0)
X = d.T
mad = (X - X.median()).abs().median()
universe = list(X.columns[mad > 0])
ncp_uni = [g for g in universe if g.startswith("EuNCP_trans")]
print(f"1) MAD>0 universe = {len(universe)} (expect 21,608) {ok(len(universe) == 21608)}"
      f" | NCP = {len(ncp_uni)} (expect 267) {ok(len(ncp_uni) == 267)} | samples = {X.shape[0]} (expect 14)")

# ---------- 2) Rebuild gene_module_info ----------
gene2mod = {}
sizes = {}
for m in MODS:
    genes = open(os.path.join(W, "hub_genes", f"genes_{m}_trait_CK.txt")).read().split()
    sizes[m] = len(genes)
    for g in genes:
        gene2mod[g] = m
member = set(gene2mod)
grey = [g for g in universe if g not in member]
print(f"2) 9-module member union = {len(member)} | grey = {len(grey)} | total {len(member) + len(grey)}"
      f" {ok(len(member) + len(grey) == 21608)} | all members within universe {ok(member <= set(universe))}"
      f" | module sizes {sizes}")

rows = [(g, gene2mod.get(g, "grey"), "annotated" if g.startswith("evm") else "novel") for g in universe]
info = pd.DataFrame(rows, columns=["Gene", "Module", "Type"])
info_path = os.path.join(OUT, "gene_module_info_reconstructed.txt")
info.to_csv(info_path, sep="\t", index=False)
print(f"   rebuilt gene_module_info -> {info_path} ({len(info)} rows)")

ncp_mod = {}
for g in ncp_uni:
    ncp_mod.setdefault(gene2mod.get(g, "grey"), set()).add(g)
grey_ncp = sorted(ncp_mod.get("grey", set()))
print(f"   NCPs in grey (2): {grey_ncp}")

# ---------- 3) veen list rule back-inference ----------
veen = pd.read_excel(os.path.join(F6, "veen.xlsx"))
dtm = set(veen["dt_module"].dropna().astype(str))
stm = set(veen["st_module"].dropna().astype(str))
dt_deseq = set(veen["dt_deseq"].dropna().astype(str))
st_deseq = set(veen["st_deseq"].dropna().astype(str))

rule_dt = set().union(*(ncp_mod[m] for m in DT_MODS)) | {DT_EXTRA}
rule_st = set().union(*(ncp_mod[m] for m in ST_MODS))
print(f"3) Rule back-inference: DT(7 modules+{DT_EXTRA}) = {len(rule_dt)} == veen dt_module(231) {ok(rule_dt == dtm)}"
      f" | ST(7 modules) = {len(rule_st)} == veen st_module(181) {ok(rule_st == stm)}")

src_dt = set(pd.read_excel(os.path.join(F6, "module_trait_associate_pep.xlsx"),
                           sheet_name="DT", header=None).iloc[:, 0].astype(str))
src_st = set(pd.read_excel(os.path.join(F6, "module_trait_associate_pep.xlsx"),
                           sheet_name="ST", header=None).iloc[:, 0].astype(str))
print(f"   source table module_trait_associate_pep.xlsx: DT {len(src_dt)} == veen {ok(src_dt == dtm)}"
      f" | ST {len(src_st)} == veen {ok(src_st == stm)}")

# ---------- 4) High-TOM pair-file cross-check ----------
pair_novel = {}
n_pairs = {}
for m in MODS:
    p = pd.read_csv(os.path.join(W, "module_specific_results", f"{m}_high_conf_pairs.txt"), sep="\t")
    n_pairs[m] = len(p)
    pair_novel[m] = set(p["novel"])
u245 = set().union(*pair_novel.values())
print(f"4) High-TOM pair novel union = {len(u245)} (wt_tom_all_id convention ~245/246)"
      f" | per-module novel {({m: len(pair_novel[m]) for m in MODS})}")

# ---------- 5) Candidate chain ----------
cand_dt = dtm & dt_deseq
cand_st = stm & st_deseq
print(f"5) Candidate chain: DT {len(dt_deseq)} intersect 231 -> {len(cand_dt)} (expect 35) {ok(len(cand_dt) == 35)}"
      f" | ST {len(st_deseq)} intersect 181 -> {len(cand_st)} (expect 27) {ok(len(cand_st) == 27)}")

rep = pd.DataFrame([dict(Module=m, N_genes=sizes[m], N_NCP=len(ncp_mod.get(m, set())),
                         High_conf_pairs=n_pairs[m], N_pair_novels=len(pair_novel[m]),
                         In_DT_list=len(dtm & ncp_mod.get(m, set())),
                         In_ST_list=len(stm & ncp_mod.get(m, set()))) for m in MODS]
                   + [dict(Module="grey", N_genes=len(grey), N_NCP=len(grey_ncp),
                           High_conf_pairs=0, N_pair_novels=0, In_DT_list=0, In_ST_list=0)])
rep_path = os.path.join(OUT, "stress_wgcna_reconstruction_report.csv")
rep.to_csv(rep_path, index=False, encoding="utf-8-sig")
print(f"saved: {rep_path}")
