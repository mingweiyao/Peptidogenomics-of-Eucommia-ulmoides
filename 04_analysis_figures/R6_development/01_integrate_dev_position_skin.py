from pathlib import Path
import pandas as pd, numpy as np

B = str(Path(__file__).resolve().parents[2])
MS = f"{B}/output/03_ms_trans_analysis"
A = f"{B}/output/04_article/R6_development"
D = f"{A}/data"

def sig(panel, comp):
    d = pd.read_csv(f"{MS}/01_deseq2/{panel}/{comp}_results.tsv", sep="\t")
    d = d.rename(columns={d.columns[0]: "tid"})
    return set(d.loc[(d["padj"] < 0.05) & (d["log2FoldChange"].abs() >= 1), "tid"])

D_SET = set()
for c in ["growing_leaves_vs_leaf_buds", "young_leaves_vs_leaf_buds",
          "mature_leaves_vs_leaf_buds"]:
    D_SET |= sig("leaf_dev", c)
P_SET = set()
for c in ["close_to_tip_vs_top", "middle_and_lower_vs_top", "lower_vs_top"]:
    P_SET |= sig("leaf_position", c)
S_SET = sig("leaf_position", "skin_vs_top")

checks = {
    "|D|": (len(D_SET), 695), "|P|": (len(P_SET), 458), "|S|": (len(S_SET), 345),
    "D∩P": (len(D_SET & P_SET), 296), "D∩S": (len(D_SET & S_SET), 231),
    "P∩S": (len(P_SET & S_SET), 195), "D∩P∩S": (len(D_SET & P_SET & S_SET), 149),
    "D∩(P∪S)": (len(D_SET & (P_SET | S_SET)), 378),
    "D∩P∖S": (len((D_SET & P_SET) - S_SET), 147),
    "D∩S∖P": (len((D_SET & S_SET) - P_SET), 82),
}
ok = True
for k, (got, tgt) in checks.items():
    flag = "OK" if got == tgt else "MISMATCH"
    if got != tgt:
        ok = False
    print(f"{k:12s} got={got:4d} target={tgt:4d} {flag}")
assert ok, "Three-set reconciliation failed"


rows = [
    ("D only", len(D_SET - P_SET - S_SET)), ("P only", len(P_SET - D_SET - S_SET)),
    ("S only", len(S_SET - D_SET - P_SET)),
    ("D∩P∖S", len((D_SET & P_SET) - S_SET)),
    ("D∩S∖P", len((D_SET & S_SET) - P_SET)),
    ("P∩S∖D", len((P_SET & S_SET) - D_SET)),
    ("D∩P∩S", len(D_SET & P_SET & S_SET)),
]
pd.DataFrame(rows, columns=["exclusive_region", "n"]).to_csv(
    f"{D}/R6_upset_exclusive_counts.tsv", sep="\t", index=False)
print("Exclusive regions:", rows)

targets = pd.read_csv(f"{MS}/00_input/ncp_unique_transcripts.tsv", sep="\t",
                      header=None, names=["transcript_id"])
m = targets.copy()
m["in_D"] = m["transcript_id"].isin(D_SET)
m["in_P"] = m["transcript_id"].isin(P_SET)
m["in_S"] = m["transcript_id"].isin(S_SET)
def cls(r):
    parts = [x for x, f in [("D", r["in_D"]), ("P", r["in_P"]), ("S", r["in_S"])] if f]
    return "∩".join(parts) if parts else "none"
m["evidence_class"] = m.apply(cls, axis=1)
m.to_csv(f"{D}/R6_sets_D_P_S_membership.tsv", sep="\t", index=False)
print("Wrote membership:", len(m), "| class distribution:", m["evidence_class"].value_counts().to_dict())
print("DONE 01")
