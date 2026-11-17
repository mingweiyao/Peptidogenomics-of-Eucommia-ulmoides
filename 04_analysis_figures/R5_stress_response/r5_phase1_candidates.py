from pathlib import Path
import pandas as pd, numpy as np

B = str(Path(__file__).resolve().parents[2])
MS = f"{B}/output/03_ms_trans_analysis"
Q = f"{B}/output/02filter_trans/07new_transcript_quant"
ART = f"{B}/output/04_article"
D = f"{ART}/R5_stress_response/data"

COMPS = ["DT_48_vs_CK", "DT_144_vs_CK", "ST_48_vs_CK", "ST_144_vs_CK"]
COND_OF_COMP = {"DT_48_vs_CK": "DT_48", "DT_144_vs_CK": "DT_144",
                "ST_48_vs_CK": "ST_48", "ST_144_vs_CK": "ST_144",
                "heat_vs_heat_CK": "heat"}
PROJECT_COMPS = {"DT_ST": COMPS, "heat": ["heat_vs_heat_CK"]}


def load_results(panel, comp):
    df = pd.read_csv(f"{MS}/01_deseq2/{panel}/{comp}_results.tsv", sep="\t")
    return df.rename(columns={df.columns[0]: "transcript_id"}).set_index("transcript_id")


def deg_set(df):
    return set(df.index[(df["padj"] < 0.05) & (df["log2FoldChange"].abs() >= 1)])



targets = pd.read_csv(f"{MS}/00_input/ncp_unique_transcripts.tsv", sep="\t",
                      header=None, names=["transcript_id"])
pairs = pd.read_csv(f"{MS}/00_input/ncp_to_transcript.tsv", sep="\t")
pairs = pairs[pairs["transcript_id"].isin(set(targets["transcript_id"]))]
ncp_ev = pairs.groupby("transcript_id").agg(
    n_ncp_peptides=("ncp_peptide", "nunique"),
    psms_sum=("psms", "sum"), psms_max=("psms", "max"),
    ncp_class_codes=("class_code", lambda s: ";".join(sorted(set(s)))))


res = {c: load_results("ST_DT", c) for c in COMPS}
res["heat_vs_heat_CK"] = load_results("heat", "heat_vs_heat_CK")
de = pd.DataFrame(index=targets["transcript_id"])
for c, df in res.items():
    tag = COND_OF_COMP[c]
    de[f"{tag}_LFC"] = df["log2FoldChange"].reindex(de.index)
    de[f"{tag}_padj"] = df["padj"].reindex(de.index)
    de[f"{tag}_DEG"] = de.index.isin(deg_set(df))


sd = pd.read_csv(f"{MS}/00_input/sample_design.tsv", sep="\t")
st_samp = sd.loc[sd["panel"] == "ST_DT", "sample"].tolist()
ht_samp = sd.loc[sd["panel"] == "heat", "sample"].tolist()
use = ["transcript_id", "class_code"] + st_samp + ht_samp
tpm = pd.read_csv(f"{Q}/TPM_matrix.tsv", sep="\t", usecols=use).set_index("transcript_id")
cond_of = sd.set_index("sample")["condition"].to_dict()

def cond_groups(tid, samples):
    """Cond groups."""
    row = tpm.loc[tid, samples]
    g = {}
    for s, v in zip(samples, row.values):
        g.setdefault(cond_of[s], []).append(float(v))
    return {k: pd.Series(v) for k, v in g.items()}

def display_pass(tid, samples):
    g = cond_groups(tid, samples)
    ok = {c: int((v >= 1).sum()) >= 2 for c, v in g.items()}
    return any(ok.values()), ok

def consistency(tid, samples, comp, direction, ref_cond="CK"):
    """Consistency."""
    top = COND_OF_COMP[comp]
    g = cond_groups(tid, samples)
    ck = g.get(ref_cond)
    tv = g.get(top)
    if ck is None or tv is None or len(ck) == 0 or len(tv) == 0:
        return np.nan
    thr = ck.max() if direction > 0 else ck.min()
    hit = (tv > thr) if direction > 0 else (tv < thr)
    return float(hit.sum()) / len(tv)


iso = pd.read_csv(f"{ART}/tables/TableS6_locus_isoforms.tsv", sep="\t").set_index("gene_id")

# ---- ST_DT WGCNA(v4) ----
gm = pd.read_csv(f"{MS}/02_wgcna/ST_DT/gene_module_assignment.tsv", sep="\t").set_index("transcript_id")
me_c = pd.read_csv(f"{MS}/02_wgcna/ST_DT/ME_trait_correlation.tsv", sep="\t", index_col=0)
me_padj = pd.read_csv(f"{MS}/02_wgcna/ST_DT/ME_trait_padj.tsv", sep="\t", index_col=0)
mod_link = {}   
for m in me_padj.index:
    hit = [t for t in me_padj.columns if me_padj.loc[m, t] < 0.05 and abs(me_c.loc[m, t]) >= 0.6]
    if hit:
        best = max(hit, key=lambda t: abs(me_c.loc[m, t]))
        mod_link[m] = f"{best}({me_c.loc[m, best]:+.2f},FDR={me_padj.loc[m, best]:.3f})"

cc = pd.read_csv(f"{MS}/01_deseq2/cross_candidates.tsv", sep="\t")
cc_st = cc[cc["evidence_panel"] == "ST_DT"]
cc_agg = cc_st.groupby("transcript_id").agg(
    cross_tier=("tier", "min"),
    cross_traits=("trait", lambda s: ";".join(sorted(set(s)))))


t = targets.set_index("transcript_id").copy()
t = t.join(ncp_ev).join(de)
t["locus"] = pairs.groupby("transcript_id")["locus"].first().reindex(t.index)
t["gffcompare_class"] = tpm["class_code"].reindex(t.index)
t["locus_total_isoforms"] = t["locus"].map(iso["total_assembled_isoforms"])
t["locus_single_isoform"] = (t["locus_total_isoforms"] == 1)

t["in_DT_union"] = t["DT_48_DEG"].fillna(False) | t["DT_144_DEG"].fillna(False)
t["in_ST_union"] = t["ST_48_DEG"].fillna(False) | t["ST_144_DEG"].fillna(False)
t["in_union393"] = t["in_DT_union"] | t["in_ST_union"]
t["in_common150"] = t["in_DT_union"] & t["in_ST_union"]
t["in_heat206"] = t["heat_DEG"].fillna(False)
t["in_cross88"] = t["in_union393"] & t["in_heat206"]


for samples, tag in [(st_samp, "st"), (ht_samp, "ht")]:
    means, pdisp = {}, {}
    for tid in t.index:
        g = cond_groups(tid, samples)
        means[tid] = {c: v.mean() for c, v in g.items()}
        ok, _ = display_pass(tid, samples)
        pdisp[tid] = ok
    t[f"pass_display_{tag}"] = pd.Series(pdisp)
    for cond, vals in pd.DataFrame(means).T.items():
        t[f"tpm_mean_{tag}_{cond}"] = vals.round(2)


def driver_of(row, comps):
    sig = [(c, abs(row[f"{COND_OF_COMP[c]}_LFC"]), row[f"{COND_OF_COMP[c]}_LFC"])
           for c in comps if row[f"{COND_OF_COMP[c]}_DEG"]]
    if not sig:
        return None
    return max(sig, key=lambda x: x[1])  # (comp, |lfc|, lfc)

drv_st, dir_st, cons_st, drv_ht, dir_ht, cons_ht = [], [], [], [], [], []
for tid, row in t.iterrows():
    d = driver_of(row, COMPS)
    if d is None:
        drv_st.append(""); dir_st.append(np.nan); cons_st.append(np.nan)
    else:
        drv_st.append(d[0]); dir_st.append(np.sign(d[2]))
        cons_st.append(consistency(tid, st_samp, d[0], np.sign(d[2]), "CK"))
    d = driver_of(row, ["heat_vs_heat_CK"])
    if d is None:
        drv_ht.append(""); dir_ht.append(np.nan); cons_ht.append(np.nan)
    else:
        drv_ht.append(d[0]); dir_ht.append(np.sign(d[2]))
        cons_ht.append(consistency(tid, ht_samp, d[0], np.sign(d[2]), "heat_CK"))
t["driver_comp_dtst"] = drv_st
t["driver_direction_dtst"] = dir_st
t["consistency_dtst"] = np.round(cons_st, 2)
t["driver_comp_heat"] = drv_ht
t["driver_direction_heat"] = dir_ht
t["consistency_heat"] = np.round(cons_ht, 2)


t["ST_DT_module"] = gm["module"].reindex(t.index).fillna("")
t["ST_DT_kME"] = gm["kME_own_module"].reindex(t.index)
t["module_FDR_link"] = t["ST_DT_module"].map(mod_link).fillna("")
t = t.join(cc_agg)


def effect_s(row, comps):
    vals = [abs(row[f"{COND_OF_COMP[c]}_LFC"]) for c in comps
            if row[f"{COND_OF_COMP[c]}_DEG"]]
    return min(max(vals) / 5.0, 1.0) if vals else 0.0

def ncp_s(row):
    return min(1.0, 0.4 + 0.2 * min(row["n_ncp_peptides"] - 1, 2) / 2
               + 0.4 * min(row["psms_sum"], 5) / 5)

def module_s(row):
    if not row["ST_DT_module"] or row["ST_DT_module"] not in mod_link:
        return 0.0
    k = row["ST_DT_kME"]
    return 0.15 if (pd.notna(k) and k >= 0.7) else 0.05

def score(row, comps, tag):
    eff = effect_s(row, comps)
    cons = row["consistency_dtst" if tag == "st" else "consistency_heat"]
    con = 0.0 if pd.isna(cons) else cons
    condcols = [c for c in t.columns if c.startswith(f"tpm_mean_{tag}_")]
    abund = min(np.log10(max(row[condcols].max(), 0) + 1) / 3.0, 1.0)
    iso = 0.1 if row["locus_single_isoform"] else 0.0
    mod = module_s(row) if tag == "st" else 0.0   
    cross = 0.1 if row["in_cross88"] else 0.0
    return round(0.30 * eff + 0.20 * con + 0.15 * abund + 0.20 * ncp_s(row) + iso + mod + cross, 3)

t["score_dtst"] = [score(r, COMPS, "st") if r["in_union393"] else np.nan for _, r in t.iterrows()]
t["score_heat"] = [score(r, ["heat_vs_heat_CK"], "ht") if r["in_heat206"] else np.nan for _, r in t.iterrows()]






DT_TAGS = ["DT_48", "DT_144"]
ST_TAGS = ["ST_48", "ST_144"]
ARROW = {1.0: "▲", -1.0: "▼"}

def sig_pat(row):
    return {tg: np.sign(row[f"{tg}_LFC"]) for tg in DT_TAGS + ST_TAGS
            if row[f"{tg}_DEG"]}

def tc_label(pat):
    if not pat:
        return ""
    nm = {"DT_48": "DT48", "DT_144": "DT144", "ST_48": "ST48", "ST_144": "ST144"}
    if len(pat) == 1:
        k, v = next(iter(pat.items()))
        return f"{nm[k]} {ARROW[v]}"
    dts = [v for k, v in pat.items() if k.startswith("DT")]
    sts = [v for k, v in pat.items() if k.startswith("ST")]
    if len(set(pat.values())) == 1:
        a = ARROW[next(iter(pat.values()))]
        if dts and sts:
            return f"DT+ST sustained {a}"
        return f"{'DT' if dts else 'ST'} sustained {a}"
    return "time-dep. " + "·".join(f"{nm[k]}{ARROW[v]}" for k, v in pat.items())

pats = {tid: sig_pat(r) for tid, r in t.iterrows()}
t["sig_pattern"] = pd.Series({tid: "·".join(
    f"{k.replace('_','')}{ARROW[v]}" for k, v in pats[tid].items()) for tid in t.index})
t["tc_label"] = pd.Series({tid: tc_label(pats[tid]) for tid in t.index})
t["n_sig_comps"] = pd.Series({tid: len(pats[tid]) for tid in t.index})
t["tc_class"] = pd.Series({tid: (
    "mixed" if len(set(pats[tid].values())) > 1 else
    ("multi" if len(pats[tid]) > 1 else "single")) for tid in t.index})
t["psm_tier"] = np.where(t["psms_sum"] >= 3, "strong", "low")
t["heat_tc_label"] = [
    ("heat " + ARROW[np.sign(r["heat_LFC"])]) if r["in_heat206"] else ""
    for _, r in t.iterrows()]

t.loc[t["in_heat206"] & t["in_union393"], "heat_tc_label"] = [
    f"heat {ARROW[np.sign(r['heat_LFC'])]} | {r['sig_pattern']}"
    for _, r in t[t["in_heat206"] & t["in_union393"]].iterrows()]
t["selected_dtst"] = ""
t["selected_heat"] = ""
sel_st = []

def pick_one(pool, col="score_dtst", prefer_strong=True):
    pool = pool[~pool.index.isin(sel_st)]
    if not len(pool):
        return None

    if prefer_strong and (pool["psm_tier"] == "strong").any():
        pool = pool[pool["psm_tier"] == "strong"]
    return pool.sort_values([col, "DT_144_padj" if col == "score_dtst" else "heat_padj"],
                            ascending=[False, True]).index[0]

pool_st = t[t["in_union393"] & t["pass_display_st"]]

def has_stress_sig(row, kind, dsign):
    tags = DT_TAGS if kind == "DT" else ST_TAGS
    return any(row[f"{tg}_DEG"] and np.sign(row[f"{tg}_LFC"]) == dsign for tg in tags)

def no_opposite_sig(row, dsign):
    """No opposite sig."""
    return not any(row[f"{tg}_DEG"] and np.sign(row[f"{tg}_LFC"]) == -dsign
                   for tg in DT_TAGS + ST_TAGS)

def sustained_pool(kind, dsign):
    m = pool_st.apply(lambda r: has_stress_sig(r, kind, dsign) and no_opposite_sig(r, dsign),
                      axis=1)
    return pool_st[m.astype(bool)]

STRATA = [("drought_up", "DT", 1), ("drought_down", "DT", -1),
          ("salt_up", "ST", 1), ("salt_down", "ST", -1)]
sel_st = []
for name, kind, dsign in STRATA:
    tid = pick_one(sustained_pool(kind, dsign))
    if tid:
        sel_st.append(tid); t.loc[tid, "selected_dtst"] = name

mixed_pool = pool_st[(pool_st["tc_class"] == "mixed").astype(bool)
                     & ~pool_st.index.isin(sel_st)]
tid = pick_one(mixed_pool)
if tid:
    sel_st.append(tid); t.loc[tid, "selected_dtst"] = "time_dependent"
rest = pool_st[~pool_st.index.isin(sel_st)]
prefer = rest[rest["in_cross88"] | rest["in_common150"]]
tid = pick_one(prefer) if len(prefer) else pick_one(rest)
if tid:
    sel_st.append(tid); t.loc[tid, "selected_dtst"] = "cross_or_common_top"
while len(sel_st) < 6:
    tid = pick_one(rest)
    if not tid:
        break
    sel_st.append(tid); t.loc[tid, "selected_dtst"] = "top_score_fill"
assert len(set(sel_st)) == len(sel_st), "Selection contains duplicate IDs"

pool_ht = t[t["in_heat206"] & t["pass_display_ht"]]
sel_ht = []
for name, dsign in [("heat_up", 1), ("heat_down", -1)]:
    pool = pool_ht[pool_ht["driver_direction_heat"] == dsign]
    if len(pool):
        tid = pool.sort_values(["score_heat", "heat_padj"], ascending=[False, True]).index[0]
        sel_ht.append(tid); t.loc[tid, "selected_heat"] = name

rest_ht = pool_ht[~pool_ht.index.isin(sel_ht)]
while len(sel_ht) < 4 and len(rest_ht):
    n_up = sum(1 for x in sel_ht if t.loc[x, "driver_direction_heat"] > 0)
    n_dn = sum(1 for x in sel_ht if t.loc[x, "driver_direction_heat"] < 0)
    want = 1 if n_up <= n_dn else -1
    sub = rest_ht[rest_ht["driver_direction_heat"] == want]
    pool = sub if len(sub) else rest_ht
    tid = pool.sort_values(["score_heat", "heat_padj"], ascending=[False, True]).index[0]
    sel_ht.append(tid)
    t.loc[tid, "selected_heat"] = ("heat_up" if t.loc[tid, "driver_direction_heat"] > 0 else "heat_down")
    rest_ht = rest_ht.drop(index=tid)
assert len(set(sel_ht)) == len(sel_ht), "Heat selection contains duplicate IDs"


def plot_table(tids, samples, project):
    rows = []
    rep_of = {}
    for s in samples:
        rep_of[s] = rep_of.get(cond_of[s], 0) + 1
    for tid in tids:
        for s in samples:
            rows.append({"transcript_id": tid, "project": project,
                         "condition": cond_of[s], "sample": s,
                         "replicate": rep_of[s], "TPM": round(float(tpm.loc[tid, s]), 3)})
    return pd.DataFrame(rows)

plot_table(sel_st, st_samp, "DT_ST").to_csv(f"{D}/R5_DT_ST_expression_plot.tsv", sep="\t", index=False)
plot_table(sel_ht, ht_samp, "heat").to_csv(f"{D}/R5_heat_expression_plot.tsv", sep="\t", index=False)


t.reset_index().to_csv(f"{D}/R5_candidate_ranking.tsv", sep="\t", index=False)
sel_cols = ["selected_dtst", "selected_heat", "score_dtst", "score_heat",
            "tc_label", "sig_pattern", "tc_class", "psm_tier",
            "driver_comp_dtst", "driver_direction_dtst", "consistency_dtst",
            "driver_comp_heat", "driver_direction_heat", "consistency_heat",
            "gffcompare_class", "locus_total_isoforms", "n_ncp_peptides", "psms_sum",
            "ST_DT_module", "ST_DT_kME", "module_FDR_link", "cross_tier", "cross_traits",
            "in_cross88", "tpm_mean_st_CK", "tpm_mean_st_DT_48", "tpm_mean_st_DT_144",
            "tpm_mean_st_ST_48", "tpm_mean_st_ST_144",
            "tpm_mean_ht_heat_CK", "tpm_mean_ht_heat"]
t.loc[sel_st + sel_ht, sel_cols].to_csv(f"{D}/R5_selected_candidates.tsv", sep="\t")

print("== Pool after display filter ==")
print("DT/ST pool:", len(pool_st), "| heat pool:", len(pool_ht))
print("  DT∪ & dir>0:", int((pool_st["in_DT_union"] & (pool_st["driver_direction_dtst"] > 0)).sum()),
      "| DT∪ & dir<0:", int((pool_st["in_DT_union"] & (pool_st["driver_direction_dtst"] < 0)).sum()),
      "| ST∪ & dir>0:", int((pool_st["in_ST_union"] & (pool_st["driver_direction_dtst"] > 0)).sum()),
      "| ST∪ & dir<0:", int((pool_st["in_ST_union"] & (pool_st["driver_direction_dtst"] < 0)).sum()))
print("  heat_up:", int((pool_ht["driver_direction_heat"] > 0).sum()),
      "| heat_down:", int((pool_ht["driver_direction_heat"] < 0).sum()))
print("\n== DT/ST selected ==")
import pandas as _pd; _pd.set_option("display.max_columns",50); print(t.loc[sel_st, [c for c in sel_cols if c in t.columns]].to_string())
print("\n== Heat selected ==")
print(t.loc[sel_ht, [c for c in sel_cols if c in t.columns]].to_string())
print("\nDONE phase1")
