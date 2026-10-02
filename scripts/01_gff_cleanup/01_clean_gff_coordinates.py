# -*- coding: utf-8 -*-
"""GFF coordinate cleanup in three legacy steps (faithful replay): 336,380 -> 298,877 -> 295,224 -> 223,870.

Steps (replaying ac_gff_make.py's remove_fully_contained / process_junctions / filter_overlapping_regions):
  1) deduplicate the five non-CDS classes by (chrom,start,end,strand) = 336,380
  2) remove_fully_contained = 298,877 (legacy behaviour: drops outer intervals that contain others)
  3) process_junctions (drop junction intervals overlapping CDS) = 295,224
  4) filter_overlapping_regions = 223,870 (legacy behaviour: groups by chrom only, not strand - kept for reproduction)

Inputs: output/_source_cache/par_Annotated_Peptides.pkl, par_Genomic_Features.pkl (from 01_prepare_source_tables.py)
Outputs: output/_base_data_chain/gff_223870.pkl / gff_223870.csv
Requires: pip install intervaltree
Usage: py 01_clean_gff_coordinates.py --out ./output   (about 3-10 min)
"""
import pandas as pd
import numpy as np
from collections import defaultdict
from intervaltree import IntervalTree
import time
import os as _os
import argparse
import os

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))

DATA = OUTROOT + "/_source_cache"
OUT = OUTROOT + "/_base_data_chain"

ap = pd.read_pickle(f"{DATA}/par_Annotated_Peptides.pkl")

TYPE = {"three_prime_UTR": "3'UTR", "five_prime_UTR": "5'UTR"}

# ---- Step 1: GFF sheet 336,380 ----
noncds = ap["type"].isin(["intergenic", "intron", "junction", "three_prime_UTR", "five_prime_UTR"])
gff = ap.loc[noncds, ["accessions", "start", "end", "strand", "chrom", "type"]].copy()
gff["type"] = gff["type"].replace(TYPE)
gff = gff.drop_duplicates(subset=["chrom", "start", "end", "strand"], keep="first").reset_index(drop=True)
gff["start"] = gff["start"].astype(int)
gff["end"] = gff["end"].astype(int)
print(f"Step 1 GFF: {len(gff)} (target 336,380)", flush=True)


def remove_fully_contained(df):
    df = df.copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    out = []
    for (chrom, strand), group in df.groupby(["chrom", "strand"], sort=False):
        tree = IntervalTree()
        for idx, row in group.iterrows():
            s, e = row["start"], row["end"]
            if s > e:
                s, e = e, s
            if s < e:
                tree[s:e] = idx
        to_remove = set()
        for interval in tree:
            containing = tree.envelop(interval.begin, interval.end)
            for other in containing:
                if other.data != interval.data:
                    to_remove.add(interval.data)
        out.append(group.drop(index=to_remove))
    return pd.concat(out, ignore_index=True)


def process_junctions(df, gf):
    cds_trees = defaultdict(IntervalTree)
    for _, row in gf.iterrows():
        key = (row["chrom"], row["strand"])
        if int(row["start"]) < int(row["end"]) and row["type"] == "CDS":
            cds_trees[key].addi(int(row["start"]), int(row["end"]))
    to_drop = []
    mask = df["type"].eq("junction")
    for idx, row in df[mask].iterrows():
        s, e = int(row["start"]), int(row["end"])
        if s >= e:
            continue
        key = (row["chrom"], row["strand"])
        if key in cds_trees and cds_trees[key].overlaps(s, e):
            to_drop.append(idx)
    return df.drop(index=to_drop)


def filter_overlapping_regions(df):
    filtered = []
    for chrom, group in df.groupby("chrom", sort=False):
        group = group.sort_values("start").reset_index(drop=True)
        n = len(group)
        starts = group["start"].to_numpy()
        ends = group["end"].to_numpy()
        tree = IntervalTree()
        for idx in range(n):
            tree.addi(int(starts[idx]), int(ends[idx]), idx)
        to_remove = set()
        starts_sorted = sorted(tree)
        for idx in range(n):
            if idx in to_remove:
                continue
            s, e = int(starts[idx]), int(ends[idx])
            overlaps = tree.overlap(s, e)
            overlapping_indices = {iv.data for iv in overlaps}
            if len(overlapping_indices) == 1:
                continue
            overlapping_indices = sorted(overlapping_indices)
            if len(overlapping_indices) == 2:
                other_idx = overlapping_indices[1] if overlapping_indices[0] == idx else overlapping_indices[0]
                if (ends[other_idx] - starts[other_idx]) > (e - s):
                    to_remove.add(idx)
                else:
                    to_remove.add(other_idx)
            else:
                best_remove = None
                max_non_overlap = 0
                for candidate in overlapping_indices:
                    temp = to_remove | {candidate}
                    remaining = [i for i in overlapping_indices if i not in temp]
                    has_overlap = False
                    for a in range(len(remaining)):
                        for b in range(a + 1, len(remaining)):
                            if ends[remaining[a]] > starts[remaining[b]]:
                                has_overlap = True
                                break
                        if has_overlap:
                            break
                    if not has_overlap:
                        if len(remaining) > max_non_overlap:
                            max_non_overlap = len(remaining)
                            best_remove = candidate
                if best_remove is not None:
                    to_remove.add(best_remove)
                else:
                    longest_idx = overlapping_indices[0]
                    max_len = ends[longest_idx] - starts[longest_idx]
                    for i in overlapping_indices[1:]:
                        L = ends[i] - starts[i]
                        if L > max_len:
                            longest_idx = i
                            max_len = L
                    for i in overlapping_indices:
                        if i != longest_idx:
                            to_remove.add(i)
        filtered.append(group[~group.index.isin(to_remove)])
    return pd.concat(filtered, ignore_index=True)


t0 = time.time()
print("Step 2 remove_fully_contained ...", flush=True)
df1 = remove_fully_contained(gff)
print(f"  -> {len(df1)} (target 298,877) [{time.time()-t0:.0f}s]", flush=True)

print("Step 3 process_junctions ...", flush=True)
gf = pd.read_pickle(f"{DATA}/par_Genomic_Features.pkl")
df2 = process_junctions(df1, gf)
print(f"  -> {len(df2)} (target 295,224)", flush=True)

print("Step 4 filter_overlapping_regions (legacy: groups by chrom only) ...", flush=True)
df3 = filter_overlapping_regions(df2)
print(f"  -> {len(df3)} (target 223,870)", flush=True)

df3.to_pickle(f"{OUT}/gff_223870.pkl")
df3.to_csv(f"{OUT}/gff_223870.csv", index=False)
print("Done: gff_223870", flush=True)
