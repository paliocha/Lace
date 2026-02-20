#!/usr/bin/env python3
"""Downsample BMAX to ~10% of clusters for a quick Lace example."""

import random
import sys
from pathlib import Path

SEED = 42
FRACTION = 0.01

CLUSTERS_IN = Path(
    "/mnt/project/FjellheimLab/martpali/AnnualPerennial/"
    "nf-denovoslim/BMAX/clustering/corset-clusters.txt"
)
FASTA_IN = Path(
    "/mnt/project/FjellheimLab/martpali/AnnualPerennial/"
    "assemblies/BMAX-Trinity1/BMAX-Trinity.fasta"
)

OUT_DIR = Path(__file__).resolve().parent
CLUSTERS_OUT = OUT_DIR / "bmax_1pct_clusters.txt"
FASTA_OUT = OUT_DIR / "bmax_1pct.fasta"

# --- 1. Sample 10% of unique cluster names ---
all_lines = CLUSTERS_IN.read_text().strip().split("\n")
cluster_to_lines: dict[str, list[str]] = {}
for line in all_lines:
    tid, cid = line.split("\t")
    cluster_to_lines.setdefault(cid, []).append(line)

all_clusters = sorted(cluster_to_lines.keys())
random.seed(SEED)
n_keep = max(1, int(len(all_clusters) * FRACTION))
keep_clusters = set(random.sample(all_clusters, n_keep))

# Write subset clusters file
kept_tids: set[str] = set()
with open(CLUSTERS_OUT, "w", encoding="utf-8") as fh:
    for cid in sorted(keep_clusters):
        for line in cluster_to_lines[cid]:
            fh.write(line + "\n")
            kept_tids.add(line.split("\t")[0])

print(f"Clusters kept: {len(keep_clusters):,} / {len(all_clusters):,}")
print(f"Transcripts kept: {len(kept_tids):,} / {len(all_lines):,}")
print(f"Wrote {CLUSTERS_OUT}")

# --- 2. Extract matching FASTA records ---
written = 0
with open(FASTA_IN, encoding="utf-8") as fin, \
     open(FASTA_OUT, "w", encoding="utf-8") as fout:
    keep = False
    for line in fin:
        if line.startswith(">"):
            tid = line[1:].split()[0]
            keep = tid in kept_tids
            if keep:
                written += 1
        if keep:
            fout.write(line)

print(f"FASTA records written: {written:,}")
print(f"Wrote {FASTA_OUT}")
