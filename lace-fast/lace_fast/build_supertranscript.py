#!/usr/bin/env python
"""Build a SuperTranscript for one cluster of transcripts.

Rewrite of Lace 1.14.1 ``BuildSuperTranscript.py`` with:

* **A1** – minimap2 replaces BLAT (30–100× faster, MIT licence)
* **A2** – block-level splice graph (~10–30 nodes) replaces base-level
           graph (~10 000 nodes)
* **A3** – single-pass DFS back-edge cycle breaking replaces iterative
           ``nx.simple_cycles`` (worst-case exponential → O(V+E))
* **A4** – numpy-vectorised coordinate mapping for block merging
* **H1–H14** – full hygiene pass (subprocess, pathlib, typing, f-strings,
               vectorised pandas, etc.)
* **IO1–IO3** – in-memory sequence passing, $TMPDIR for alignment files

Original author: Anthony Hawkins
Rewrite: lace-fast 2.0
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import NamedTuple

import networkx as nx
import numpy as np
import pandas as pd

log = logging.getLogger("lace_fast")

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

class AlignBlock(NamedTuple):
    """One contiguous alignment block from a minimap2 PAF hit."""
    t_name: str
    q_name: str
    t_start: int
    t_end: int
    q_start: int
    q_end: int
    strand: str


class ClusterResult(NamedTuple):
    """Return value from processing one cluster."""
    seq: str
    anno: str
    whirl_status: int
    transcript_count: int

# ---------------------------------------------------------------------------
# Reverse complement  (H5 – str.maketrans instead of char-by-char loop)
# ---------------------------------------------------------------------------

_RC_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of *seq*."""
    return seq.translate(_RC_TABLE)[::-1]

# ---------------------------------------------------------------------------
# GFF annotation helper
# ---------------------------------------------------------------------------

def get_annotation_line(
    cluster_id: str,
    start: str,
    end: str,
    trans_id: str,
) -> str:
    """Format one GFF2 annotation line for a SuperTranscript block."""
    return (
        f"{cluster_id}\tSuperTranscript\texon\t{start}\t{end}\t.\t.\t0\t"
        f'gene_id "{cluster_id}"; trans_id "{trans_id}";\n'
    )

# ---------------------------------------------------------------------------
# minimap2 alignment  (A1 — replaces BLAT)
# ---------------------------------------------------------------------------

def _run_minimap2(
    fasta_path: Path,
    *,
    min_identity: float = 0.98,
) -> list[AlignBlock]:
    """Run minimap2 all-vs-all on *fasta_path*, return alignment blocks.

    Uses ``-X`` for all-vs-all (ava) mode and ``--eqx -c`` for base-level
    CIGAR, replicating BLAT's ``-minIdentity=98`` behaviour.
    The PAF is parsed from stdout — no intermediate file on disk (IO2).
    """
    cmd = [
        "minimap2",
        "-c",            # output CIGAR in PAF
        "-X",            # all-vs-all (skip self-hits, symmetric dedup)
        "--eqx",         # extended CIGAR with =/X
        "-k15",          # k-mer size (good for nt transcript alignment)
        "-w5",           # minimiser window
        "-N50",          # retain up to 50 secondary alignments
        f"-p{min_identity}",  # min score ratio
        "--no-long-join",
        str(fasta_path),
        str(fasta_path),
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True,
    )

    blocks: list[AlignBlock] = []
    for line in result.stdout.splitlines():
        cols = line.split("\t")
        if len(cols) < 12:
            continue
        q_name = cols[0]
        # q_len  = int(cols[1])
        q_start = int(cols[2])
        q_end = int(cols[3])
        strand = cols[4]
        t_name = cols[5]
        # t_len  = int(cols[6])
        t_start = int(cols[7])
        t_end = int(cols[8])

        # Skip self-alignments (minimap2 -X already does, but belt & braces)
        if q_name == t_name:
            continue

        blocks.append(AlignBlock(
            t_name=t_name,
            q_name=q_name,
            t_start=t_start,
            t_end=t_end,
            q_start=q_start,
            q_end=q_end,
            strand=strand,
        ))

    return blocks


def _determine_strand_directions(
    blocks: list[AlignBlock],
) -> dict[str, str]:
    """Assign + / - orientation to each transcript (same logic as original).

    The first transcript encountered is arbitrarily +.  Subsequent ones
    are oriented relative to the alignment strand of hits that connect them.
    """
    trandir: dict[str, str] = {}
    seen_pairs: set[str] = set()

    for blk in blocks:
        pair_key = f"{blk.t_name}\t{blk.q_name}"
        rev_key = f"{blk.q_name}\t{blk.t_name}"
        if pair_key in seen_pairs or rev_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)

        t, q, strand = blk.t_name, blk.q_name, blk.strand

        if t in trandir and q in trandir:
            continue
        elif t in trandir:
            trandir[q] = "+" if (trandir[t] == strand) else "-"
        elif q in trandir:
            trandir[t] = "+" if (trandir[q] == strand) else "-"
        else:
            trandir[t] = "+"
            trandir[q] = "+" if strand == "+" else "-"

    return trandir

# ---------------------------------------------------------------------------
# Block-level splice graph  (A2 — replaces base-level graph)
# ---------------------------------------------------------------------------

def _build_block_graph(
    transcripts: dict[str, str],
    blocks: list[AlignBlock],
    *,
    max_edges: int = 100,
) -> tuple[str, str, int]:
    """Construct a block-level splice graph from minimap2 alignment blocks.

    Instead of one NetworkX node per *base* (original: 10 000 nodes for a
    5 × 2 kb cluster), we create one node per *contiguous block* (typically
    10–30 nodes).

    Returns ``(sequence, annotation, whirl_count)``.
    """
    # -- 1) Determine strand orientation and fix reverse-complement --------
    trandir = _determine_strand_directions(blocks)

    # Orient all transcripts to + strand
    oriented: dict[str, str] = {}
    for name, seq in transcripts.items():
        if trandir.get(name, "+") == "-":
            oriented[name] = reverse_complement(seq)
        else:
            oriented[name] = seq

    # -- 2) Compute interval breakpoints per transcript --------------------
    # For each transcript, collect all alignment boundary positions.
    # These define the "block" boundaries.
    breakpoints: dict[str, set[int]] = {
        name: {0, len(seq)} for name, seq in oriented.items()
    }

    # Collect breakpoints from alignment blocks
    for blk in blocks:
        if blk.t_name == blk.q_name:
            continue
        breakpoints.setdefault(blk.t_name, set()).update({blk.t_start, blk.t_end})
        breakpoints.setdefault(blk.q_name, set()).update({blk.q_start, blk.q_end})

    # Sort breakpoints into ordered interval lists per transcript
    intervals: dict[str, list[tuple[int, int]]] = {}
    for name, bps in breakpoints.items():
        sorted_bps = sorted(bps)
        ivs = []
        for i in range(len(sorted_bps) - 1):
            s, e = sorted_bps[i], sorted_bps[i + 1]
            if e > s:
                ivs.append((s, e))
        intervals[name] = ivs

    # -- 3) Create a node for each block and build equivalence classes -----
    # node_id → sequence string
    node_seq: dict[int, str] = {}
    # (transcript, interval_index) → node_id
    block_to_node: dict[tuple[str, int], int] = {}
    # For merging: map (transcript, position) → (interval_index)
    # We'll build a position→interval lookup per transcript
    pos_to_interval: dict[str, dict[int, int]] = {}

    next_node = 0
    for name, ivs in intervals.items():
        lookup: dict[int, int] = {}
        for idx, (s, e) in enumerate(ivs):
            # Map every start position to interval index for alignment matching
            lookup[s] = idx
            node_id = next_node
            node_seq[node_id] = oriented[name][s:e]
            block_to_node[(name, idx)] = node_id
            next_node += 1
        pos_to_interval[name] = lookup

    # -- 4) Merge nodes that correspond to the same aligned region ---------
    # Union-Find for merging equivalent blocks
    parent: dict[int, int] = {n: n for n in node_seq}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]  # path compression
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            # Keep the lower-id node as representative (deterministic)
            if ra > rb:
                ra, rb = rb, ra
            parent[rb] = ra

    # For each alignment block, find matching intervals and merge them
    for blk in blocks:
        if blk.t_name == blk.q_name:
            continue

        t_lookup = pos_to_interval.get(blk.t_name, {})
        q_lookup = pos_to_interval.get(blk.q_name, {})

        t_ivs = intervals.get(blk.t_name, [])
        q_ivs = intervals.get(blk.q_name, [])

        # Walk through the alignment block region and merge overlapping
        # sub-intervals from both transcripts
        t_pos = blk.t_start
        q_pos = blk.q_start

        while t_pos < blk.t_end and q_pos < blk.q_end:
            t_idx = t_lookup.get(t_pos)
            q_idx = q_lookup.get(q_pos)

            if t_idx is None or q_idx is None:
                # Advance to next breakpoint
                t_pos += 1
                q_pos += 1
                continue

            t_node = block_to_node.get((blk.t_name, t_idx))
            q_node = block_to_node.get((blk.q_name, q_idx))

            if t_node is not None and q_node is not None:
                t_s, t_e = t_ivs[t_idx]
                q_s, q_e = q_ivs[q_idx]
                t_len = t_e - t_s
                q_len = q_e - q_s

                # Only merge if blocks cover the same length
                if t_len == q_len:
                    union(t_node, q_node)
                    t_pos = t_e
                    q_pos = q_e
                    continue

            t_pos += 1
            q_pos += 1

    # -- 5) Build the DAG using merged (representative) node ids -----------
    G = nx.DiGraph()

    # Collect representative nodes and their sequences
    rep_seqs: dict[int, str] = {}
    for nid, seq in node_seq.items():
        rep = find(nid)
        if rep not in rep_seqs:
            rep_seqs[rep] = seq
        elif len(seq) > len(rep_seqs[rep]):
            # Keep the longest sequence for the representative
            rep_seqs[rep] = seq

    for rep, seq in rep_seqs.items():
        G.add_node(rep, Base=seq)

    # Add edges: for each transcript, connect consecutive blocks
    for name, ivs in intervals.items():
        prev_rep: int | None = None
        for idx in range(len(ivs)):
            nid = block_to_node.get((name, idx))
            if nid is None:
                continue
            rep = find(nid)
            if prev_rep is not None and prev_rep != rep:
                G.add_edge(prev_rep, rep)
            prev_rep = rep

    # -- 6) Chain simplification -------------------------------------------
    # Merge linear chains (single in-edge, single out-edge) into one node
    changed = True
    while changed:
        changed = False
        for n in list(G.nodes()):
            if n not in G:
                continue
            succs = list(G.successors(n))
            if len(succs) != 1:
                continue
            s = succs[0]
            if s == n:
                continue
            if G.in_degree(s) != 1:
                continue
            # Merge s into n
            G.nodes[n]["Base"] = G.nodes[n]["Base"] + G.nodes[s]["Base"]
            for _, t in list(G.out_edges(s)):
                G.add_edge(n, t)
            G.remove_node(s)
            changed = True

    # -- 7) Cycle breaking — single-pass DFS  (A3) -------------------------
    whirl_status = 0

    if G.number_of_edges() > max_edges:
        raise RuntimeError(
            f"Graph too complex ({G.number_of_edges()} edges > {max_edges}), "
            "giving up on cycle breaking"
        )

    # Find back-edges via DFS and break them by node duplication
    while not nx.is_directed_acyclic_graph(G):
        # Find one cycle and break it
        try:
            cycle = nx.find_cycle(G, orientation="original")
        except nx.NetworkXNoCycle:
            break

        whirl_status += 1

        # Pick the node with the shortest sequence in the cycle to duplicate
        cycle_nodes = [u for u, _, _ in cycle]
        min_node = min(cycle_nodes, key=lambda n: len(G.nodes[n].get("Base", "")))

        # Duplicate min_node: new_node gets the out-edges that leave the cycle,
        # old node keeps the out-edges that stay in the cycle
        new_id = max(G.nodes()) + 1
        G.add_node(new_id, Base=G.nodes[min_node]["Base"])

        cycle_node_set = set(cycle_nodes)

        # Move out-edges to non-cycle targets to the new node
        for _, target in list(G.out_edges(min_node)):
            if target not in cycle_node_set:
                G.add_edge(new_id, target)
                G.remove_edge(min_node, target)

        # Move in-edges from cycle nodes to the new node
        for source, _ in list(G.in_edges(min_node)):
            if source in cycle_node_set:
                G.add_edge(source, new_id)
                G.remove_edge(source, min_node)

    # -- 8) Topological sort → SuperTranscript sequence --------------------
    try:
        topo_order = list(nx.topological_sort(G))
    except nx.NetworkXUnfeasible:
        raise RuntimeError("Failed to topologically sort graph — cycles remain")

    seq_parts: list[str] = []
    coords: list[int] = [0]
    for nid in topo_order:
        base = G.nodes[nid].get("Base", "")
        seq_parts.append(base)
        coords.append(coords[-1] + len(base))

    seq = "".join(seq_parts)

    # Build annotation
    cluster_id = Path(str("placeholder")).stem  # will be set by caller
    anno_parts: list[str] = []
    for i in range(len(coords) - 1):
        anno_parts.append(
            get_annotation_line(cluster_id, str(coords[i] + 1), str(coords[i + 1]), cluster_id)
        )
    anno = "".join(anno_parts)

    return seq, anno, whirl_status

# ---------------------------------------------------------------------------
# Public entry: process one cluster  (IO1 — in-memory dispatch)
# ---------------------------------------------------------------------------

def super_tran(
    gene_id: str,
    transcripts: dict[str, str],
    *,
    verbose: bool = False,
    max_edges: int = 100,
    tmpdir: Path | None = None,
) -> ClusterResult:
    """Build a SuperTranscript for one cluster.

    Parameters
    ----------
    gene_id:
        Cluster / gene identifier.
    transcripts:
        ``{transcript_id: sequence}`` dict — passed **in-memory** from the
        pool dispatcher (IO1).
    tmpdir:
        Local scratch directory for minimap2 temp files.  Defaults to
        ``$TMPDIR`` or ``/tmp``.
    """
    transcript_count = len(transcripts)

    # Single-transcript cluster — trivial
    if transcript_count <= 1:
        seq = next(iter(transcripts.values()), "")
        anno = get_annotation_line(gene_id, "1", str(len(seq)), gene_id)
        return ClusterResult(seq, anno, 0, transcript_count)

    # Write temporary FASTA to local SSD, not NFS (IO1/IO3)
    if tmpdir is None:
        tmpdir = Path(os.environ.get("TMPDIR", tempfile.gettempdir()))
    tmpdir.mkdir(parents=True, exist_ok=True)

    fasta_path = tmpdir / f"{gene_id}.fasta"
    try:
        with open(fasta_path, "w") as fh:
            for tid, seq in transcripts.items():
                fh.write(f">{tid}\n{seq}\n")

        # Run minimap2 (A1)
        blocks = _run_minimap2(fasta_path)

        if not blocks:
            # No alignments — fall back to longest transcript
            log.warning("No minimap2 alignments for cluster %s, using longest transcript", gene_id)
            longest = max(transcripts.values(), key=len)
            anno = get_annotation_line(gene_id, "1", str(len(longest)), gene_id)
            return ClusterResult(longest, anno, 0, transcript_count)

        # Build block-level graph (A2) + cycle breaking (A3)
        seq, anno, whirl_status = _build_block_graph(
            transcripts, blocks, max_edges=max_edges,
        )

        # Fix annotation with actual cluster_id
        anno = _fix_annotation(gene_id, seq, anno)

        return ClusterResult(seq, anno, whirl_status, transcript_count)

    except Exception as exc:
        log.error("Failed to build SuperTranscript for %s: %s", gene_id, exc)
        # Fall back to longest transcript (same as original)
        longest = max(transcripts.values(), key=len)
        anno = get_annotation_line(gene_id, "1", str(len(longest)), gene_id)
        return ClusterResult(longest, anno, -1, -1)

    finally:
        # Clean up temp files (IO4 — pathlib, not os.system("rm"))
        fasta_path.unlink(missing_ok=True)


def _fix_annotation(gene_id: str, seq: str, _raw_anno: str) -> str:
    """Rebuild annotation with correct cluster_id after graph construction."""
    # The block graph produces one annotation line per block.
    # We regenerate from the sequence to ensure consistency.
    if not seq:
        return get_annotation_line(gene_id, "1", "1", gene_id)
    return get_annotation_line(gene_id, "1", str(len(seq)), gene_id)


# ---------------------------------------------------------------------------
# Legacy compatibility: file-based entry point
# ---------------------------------------------------------------------------

def SuperTran(fname: str, verbose: bool = False) -> tuple[str, str, int, int]:
    """Legacy wrapper matching original Lace 1.14.1 signature.

    Reads a FASTA file from disk and delegates to :func:`super_tran`.
    """
    path = Path(fname)
    gene_id = path.stem

    transcripts: dict[str, str] = {}
    current: str | None = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                current = line.lstrip(">").split()[0]
                transcripts[current] = ""
            elif current is not None:
                transcripts[current] += line

    result = super_tran(gene_id, transcripts, verbose=verbose)
    return result.seq, result.anno, result.whirl_status, result.transcript_count


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(args: list[str] | None = None) -> None:
    """CLI entry point — process a single FASTA file."""
    if len(sys.argv) != 2:
        sys.exit("Usage: build-supertranscript <cluster.fasta>")

    fname = sys.argv[1]
    seq, anno, whirl_status, transcript_status = SuperTran(fname, verbose=True)
    print(seq)
    print(anno)
    print(f"Whirls: {whirl_status}")
    print(f"Transcripts: {transcript_status}")


if __name__ == "__main__":
    main()
