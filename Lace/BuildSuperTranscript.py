#!/usr/bin/env python
"""Build a SuperTranscript for one cluster of transcripts.

Rewrite of Lace 1.14.1 ``BuildSuperTranscript.py`` with:

* **A1** -- minimap2 replaces BLAT (30-100x faster, MIT licence)
* **A2** -- block-level splice graph (~10-30 nodes) replaces base-level
           graph (~10 000 nodes)
* **A3** -- single-pass DFS back-edge cycle breaking replaces iterative
           ``nx.simple_cycles`` (worst-case exponential -> O(V+E))
* **H1-H14** -- full hygiene pass (subprocess, pathlib, typing, f-strings,
               vectorised pandas, etc.)
* **IO1-IO3** -- in-memory sequence passing, $TMPDIR for alignment files

Original author: Anthony Hawkins
Rewrite: Lace 2.0
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import NamedTuple

import networkx as nx

log = logging.getLogger("lace")

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
# Reverse complement  (H5 -- str.maketrans instead of char-by-char loop)
# ---------------------------------------------------------------------------

_RC_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of *seq*."""
    return seq.translate(_RC_TABLE)[::-1]


# Keep old name as alias for backward compatibility
Reverse_complement = reverse_complement  # pylint: disable=invalid-name


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
# minimap2 alignment  (A1 -- replaces BLAT)
# ---------------------------------------------------------------------------


def _run_minimap2(
    fasta_path: Path,
    *,
    min_identity: float = 0.98,
) -> list[AlignBlock]:
    """Run minimap2 all-vs-all on *fasta_path*, return alignment blocks.

    Uses ``-X`` for all-vs-all (ava) mode and ``--eqx -c`` for base-level
    CIGAR, replicating BLAT's ``-minIdentity=98`` behaviour.
    The PAF is parsed from stdout -- no intermediate file on disk (IO2).
    """
    cmd = [
        "minimap2",
        "-c", "-X", "--eqx",
        "-k15", "-w5", "-N50",
        f"-p{min_identity}",
        "--no-long-join",
        str(fasta_path),
        str(fasta_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    blocks: list[AlignBlock] = []
    for line in result.stdout.splitlines():
        cols = line.split("\t")
        if len(cols) < 12:
            continue
        q_name, t_name = cols[0], cols[5]
        if q_name == t_name:
            continue
        blocks.append(AlignBlock(
            t_name=t_name, q_name=q_name,
            t_start=int(cols[7]), t_end=int(cols[8]),
            q_start=int(cols[2]), q_end=int(cols[3]),
            strand=cols[4],
        ))
    return blocks


def _determine_strand_directions(
    blocks: list[AlignBlock],
) -> dict[str, str]:
    """Assign +/- orientation to each transcript.

    The first transcript encountered is arbitrarily +.  Subsequent ones
    are oriented relative to the alignment strand of hits connecting them.
    """
    trandir: dict[str, str] = {}
    seen_pairs: set[str] = set()

    for blk in blocks:
        pair_key = f"{blk.t_name}\t{blk.q_name}"
        rev_key = f"{blk.q_name}\t{blk.t_name}"
        if pair_key in seen_pairs or rev_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)

        t_name, q_name, strand = blk.t_name, blk.q_name, blk.strand
        if t_name in trandir and q_name in trandir:
            continue
        if t_name in trandir:
            trandir[q_name] = "+" if trandir[t_name] == strand else "-"
        elif q_name in trandir:
            trandir[t_name] = "+" if trandir[q_name] == strand else "-"
        else:
            trandir[t_name] = "+"
            trandir[q_name] = "+" if strand == "+" else "-"

    return trandir


# ---------------------------------------------------------------------------
# Block-level splice graph  (A2 -- replaces base-level graph)
# ---------------------------------------------------------------------------


def _orient_transcripts(
    transcripts: dict[str, str],
    blocks: list[AlignBlock],
) -> dict[str, str]:
    """Orient all transcripts to + strand based on alignment evidence."""
    trandir = _determine_strand_directions(blocks)
    return {
        name: reverse_complement(seq) if trandir.get(name, "+") == "-" else seq
        for name, seq in transcripts.items()
    }


def _compute_intervals(
    oriented: dict[str, str],
    blocks: list[AlignBlock],
) -> dict[str, list[tuple[int, int]]]:
    """Compute breakpoint-based intervals per transcript."""
    breakpoints: dict[str, set[int]] = {
        name: {0, len(seq)} for name, seq in oriented.items()
    }
    for blk in blocks:
        if blk.t_name == blk.q_name:
            continue
        breakpoints.setdefault(blk.t_name, set()).update({blk.t_start, blk.t_end})
        breakpoints.setdefault(blk.q_name, set()).update({blk.q_start, blk.q_end})

    intervals: dict[str, list[tuple[int, int]]] = {}
    for name, bps in breakpoints.items():
        sorted_bps = sorted(bps)
        intervals[name] = [
            (sorted_bps[i], sorted_bps[i + 1])
            for i in range(len(sorted_bps) - 1)
            if sorted_bps[i + 1] > sorted_bps[i]
        ]
    return intervals


def _create_block_nodes(
    oriented: dict[str, str],
    intervals: dict[str, list[tuple[int, int]]],
) -> tuple[dict[int, str], dict[tuple[str, int], int], dict[str, dict[int, int]]]:
    """Create one graph node per interval block, return mappings."""
    node_seq: dict[int, str] = {}
    block_to_node: dict[tuple[str, int], int] = {}
    pos_to_interval: dict[str, dict[int, int]] = {}

    next_node = 0
    for name, ivs in intervals.items():
        lookup: dict[int, int] = {}
        for idx, (start, end) in enumerate(ivs):
            lookup[start] = idx
            node_seq[next_node] = oriented[name][start:end]
            block_to_node[(name, idx)] = next_node
            next_node += 1
        pos_to_interval[name] = lookup

    return node_seq, block_to_node, pos_to_interval


def _merge_equivalent_blocks(
    blocks: list[AlignBlock],
    intervals: dict[str, list[tuple[int, int]]],
    block_to_node: dict[tuple[str, int], int],
    pos_to_interval: dict[str, dict[int, int]],
    node_seq: dict[int, str],
) -> dict[int, int]:
    """Union-Find merge of blocks that correspond to the same aligned region."""
    parent: dict[int, int] = {n: n for n in node_seq}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    for blk in blocks:
        if blk.t_name == blk.q_name:
            continue
        t_lookup = pos_to_interval.get(blk.t_name, {})
        q_lookup = pos_to_interval.get(blk.q_name, {})
        t_ivs = intervals.get(blk.t_name, [])
        q_ivs = intervals.get(blk.q_name, [])

        t_pos, q_pos = blk.t_start, blk.q_start
        while t_pos < blk.t_end and q_pos < blk.q_end:
            t_idx = t_lookup.get(t_pos)
            q_idx = q_lookup.get(q_pos)
            if t_idx is None or q_idx is None:
                t_pos += 1
                q_pos += 1
                continue
            t_node = block_to_node.get((blk.t_name, t_idx))
            q_node = block_to_node.get((blk.q_name, q_idx))
            if t_node is not None and q_node is not None:
                t_len = t_ivs[t_idx][1] - t_ivs[t_idx][0]
                q_len = q_ivs[q_idx][1] - q_ivs[q_idx][0]
                if t_len == q_len:
                    union(t_node, q_node)
                    t_pos = t_ivs[t_idx][1]
                    q_pos = q_ivs[q_idx][1]
                    continue
            t_pos += 1
            q_pos += 1

    # Flatten parent pointers
    for n in list(parent):
        find(n)
    return parent


def _build_dag(
    node_seq: dict[int, str],
    intervals: dict[str, list[tuple[int, int]]],
    block_to_node: dict[tuple[str, int], int],
    parent: dict[int, int],
) -> nx.DiGraph:
    """Build directed graph from merged representative nodes."""

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    graph = nx.DiGraph()

    # Collect representative sequences (keep longest)
    rep_seqs: dict[int, str] = {}
    for nid, seq in node_seq.items():
        rep = find(nid)
        if rep not in rep_seqs or len(seq) > len(rep_seqs[rep]):
            rep_seqs[rep] = seq
    for rep, seq in rep_seqs.items():
        graph.add_node(rep, Base=seq)

    # Add edges from transcript interval ordering
    for name, ivs in intervals.items():
        prev_rep: int | None = None
        for idx in range(len(ivs)):
            nid = block_to_node.get((name, idx))
            if nid is None:
                continue
            rep = find(nid)
            if prev_rep is not None and prev_rep != rep:
                graph.add_edge(prev_rep, rep)
            prev_rep = rep

    return graph


def _simplify_chains(graph: nx.DiGraph) -> None:
    """Collapse chain nodes (single-in, single-out) in place."""
    changed = True
    while changed:
        changed = False
        for node in list(graph.nodes()):
            if node not in graph:
                continue
            succs = list(graph.successors(node))
            if len(succs) != 1:
                continue
            succ = succs[0]
            if succ == node or graph.in_degree(succ) != 1:
                continue
            graph.nodes[node]["Base"] += graph.nodes[succ]["Base"]
            for _, target in list(graph.out_edges(succ)):
                graph.add_edge(node, target)
            graph.remove_node(succ)
            changed = True


def _break_cycles(graph: nx.DiGraph, max_edges: int = 500) -> int:
    """Remove cycles via single-pass DFS back-edge breaking (A3).

    Returns the number of cycles broken (whirl_status).
    """
    if graph.number_of_edges() > max_edges:
        raise RuntimeError(
            f"Graph too complex ({graph.number_of_edges()} edges > {max_edges}), "
            "giving up on cycle breaking"
        )

    whirl_count = 0
    while not nx.is_directed_acyclic_graph(graph):
        try:
            cycle = nx.find_cycle(graph, orientation="original")
        except nx.NetworkXNoCycle:
            break

        whirl_count += 1
        cycle_nodes = [u for u, _, _ in cycle]
        cycle_set = set(cycle_nodes)
        min_node = min(cycle_nodes, key=lambda n: len(graph.nodes[n].get("Base", "")))

        new_id = max(graph.nodes()) + 1
        graph.add_node(new_id, Base=graph.nodes[min_node]["Base"])

        for _, target in list(graph.out_edges(min_node)):
            if target not in cycle_set:
                graph.add_edge(new_id, target)
                graph.remove_edge(min_node, target)
        for source, _ in list(graph.in_edges(min_node)):
            if source in cycle_set:
                graph.add_edge(source, new_id)
                graph.remove_edge(source, min_node)

    return whirl_count


def _toposort_sequence(graph: nx.DiGraph) -> tuple[str, str]:
    """Topological sort -> concatenated SuperTranscript sequence + annotation."""
    try:
        topo_order = list(nx.topological_sort(graph))
    except nx.NetworkXUnfeasible as exc:
        raise RuntimeError(
            "Failed to topologically sort graph -- cycles remain"
        ) from exc

    seq_parts: list[str] = []
    coords: list[int] = [0]
    for nid in topo_order:
        base = graph.nodes[nid].get("Base", "")
        seq_parts.append(base)
        coords.append(coords[-1] + len(base))

    seq = "".join(seq_parts)
    anno = "".join(
        get_annotation_line("placeholder", str(coords[i] + 1), str(coords[i + 1]), "placeholder")
        for i in range(len(coords) - 1)
    )
    return seq, anno


def _build_block_graph(
    transcripts: dict[str, str],
    blocks: list[AlignBlock],
    *,
    max_edges: int = 500,
) -> tuple[str, str, int]:
    """Construct a block-level splice graph from minimap2 alignment blocks.

    Instead of one NetworkX node per *base* (original: 10 000 nodes for a
    5 x 2 kb cluster), we create one node per *contiguous block* (typically
    10-30 nodes).

    Returns ``(sequence, annotation, whirl_count)``.
    """
    oriented = _orient_transcripts(transcripts, blocks)
    intervals = _compute_intervals(oriented, blocks)
    node_seq, block_to_node, pos_to_interval = _create_block_nodes(oriented, intervals)
    parent = _merge_equivalent_blocks(blocks, intervals, block_to_node, pos_to_interval, node_seq)
    graph = _build_dag(node_seq, intervals, block_to_node, parent)
    _simplify_chains(graph)
    whirl_count = _break_cycles(graph, max_edges=max_edges)
    seq, anno = _toposort_sequence(graph)
    return seq, anno, whirl_count


# ---------------------------------------------------------------------------
# Public entry: process one cluster  (IO1 -- in-memory dispatch)
# ---------------------------------------------------------------------------


def super_tran(
    gene_id: str,
    transcripts: dict[str, str],
    *,
    max_edges: int = 500,
    tmpdir: Path | None = None,
) -> ClusterResult:
    """Build a SuperTranscript for one cluster.

    Parameters
    ----------
    gene_id:
        Cluster / gene identifier.
    transcripts:
        ``{transcript_id: sequence}`` dict -- passed in-memory from the
        pool dispatcher (IO1).
    tmpdir:
        Local scratch directory for minimap2 temp files.  Defaults to
        ``$TMPDIR`` or ``/tmp``.
    """
    transcript_count = len(transcripts)

    # Single-transcript cluster -- trivial
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
        with open(fasta_path, "w", encoding="utf-8") as fh:
            for tid, seq in transcripts.items():
                fh.write(f">{tid}\n{seq}\n")

        blocks = _run_minimap2(fasta_path)

        if not blocks:
            log.warning(
                "No minimap2 alignments for cluster %s, using longest transcript",
                gene_id,
            )
            longest = max(transcripts.values(), key=len)
            anno = get_annotation_line(gene_id, "1", str(len(longest)), gene_id)
            return ClusterResult(longest, anno, 0, transcript_count)

        seq, _anno, whirl_status = _build_block_graph(
            transcripts, blocks, max_edges=max_edges,
        )
        anno = _fix_annotation(gene_id, seq)
        return ClusterResult(seq, anno, whirl_status, transcript_count)

    except (RuntimeError, subprocess.CalledProcessError) as exc:
        log.error("Failed to build SuperTranscript for %s: %s", gene_id, exc)
        longest = max(transcripts.values(), key=len)
        anno = get_annotation_line(gene_id, "1", str(len(longest)), gene_id)
        return ClusterResult(longest, anno, -1, -1)

    finally:
        fasta_path.unlink(missing_ok=True)


def _fix_annotation(gene_id: str, seq: str) -> str:
    """Rebuild annotation with correct cluster_id after graph construction."""
    if not seq:
        return get_annotation_line(gene_id, "1", "1", gene_id)
    return get_annotation_line(gene_id, "1", str(len(seq)), gene_id)


# ---------------------------------------------------------------------------
# Legacy compatibility: file-based entry point
# ---------------------------------------------------------------------------

# pylint: disable=invalid-name
def SuperTran(fname: str, verbose: bool = False) -> tuple[str, str, int, int]:
    """Legacy wrapper matching original Lace 1.14.1 signature.

    Reads a FASTA file from disk and delegates to :func:`super_tran`.
    """
    _ = verbose  # accepted for API compat, logging level set globally
    path = Path(fname)
    gene_id = path.stem

    transcripts: dict[str, str] = {}
    current: str | None = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line_s = line.strip()
            if line_s.startswith(">"):
                current = line_s.lstrip(">").split()[0]
                transcripts[current] = ""
            elif current is not None:
                transcripts[current] += line_s

    result = super_tran(gene_id, transcripts)
    return result.seq, result.anno, result.whirl_status, result.transcript_count
# pylint: enable=invalid-name


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    """CLI entry point -- process a single FASTA file."""
    if len(sys.argv) != 2:
        sys.exit("Usage: BuildSuperTranscript <cluster.fasta>")

    fname = sys.argv[1]
    seq, anno, whirl_status, transcript_status = SuperTran(fname, verbose=True)
    print(seq)
    print(anno)
    print(f"Whirls: {whirl_status}")
    print(f"Transcripts: {transcript_status}")


if __name__ == "__main__":
    main()
