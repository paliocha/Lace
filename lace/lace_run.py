#!/usr/bin/env python
"""Lace entry point -- parallelised SuperTranscript construction.

Rewrite of Lace 1.14.1 ``Lace_run.py`` with:

* **H1-H14** -- full hygiene pass (subprocess, pathlib, typing, f-strings,
               structured logging, modern exception handling)
* **IO1**    -- in-memory cluster dispatch via ``ProcessPoolExecutor``
               (eliminates ~376K NFS file operations)
* **L3**     -- ``concurrent.futures.ProcessPoolExecutor`` replaces
               ``multiprocessing.Pool`` for cleaner exception handling
* **L4**     -- largest-first scheduling for natural work-stealing
* **L5**     -- ``tqdm`` progress bar + structured ``logging``

Original author: Anthony Hawkins
Rewrite: Lace 2.0
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from tqdm import tqdm

from lace import __version__
from lace.build_supertranscript import (
    ClusterResult,
    get_annotation_line,
    super_tran,
)

log = logging.getLogger("lace")

# ---------------------------------------------------------------------------
# Banner -- resembles the Lace logo (double-helix + "Lace" in Century Gothic)
# https://github.com/Oshlack/Lace/blob/master/WikiFigs/logo.png
# Shading: + = solid, . = edge/anti-alias
# ---------------------------------------------------------------------------

BANNER = (
    "\n"
    "                      .++++++.\n"
    "                  .++++++++++++.\n"
    "                .++++++.  .++++++.\n"
    "              .++++.       .++++++.                .+++.\n"
    "             .++++.         .+++++.            .+++++++++.\n"
    "             +++++.         .+++++.          .+++++++++++++.\n"
    "            .+++++.        .+++++.          .+++++. .+++++.\n"
    "            .++++++.      .+++++.          .+++++. .+++++.        .+++++.\n"
    "             .+++++++.  .+++++.            .++++++++++.+.      .+++++++++++.\n"
    "              .++++++++++++++.              .+++++++++.        .++++. .++++.\n"
    "                .++++++++++.                 .+++++++++++.     .++++++++++.\n"
    "                 ...++++++++++++.              ...+++++++++++++++++.+++++++++.+.\n"
    "                .+++++..++++++++++++++++++++++.. .+++++++++.      .+++++++.\n"
    "              .++++.    .+++++++++++++++.\n"
    "         .+++++++.           .+++++.\n"
    "     .++++++++.\n"
    "  .++++++.\n"
    "\n"
    "  ++\n"
    "  ++\n"
    "  ++           .++++++.        .+++++.        .++++++.\n"
    "  ++         .++.    .++.    .++.   .++.    .++.    .++.\n"
    "  ++         ++.      .++    ++.       .    .++      .++\n"
    "  ++         ++.      .++    ++.            .++++++++++.\n"
    "  ++         ++.      .++    ++.            .++.\n"
    "  ++         .++.    .+++    .++.   .++.    .++.    .+.\n"
    "  +++++++++    .++++++.++      .+++++.        .++++++.\n"
    "\n"
    "  -----------------------------------------------------\n"
    f"  Version {__version__}  minimap2 \u00b7 block-graph \u00b7 Python "
    f"{sys.version_info.major}.{sys.version_info.minor}\n"
)


# ---------------------------------------------------------------------------
# Worker -- called in child processes  (IO1 -- receives data in-memory)
# ---------------------------------------------------------------------------


def _worker(
    gene_id: str,
    transcripts: dict[str, str],
    max_edges: int,
    tmpdir: str,
) -> tuple[str, ClusterResult]:
    """Process one cluster entirely in-memory.

    Receives the transcript dict via IPC pickle (cheaper than NFS I/O).
    Returns ``(gene_id, ClusterResult)``.
    """
    result = super_tran(
        gene_id,
        transcripts,
        max_edges=max_edges,
        tmpdir=Path(tmpdir),
    )
    return gene_id, result


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------


def _parse_clusters(corset_path: Path) -> dict[str, str]:
    """Parse Corset cluster file -> {transcript_id: cluster_id}."""
    cluster_map: dict[str, str] = {}
    if not corset_path.is_file():
        return cluster_map
    log.info("Parsing cluster file %s", corset_path)
    with open(corset_path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                cluster_map[parts[0]] = parts[1].rstrip("\n")
    return cluster_map


def _parse_transcripts(
    genome_path: Path,
    cluster_map: dict[str, str],
    single_clusters: set[str],
) -> tuple[dict[str, str], dict[str, str]]:
    """Parse FASTA -> (all_transcripts, gene_of) dicts."""
    log.info("Parsing transcripts from %s", genome_path)
    transcripts: dict[str, str] = {}
    gene_of: dict[str, str] = {}
    current: str | None = None

    with open(genome_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                current = line.split()[0].lstrip(">")
                transcripts[current] = ""
                clust = cluster_map.get(current)
                gene_of[current] = clust if (clust and clust not in single_clusters) else "None"
            elif current is not None:
                transcripts[current] += line.strip()

    return transcripts, gene_of


def _group_by_cluster(
    gene_of: dict[str, str],
    transcripts: dict[str, str],
    max_tran: int,
) -> tuple[dict[str, dict[str, str]], int]:
    """Group transcripts by cluster, capping at *max_tran*.

    Returns ``(gene_transcripts, n_capped)``.
    """
    log.info("Grouping transcripts by cluster...")
    gene_transcripts: dict[str, dict[str, str]] = {}
    n_capped = 0

    for tag, clust in gene_of.items():
        if clust == "None":
            continue
        if clust not in gene_transcripts:
            gene_transcripts[clust] = {}
        if len(gene_transcripts[clust]) < max_tran:
            gene_transcripts[clust][tag] = transcripts[tag]
        elif len(gene_transcripts[clust]) == max_tran:
            n_capped += 1
            log.debug(
                "Cluster %s capped at %d transcripts",
                clust, max_tran,
            )
    return gene_transcripts, n_capped


# ---------------------------------------------------------------------------
# Core pipeline
# ---------------------------------------------------------------------------


def _log_input_summary(
    n_transcripts: int,
    total_bases: int,
    n_total_clusters: int,
    n_multi: int,
    n_singletons: int,
    max_tran: int,
    n_capped: int,
) -> None:
    """Print the input-side statistics table."""
    avg_len = total_bases / n_transcripts if n_transcripts else 0
    log.info("")
    log.info("  %-34s %s", "Transcripts in FASTA:", f"{n_transcripts:,}")
    log.info("  %-34s %s bp", "  Total bases:", f"{total_bases:,}")
    log.info("  %-34s %s bp", "  Avg transcript length:", f"{avg_len:,.0f}")
    log.info("  %-34s %s", "Corset clusters (total):", f"{n_total_clusters:,}")
    log.info("  %-34s %s", "  Multi-transcript clusters:", f"{n_multi:,}")
    log.info("  %-34s %s", "  Singleton clusters:", f"{n_singletons:,}")
    log.info("  %-34s %d", "  Max transcripts/cluster (--maxTran):", max_tran)
    if n_capped:
        log.warning(
            "  %d cluster(s) hit the --maxTran=%d cap "
            "(only first %d transcripts kept)",
            n_capped, max_tran, max_tran,
        )
    log.info("")


def _dispatch_workers(
    sorted_genes: list[tuple[str, dict[str, str]]],
    n_cores: int,
) -> tuple[dict[str, ClusterResult], int]:
    """Run SuperTranscript construction across worker processes.

    Returns ``(results_map, n_failed)``.
    """
    tmpdir = os.environ.get("TMPDIR", "/tmp")
    max_edges = 500
    n_multi = len(sorted_genes)
    results_map: dict[str, ClusterResult] = {}
    n_failed = 0

    log.info(
        "Building %s multi-transcript clusters with %d workers "
        "(minimap2 + block-graph)...",
        f"{n_multi:,}", n_cores,
    )

    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        futures = {
            executor.submit(_worker, gene_id, tdict, max_edges, tmpdir): gene_id
            for gene_id, tdict in sorted_genes
        }

        with tqdm(
            total=n_multi,
            desc="Building SuperTranscripts",
            unit="cluster",
            file=sys.stderr,
        ) as pbar:
            for future in as_completed(futures):
                gene_id = futures[future]
                try:
                    _, result = future.result()
                    results_map[gene_id] = result
                except (RuntimeError, OSError) as exc:
                    log.error(
                        "Processing cluster %s: FAILED to construct: %s",
                        gene_id, exc,
                    )
                    results_map[gene_id] = ClusterResult("", "", -1, -1)
                    n_failed += 1
                pbar.update(1)

    if n_failed:
        log.warning("%d cluster(s) failed to construct", n_failed)

    return results_map, n_failed


def _log_output_summary(
    n_written: int,
    n_assembled: int,
    n_singletons: int,
    n_failed: int,
    n_transcripts: int,
    total_out_bases: int,
    elapsed: float,
) -> None:
    """Print the output-side statistics table."""
    reduction_pct = (1 - n_written / n_transcripts) * 100 if n_transcripts else 0
    avg_st_len = total_out_bases / n_written if n_written else 0
    log.info("")
    log.info("  %-34s %s", "SuperTranscripts written:", f"{n_written:,}")
    log.info("  %-34s %s", "  Assembled (multi-transcript):", f"{n_assembled:,}")
    log.info("  %-34s %s", "  Pass-through (singletons):", f"{n_singletons:,}")
    if n_failed:
        log.warning("  %-34s %d", "  Failed:", n_failed)
    log.info("  %-34s %s bp", "  Total bases:", f"{total_out_bases:,}")
    log.info("  %-34s %s bp", "  Avg SuperTranscript length:", f"{avg_st_len:,.0f}")
    log.info(
        "  %-34s %.1f%% (%s \u2192 %s sequences)",
        "Redundancy reduction:",
        reduction_pct, f"{n_transcripts:,}", f"{n_written:,}",
    )
    log.info("")
    log.info("BUILT SUPERTRANSCRIPTS ---- %.1f seconds ----", elapsed)
    log.info("Done")


def split_and_build(
    genome_path: Path,
    corset_path: Path,
    n_cores: int,
    max_tran: int,
    out_dir: Path,
) -> None:
    """Parse inputs, dispatch clusters to workers, write outputs."""
    start_time = time.time()

    # -- 1) Parse Corset cluster file --------------------------------------
    log.info("Parsing cluster assignments...")
    cluster_map = _parse_clusters(corset_path)

    cluster_counts: dict[str, int] = {}
    for clust in cluster_map.values():
        cluster_counts[clust] = cluster_counts.get(clust, 0) + 1
    single_clusters = {c for c, n in cluster_counts.items() if n == 1}

    # -- 2) Parse FASTA file -----------------------------------------------
    transcripts, gene_of = _parse_transcripts(genome_path, cluster_map, single_clusters)
    total_in_bases = sum(len(seq) for seq in transcripts.values())

    # -- 3) Group transcripts by cluster (in-memory, IO1) ------------------
    gene_transcripts, n_capped = _group_by_cluster(gene_of, transcripts, max_tran)
    n_multi = len(gene_transcripts)

    _log_input_summary(
        len(transcripts), total_in_bases, len(cluster_counts), n_multi,
        len(single_clusters), max_tran, n_capped,
    )

    # -- 4) Sort largest-first (L4) + dispatch -----------------------------
    sorted_genes = sorted(
        gene_transcripts.items(), key=lambda kv: len(kv[1]), reverse=True,
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    results_map, n_failed = _dispatch_workers(sorted_genes, n_cores)
    gene_order = [g for g, _ in sorted_genes]

    # -- 5) Write SuperDuper.fasta and SuperDuper.gff ----------------------
    n_written, total_out_bases = _write_outputs(
        out_dir, gene_order, results_map, cluster_map,
        single_clusters, transcripts,
    )

    # -- 6) Summary --------------------------------------------------------
    _log_output_summary(
        n_written, n_multi - n_failed, len(single_clusters),
        n_failed, len(transcripts),
        total_out_bases, time.time() - start_time,
    )


def _write_multi(fasta_fh, gff_fh, gene_order, results_map):
    """Write multi-transcript SuperTranscripts. Returns (count, bases)."""
    count = 0
    bases = 0
    for gene_id in gene_order:
        res = results_map.get(gene_id)
        if res is None or not res.seq:
            continue
        fasta_fh.write(
            f">{gene_id} NoTrans:{res.transcript_count},"
            f"Whirls:{res.whirl_status}\n"
        )
        fasta_fh.write(f"{res.seq}\n")
        gff_fh.write(res.anno)
        count += 1
        bases += len(res.seq)
    return count, bases


def _write_singles(fasta_fh, gff_fh, cluster_map, single_clusters, transcripts):
    """Write singleton pass-through SuperTranscripts. Returns (count, bases)."""
    count = 0
    bases = 0
    for tag, clust in cluster_map.items():
        if clust not in single_clusters:
            continue
        seq = transcripts.get(tag, "")
        anno = get_annotation_line(clust, "1", str(len(seq)), tag)
        fasta_fh.write(f">{clust} NoTrans:1,Whirls:0\n")
        fasta_fh.write(f"{seq}\n")
        gff_fh.write(anno)
        count += 1
        bases += len(seq)
    return count, bases


def _write_outputs(
    out_dir: Path,
    gene_order: list[str],
    results_map: dict[str, ClusterResult],
    cluster_map: dict[str, str],
    single_clusters: set[str],
    transcripts: dict[str, str],
) -> tuple[int, int]:
    """Write SuperDuper.fasta and SuperDuper.gff to *out_dir*.

    Returns ``(n_written, total_bases)``.
    """
    with (
        open(out_dir / "SuperDuper.fasta", "w", encoding="utf-8") as ff,
        open(out_dir / "SuperDuper.gff", "w", encoding="utf-8") as fg,
    ):
        n_multi, bases_multi = _write_multi(ff, fg, gene_order, results_map)
        n_singles, bases_singles = _write_singles(
            ff, fg, cluster_map, single_clusters, transcripts,
        )
    return n_multi + n_singles, bases_multi + bases_singles


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(args: list[str] | None = None) -> None:
    """Main entry point for Lace."""
    print(BANNER)

    parser = argparse.ArgumentParser(
        prog="Lace",
        description="Build SuperTranscripts from clustered transcript assemblies",
    )
    parser.add_argument(
        "TranscriptsFile",
        help="FASTA file containing all transcripts",
    )
    parser.add_argument(
        "ClusterFile",
        help="Tab-delimited file mapping transcripts to clusters "
             "(e.g. Corset output)",
    )
    parser.add_argument(
        "--cores", type=int, default=1,
        help="Number of parallel workers (default: 1)",
    )
    parser.add_argument(
        "--maxTran", type=int, default=50,
        help="Maximum transcripts per cluster (default: 50)",
    )
    parser.add_argument(
        "-o", "--outputDir", default=".",
        help="Output directory (default: .)",
    )
    parser.add_argument(
        "-t", "--tidy", action="store_true",
        help="Remove intermediate files after running",
    )
    parser.add_argument(
        "-a", "--alternate", action="store_true",
        help="Create alternate annotations and metrics (requires Checker)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable debug logging",
    )

    parsed = parser.parse_args(args)

    # Configure logging (L5)
    level = logging.DEBUG if parsed.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    log.info("Lace Version: %s", __version__)

    out_dir = Path(parsed.outputDir)
    if out_dir.exists():
        log.info("Output directory exists")
    else:
        log.info("Creating output directory")
        out_dir.mkdir(parents=True, exist_ok=True)

    split_and_build(
        genome_path=Path(parsed.TranscriptsFile),
        corset_path=Path(parsed.ClusterFile),
        n_cores=parsed.cores,
        max_tran=parsed.maxTran,
        out_dir=out_dir,
    )

    if parsed.alternate:
        from lace.checker import Checker  # pylint: disable=import-outside-toplevel
        cwd = os.getcwd()
        os.chdir(parsed.outputDir)
        log.info("Making alternate annotation and checks")
        Checker(
            "SuperDuper.fasta", "SuperDuper.gff",
            parsed.cores, "SuperFiles",
        )
        os.chdir(cwd)


if __name__ == "__main__":
    main()
