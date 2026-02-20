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

from Lace import __version__
from Lace.BuildSuperTranscript import (
    ClusterResult,
    get_annotation_line,
    super_tran,
)

log = logging.getLogger("lace")

# ---------------------------------------------------------------------------
# Banner
# ---------------------------------------------------------------------------

BANNER = f"""\
 __      __    ____  ____
(  )    / _\\  /    )(  __)
/  (_/\\/    \\(  (__  ) _)
\\_____/\\_/\\_/\\_____)(____)\u0020
Version {__version__}  (minimap2 \u00b7 block-graph \u00b7 Python {sys.version_info.major}.{sys.version_info.minor})
"""

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
# Core pipeline
# ---------------------------------------------------------------------------

def split_and_build(
    genome_path: Path,
    corset_path: Path,
    n_cores: int,
    max_tran: int,
    out_dir: Path,
    *,
    tidy: bool = False,
) -> None:
    """Parse inputs, dispatch clusters to workers, write outputs."""
    start_time = time.time()

    # -- 1) Parse Corset cluster file -> transcript->cluster mapping -------
    cluster_map: dict[str, str] = {}
    if corset_path.is_file():
        log.info("Parsing cluster file %s", corset_path)
        with open(corset_path) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 2:
                    cluster_map[parts[0]] = parts[1].rstrip("\n")

    # Count transcripts per cluster to identify singletons
    cluster_counts: dict[str, int] = {}
    for clust in cluster_map.values():
        cluster_counts[clust] = cluster_counts.get(clust, 0) + 1
    single_clusters = {c for c, n in cluster_counts.items() if n == 1}

    # -- 2) Parse FASTA file -----------------------------------------------
    log.info("Parsing transcripts from %s", genome_path)
    transcripts: dict[str, str] = {}
    gene_of: dict[str, str] = {}
    current: str | None = None

    with open(genome_path) as fh:
        for line in fh:
            if line.startswith(">"):
                current = line.split()[0].lstrip(">")
                transcripts[current] = ""
                clust = cluster_map.get(current)
                if clust and clust not in single_clusters:
                    gene_of[current] = clust
                else:
                    gene_of[current] = "None"
            elif current is not None:
                transcripts[current] += line.strip()

    # -- 3) Group transcripts by gene/cluster (in-memory, IO1) -------------
    gene_transcripts: dict[str, dict[str, str]] = {}
    for tag, clust in gene_of.items():
        if clust == "None":
            continue
        if clust not in gene_transcripts:
            gene_transcripts[clust] = {}
        if len(gene_transcripts[clust]) < max_tran:
            gene_transcripts[clust][tag] = transcripts[tag]
        elif len(gene_transcripts[clust]) == max_tran:
            log.warning(
                "Cluster %s: capping at %d transcripts (has more)",
                clust, max_tran,
            )

    n_multi = len(gene_transcripts)
    log.info(
        "%d multi-transcript clusters, %d single-transcript clusters",
        n_multi,
        len(single_clusters),
    )

    # -- 4) Sort by transcript count descending (L4 -- largest first) ------
    sorted_genes = sorted(
        gene_transcripts.items(),
        key=lambda kv: len(kv[1]),
        reverse=True,
    )

    # -- 5) Dispatch to ProcessPoolExecutor (L3, IO1) ----------------------
    tmpdir = os.environ.get("TMPDIR", "/tmp")
    out_dir.mkdir(parents=True, exist_ok=True)
    max_edges = 100

    results_map: dict[str, ClusterResult] = {}
    gene_order: list[str] = [g for g, _ in sorted_genes]

    log.info(
        "Building SuperTranscripts with %d cores (minimap2 + block-graph)",
        n_cores,
    )

    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        futures = {
            executor.submit(
                _worker, gene_id, tdict, max_edges, tmpdir,
            ): gene_id
            for gene_id, tdict in sorted_genes
        }

        with tqdm(total=n_multi, desc="Building SuperTranscripts", unit="cluster") as pbar:
            for future in as_completed(futures):
                gene_id = futures[future]
                try:
                    _, result = future.result()
                    results_map[gene_id] = result
                except Exception as exc:
                    log.error("Cluster %s failed: %s", gene_id, exc)
                    results_map[gene_id] = ClusterResult("", "", -1, -1)
                pbar.update(1)

    # -- 6) Write SuperDuper.fasta and SuperDuper.gff ----------------------
    super_fasta = out_dir / "SuperDuper.fasta"
    super_gff = out_dir / "SuperDuper.gff"

    with (
        open(super_fasta, "w") as ff,
        open(super_gff, "w") as fg,
    ):
        # Multi-transcript clusters (in sorted order for reproducibility)
        for gene_id in gene_order:
            res = results_map.get(gene_id)
            if res is None:
                continue
            ff.write(
                f">{gene_id} NoTrans:{res.transcript_count},"
                f"Whirls:{res.whirl_status}\n"
            )
            ff.write(f"{res.seq}\n")
            fg.write(res.anno)

        # Single-transcript clusters
        for tag, clust in cluster_map.items():
            if clust not in single_clusters:
                continue
            seq = transcripts.get(tag, "")
            anno = get_annotation_line(clust, "1", str(len(seq)), tag)
            ff.write(f">{clust} NoTrans:1,Whirls:0\n")
            ff.write(f"{seq}\n")
            fg.write(anno)

    elapsed = time.time() - start_time
    log.info("BUILT SUPERTRANSCRIPTS ---- %.1f seconds ----", elapsed)

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
        help="Tab-delimited file mapping transcripts to clusters (e.g. Corset output)",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1)",
    )
    parser.add_argument(
        "--maxTran",
        type=int,
        default=50,
        help="Maximum transcripts per cluster (default: 50)",
    )
    parser.add_argument(
        "-o", "--outputDir",
        default=".",
        help="Output directory (default: .)",
    )
    parser.add_argument(
        "-t", "--tidy",
        action="store_true",
        help="Remove intermediate files after running",
    )
    parser.add_argument(
        "-a", "--alternate",
        action="store_true",
        help="Create alternate annotations and metrics (requires Checker)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
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

    split_and_build(
        genome_path=Path(parsed.TranscriptsFile),
        corset_path=Path(parsed.ClusterFile),
        n_cores=parsed.cores,
        max_tran=parsed.maxTran,
        out_dir=Path(parsed.outputDir),
        tidy=parsed.tidy,
    )

    if parsed.alternate:
        from Lace.Checker import Checker
        cwd = os.getcwd()
        os.chdir(parsed.outputDir)
        log.info("Making alternate annotation and checks")
        Checker("SuperDuper.fasta", "SuperDuper.gff", parsed.cores, "SuperFiles")
        os.chdir(cwd)

    log.info("Done")


if __name__ == "__main__":
    main()
