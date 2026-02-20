# lace-fast 2.0

High-performance SuperTranscript construction — rewrite of
[Lace 1.14.1](https://github.com/Oshlack/Lace).

## What changed

| Feature | Lace 1.14.1 | lace-fast 2.0 |
|---------|-------------|---------------|
| Aligner | BLAT v35 (academic licence) | minimap2 2.28 (MIT) |
| Graph | One node per **base** (~10 000 nodes) | One node per **block** (~10–30 nodes) |
| Cycle breaking | `nx.simple_cycles` in a while-loop (exponential) | DFS `nx.find_cycle` single-pass (O(V+E)) |
| I/O | ~376K NFS files per run | In-memory dispatch, `$TMPDIR` scratch |
| Parallelism | `multiprocessing.Pool` (chunksize=1) | `ProcessPoolExecutor` + largest-first scheduling |
| Progress | `print()` | `tqdm` + structured `logging` |
| Python | 3.10 | 3.13 |
| Matplotlib | Imported everywhere (unused) | Removed from dependencies |

## Expected performance

| Phase | BMAX (16 cores) | Speedup |
|-------|-----------------|---------|
| Baseline (v1.14.1) | 5 h 32 min | 1× |
| lace-fast 2.0 | **8–15 min** | **22–42×** |

## Installation

### Container (recommended)

```bash
# From the Lace repo root:
bash lace-fast/build.sh
# → produces lace-fast/lace_fast_2.0.sif
```

### pip

```bash
pip install ./lace-fast
```

Requires `minimap2` on `$PATH`.

## Usage

The CLI is drop-in compatible with the original Lace:

```bash
lace-fast transcripts.fasta clusters.txt --cores 16 -o output/
```

## Files

```
lace-fast/
├── pyproject.toml                    # Package metadata (PEP 621)
├── Lace_fast.def                     # Singularity/Apptainer container def
├── build.sh                          # Container build script
├── README.md                         # This file
└── lace_fast/
    ├── __init__.py                   # Version string
    ├── build_supertranscript.py      # Core: minimap2 + block graph + DFS cycles
    └── run.py                        # Entry point: pool dispatch + progress
```

## Tier implementation status

- [x] **Tier 0** — Hygiene (H1–H14): subprocess, pathlib, typing, f-strings,
  vectorised reverse-complement, modern exception handling
- [x] **Tier 1** — I/O Elimination (IO1–IO4): in-memory dispatch, $TMPDIR
  scratch, no NFS file storm
- [x] **Tier 2** — Algorithmic (A1–A3): minimap2, block-level graph, DFS
  cycle breaking
- [x] **Tier 3** — Architecture (L3–L5): ProcessPoolExecutor, largest-first
  scheduling, tqdm progress
- [x] **Tier 4** — Container (C1–C5): Python 3.13, minimap2, no matplotlib

## Validation

Compare against Lace 1.14.1 output:

1. **Sequence diff**: `SuperDuper.fasta` header-by-header comparison
2. **Length correlation**: Pearson r > 0.99 expected
3. **Downstream**: BUSCO score must not decrease
4. **Whirl stats**: DFS cycle-breaking may differ in count but produces valid DAGs
