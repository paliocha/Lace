# Changelog

## [2.0.0] — 2026-02-20

### Complete rewrite of Lace

Lace 2.0 is a ground-up rewrite that replaces the BLAT-based alignment
pipeline with **minimap2** and introduces a **block-level graph** algorithm
for SuperTranscript construction.

### Highlights

- **minimap2** replaces BLAT — no licence restrictions, dramatically faster
- **Block-level directed graph** replaces the coloured de Bruijn graph
- **ProcessPoolExecutor** parallelism (replaces GNU Parallel)
- **Python 3.12+** (tested on 3.12, 3.13 and 3.14)
- PEP 8 snake_case package layout (`lace/`)
- Informative stderr summary: transcript/base counts, mean/median lengths,
  N50, redundancy reduction %, `--maxTran` cap warnings
- Pylint 10.00/10

### Benchmark — *Briza maxima* (NCBI taxid:29665, 1.17 M transcripts)

| Metric | v1.14.1 | v2.0.0 | Change |
|---|---:|---:|---:|
| Wall-clock time | 19 938 s | 94.5 s | **211× faster** |
| SuperTranscript count | 270 870 | 270 870 | identical |
| Total bases | 299 136 517 bp | 461 233 884 bp | **+54.2 %** |
| Mean SuperTranscript length | 1 104 bp | 1 703 bp | **+54.3 %** |
| Median SuperTranscript length | 634 bp | 643 bp | **+1.4 %** |

The increase in total bases reflects the block-level algorithm retaining
more biological sequence from multi-transcript clusters compared to the
de Bruijn graph approach in v1.14.1.

### Container

- `lace_2.0.sif` — Apptainer image (~281 MB)
- Base: Python 3.14 + minimap2 2.30 via micromamba
- Dependencies: networkx, tqdm, numpy, pandas

### Breaking changes

- Minimum Python version raised to 3.12
- BLAT is no longer used or required
- Module names changed to snake_case (`build_supertranscript`, `lace_run`)
