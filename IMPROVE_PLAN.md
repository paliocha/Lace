# Lace Improvement Plan

Style-guide optimisation plan for **Lace 1.14.1** following the same
tiered principles applied to Corset → corset-omp (hygiene → algorithmic →
low-level).  The goal is an in-tree fork (`lace-fast`) that can be
swapped into the nf-denovoslim pipeline by changing one container path.

---

## Baseline

| Species | Transcripts | Corset clusters | Multi-transcript clusters | Lace wall-time | Cores |
|---------|-------------|-----------------|---------------------------|----------------|-------|
| BMAX    | 460 386     | 270 870         | 94 437                    | 5 h 32 min     | 16    |
| FPRA    | 793 420     | 244 313         | ~120 000 (est.)           | > 6 h 47 min   | 16    |

Container: micromamba, Python 3.10.19, NetworkX 3.4.2, pandas 2.3.3,
numpy 1.26.4, BLAT v35, matplotlib-base 3.5.x.

Lace is **120× slower than Corset** for the same dataset despite being
embarrassingly parallel.  The single-cluster hot path
(`BuildSuperTranscript.SuperTran`) dominates.

---

## Inventory — Current Architecture

```
Lace_run.py          229 loc   Entry point, multiprocessing.Pool
BuildSuperTranscript 574 loc   Core algorithm (BLAT + graph + topo sort)
Mobius.py            120 loc   Alternative linear weighting (unused in pipeline)
Mobius_as.py         173 loc   Alternative alt-splicing view (unused)
Checker.py           261 loc   Alternate annotation metrics (unused, --alternate)
STViewer.py          151 loc   Matplotlib viewer (unused)
```

### Hot-Path Walkthrough (`BuildGraph`)

For **each** multi-transcript cluster:

1. **BLAT alignment** — `os.system("blat %s %s -minIdentity=98 %s.psl")`
   spawns one subprocess on disk.
2. **PSL parse** — `pd.read_table()` + row-wise `filt_dir()` (`.iloc[i,j]`
   in a Python for-loop).
3. **Double-strand handling** — if mixed strands detected, reverse-complement
   the negative transcripts, rewrite FASTA to disk, re-BLAT, re-parse.
4. **One NetworkX node per BASE** — every base in every transcript becomes a
   `G.add_node(i, Base='A', T1=pos, T2=None, ...)`.  A cluster of 5
   transcripts averaging 2 000 bp creates **10 000 nodes** with O(n) edges.
5. **Base-by-base merging** — for every aligned base pair in every BLAT
   block: redirect edges, copy attributes, `G.remove_node()`.  All Python
   for-loops.
6. **Chain simplification** — `successor_check()` walks single-degree chains,
   `merge_nodes()` concatenates `Base` strings, redirects edges.
7. **Whirl removal** — `nx.simple_cycles(C)` (Johnson's algorithm, worst-case
   exponential).  **Recalculated after every single whirl fix** inside a
   while-loop.  Bail-out at `max_edges=100`.
8. **Topological sort** — `nx.topological_sort(C)`, concatenate `Base` attrs.

### Parallelism Model

`multiprocessing.Pool(processes=ncore)` with **`chunksize=1`**.  Each worker
receives a 2-element list `[filename, "X of Y"]`.  The per-cluster files
are written to disk before the pool starts (one `.fasta` per gene in
`outdir/`).

### I/O Pattern

- Write N `.fasta` files (one per cluster) before pool launch.
- Each worker reads its `.fasta`, writes `.psl`, reads `.psl`.
- Double-strand clusters write a second `_stranded.fasta` and `.psl`.
- Cleanup via `os.system("rm ...")` or `os.system("mv ...")` per cluster.

**Total file operations**: ≈ 4 files × 94 000 clusters = **376 000 files**
created + deleted on a shared NFS filesystem.

---

## Tier 0 — Hygiene (zero behavioural change)

| # | Issue | Fix | Risk |
|---|-------|-----|------|
| H1 | `os.system()` for BLAT, rm, mv, mkdir | `subprocess.run(..., check=True)` | None |
| H2 | `os.system('rm/mv ...')` for cleanup | `pathlib.Path.unlink()` / `.rename()` | None |
| H3 | Bare `except:` in `worker()` | `except Exception as exc:` + structured logging | None |
| H4 | `filt_dir()` row-wise `.iloc[i,j]` loop | Vectorised pandas / numpy column ops | None |
| H5 | `Reverse_complement` char-by-char loop | `str.maketrans` + `[::-1]` (or Biopython `reverse_complement()`) | None |
| H6 | Mutable default `transcripts={}` pattern | Explicit `None` default, assign inside function | None |
| H7 | String formatting `%s` → f-strings | Modern Python style | None |
| H8 | `#Lacing together different transcripts` comment style | PEP 257 docstrings | None |
| H9 | `conmerge=True` dead-code branch | Remove always-true guard | None |
| H10 | Unused imports (`matplotlib.pyplot`, `cm`) everywhere | Remove | None |
| H11 | `sys.setrecursionlimit(100000)` commented out | Delete | None |
| H12 | `from multiprocessing import Process` (unused) | Delete | None |
| H13 | Type hints (Python 3.10 `X \| None` syntax) | Add to all public signatures | None |
| H14 | `chunksize=1` in `Pool.map_async` | `chunksize = max(1, len(fnames) // (ncore * 4))` | None |

### Expected impact

Marginal runtime improvement (H4, H5 speed up inner loops slightly), but
cleans the codebase for safe algorithmic work.

---

## Tier 1 — I/O Elimination

The single biggest source of wall-clock waste is **hundreds of thousands of
small-file operations on NFS**.

| # | Issue | Fix | Impact |
|---|-------|-----|--------|
| IO1 | One `.fasta` per cluster on disk | Pass sequences in-memory via Pool `starmap` — worker receives `(gene_id, {tran_id: seq})` dict, never touches disk | Eliminate ~376 K file ops |
| IO2 | BLAT writes `.psl` to disk, parsed back by pandas | Pipe BLAT stdout (`-out=psl -noHead stdout`) or capture via `subprocess.PIPE` | Eliminate 2 files per cluster |
| IO3 | Double-strand: rewrite `_stranded.fasta` + re-BLAT | Reverse-complement in-memory, write tmp to `/dev/shm` or `$TMPDIR` (local SSD) if BLAT requires a file | Eliminate NFS round-trips |
| IO4 | Cleanup shells `rm`, `mv` per cluster | Already eliminated by IO1 | — |

### Detailed Design — IO1

```python
def worker(args: tuple[str, dict[str, str], int, int]) -> tuple[str, str, int, int]:
    """Process one cluster entirely in-memory."""
    gene_id, transcripts, job_idx, total = args
    # Write temp FASTA to local SSD ($TMPDIR), not NFS
    tmp = Path(os.environ.get("TMPDIR", "/tmp")) / f"{gene_id}.fasta"
    ...
```

With `multiprocessing.Pool.starmap()` the dict is pickled once (IPC) —
cheaper than writing + reading 2 files per cluster on NFS.

### Expected impact

Lace BMAX currently creates/deletes ~376 K NFS files.  Node-local SSD
(`$TMPDIR`) IOPS is 50–100× faster than NFS; eliminating the NFS round-trips
should yield **1.5–2× wall-clock reduction** from I/O alone.

---

## Tier 2 — Algorithmic

### A1: Replace BLAT with minimap2

| | BLAT v35 | minimap2 2.28 |
|---|----------|---------------|
| Licence | Academic-only | MIT |
| Speed | ~30 s for 50 transcripts (typical cluster) | ~0.2 s (30-fold faster, PAF output) |
| Index | Rebuilds per cluster (all-vs-all) | Can index once, query many |
| Output | PSL (tab-separated, 5-line header) | PAF (tab-separated, no header, richer) |
| Accuracy | minIdentity=98 (nt) | `--eqx -c -X -k15 -w5 -N50 -p0.98` matches |

**Plan**: Replace `os.system("blat ...")` with
`subprocess.run(["minimap2", "-c", "-X", "--eqx", ...], capture_output=True)`.
Parse PAF instead of PSL (simpler: no 5-line header, direct column mapping).

**Impact**: 30–100× speedup of the per-cluster alignment step, which is
≈60% of hot-path time.

### A2: Block-level graph instead of base-level

The current design creates **one node per base** then merges aligned bases
one-by-one.  A cluster with 5 × 2 kb transcripts starts with 10 000 nodes.

**Replace with block-level graph**:

1. Parse alignment blocks (contiguous matches) from minimap2 PAF.
2. Build an interval union over all transcripts (merge overlapping blocks).
3. Each **block** (contiguous aligned region) becomes one node.
   Transcript-private regions become their own nodes.
4. Edges connect blocks in transcript order.
5. Chain simplification is trivial (most chains are already single-node).

A typical cluster with 5 transcripts and 3 splice variants produces
**~10–30 block nodes** instead of 10 000 base nodes.

**Impact**: 100–1 000× fewer graph nodes → proportional speedup in all
graph operations.

### A3: Single-pass cycle breaking

Current code:
```python
while len(whirls) > 0:
    whirls = list(nx.simple_cycles(C))   # recalculate ALL cycles
    # fix one whirl
    # repeat
```

This is **O(k × |cycles|)** where `|cycles|` can be exponential.

**Replace with**: DFS-based back-edge detection + single-pass breaking.
For a DAG with occasional tandem duplications:

```python
# Identify all back-edges in one DFS
back_edges = set()
for edge in nx.dfs_edges(C):
    if C.has_edge(edge[1], edge[0]):
        back_edges.add((edge[1], edge[0]))

# Break all back-edges by node duplication (same logic as current)
for u, v in back_edges:
    duplicate_and_redirect(C, u, v)
```

**Impact**: Worst-case exponential → O(V + E).  Eliminates the main
source of pathological slowdowns on repeat-rich clusters.

### A4: Vectorise base merging with numpy

Even with block-level graphs (A2), merging block coordinates benefits
from numpy vectorisation:

```python
# Instead of Python for-loop over block_seq[i]
tpos = np.arange(tStart[i], tStart[i] + block_size)
qpos = np.arange(qStart[i], qStart[i] + block_size)
# Batch update node_dict arrays
node_dict_arr[qName_idx, qpos] = node_dict_arr[tName_idx, tpos]
```

**Impact**: 10–50× speedup of the merge step for large clusters.

---

## Tier 3 — Architecture & Low-Level

### L1: Replace NetworkX with lightweight custom DAG

NetworkX stores each node as a dict-of-dicts — massive per-node overhead.
For a DAG of ~10–30 block nodes per cluster, a simple adjacency list
(plain Python `dict[int, list[int]]` + separate `node_data: list[str]`)
is sufficient and avoids:

- Per-node `__hash__` / `__eq__` overhead
- GC pressure from millions of small dicts
- Import cost (~200 ms for `import networkx`)

**Alternative**: Keep NetworkX but operate exclusively on the block-level
graph (A2), which is small enough that overhead is negligible.

**Recommendation**: Implement A2 first.  If profiling still shows
NetworkX overhead, replace with custom DAG.

### L2: Cython / C extension for hot inner loop

If the block-merging inner loop (Tier 2, A2–A4) is still the bottleneck
after numpy vectorisation, move the merge kernel to Cython or a small C
extension compiled at container build time.  This follows the same
pattern as Corset's SSE2 intersection — a surgical replacement of the
innermost loop.

### L3: concurrent.futures.ProcessPoolExecutor

Replace `multiprocessing.Pool` with `concurrent.futures.ProcessPoolExecutor`
for:
- Cleaner exception propagation (futures vs `AsyncResult`)
- `as_completed()` progress reporting
- Better integration with structured logging

### L4: Work-stealing load balance

The current `chunksize=1` + `map_async` dispatches clusters in input order.
Large clusters (50 transcripts, complex repeats) can create stragglers.

**Fix**: Sort clusters by transcript count (descending) before dispatch.
Largest clusters start first; smaller ones fill gaps (natural work-stealing
with sorted longest-job-first scheduling).

### L5: Progress reporting

Replace `print("Processing cluster %s - %s")` with proper logging:

```python
import logging
from tqdm import tqdm

log = logging.getLogger("lace")

# In pool dispatch:
with tqdm(total=len(fnames), desc="Building SuperTranscripts") as pbar:
    for result in pool.imap_unordered(worker, fnames, chunksize=chunk):
        pbar.update(1)
        results.append(result)
```

---

## Tier 4 — Container & Environment

| # | Issue | Fix | Impact |
|---|-------|-----|--------|
| C1 | Python 3.10 | Upgrade to 3.12+ (faster interpreter, per-PEP 709 inlined comprehensions) | 5–15% baseline speedup |
| C2 | BLAT v35 (academic licence) | Replace with minimap2 (MIT) — see A1 | Licence + speed |
| C3 | matplotlib imported but unused in pipeline | Remove from container (saves ~100 MB, eliminates `MPLCONFIGDIR` workaround) | Smaller image |
| C4 | Lace installed via conda (immutable) | Fork to `lace-fast`, install from local source in container | Full control |
| C5 | Container def uses micromamba base | Keep micromamba for minimap2 + numpy; add `tqdm` | Minimal change |

### Container Blueprint

```singularity
Bootstrap: docker
From: docker.io/mambaorg/micromamba:1.5.8

%files
    lace-fast/ /opt/lace-fast

%post
    micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 \
        minimap2=2.28 \
        numpy>=1.26 \
        pandas>=2.2 \
        networkx>=3.4 \
        tqdm \
        && micromamba clean --all --yes

    pip install --no-deps /opt/lace-fast

%environment
    export PATH="/opt/conda/bin:${PATH}"
```

---

## Implementation Roadmap

### Phase 1 — Fork & Hygiene (1 day)

1. Copy Lace source from container to `lace-fast/` in nf-denovoslim tree.
2. Apply all Tier 0 fixes (H1–H14).
3. Validate: run on 10 random BMAX clusters, diff output against v1.14.1.
4. Git commit: `style: initial lace-fast fork, hygiene pass`.

### Phase 2 — I/O Elimination (1 day)

1. Implement IO1–IO4 (in-memory sequence passing, TMPDIR for alignment files).
2. Validate: full BMAX run, diff SuperDuper.fasta.
3. Benchmark: expect 1.5–2× speedup.
4. Git commit: `perf: eliminate NFS I/O, in-memory cluster dispatch`.

### Phase 3 — minimap2 + Block Graph (2–3 days)

1. Replace BLAT with minimap2 (A1).
2. Parse PAF, build block-level graph (A2).
3. Single-pass cycle breaking (A3).
4. Validate: full BMAX run, compare SuperTranscripts (sequence diff,
   length correlation, BUSCO score).
5. Benchmark: expect 10–50× total speedup vs baseline.
6. Git commit: `perf: minimap2 aligner, block-level splice graph`.

### Phase 4 — Polish & Container (1 day)

1. Python 3.12 container build (C1–C5).
2. Load-balancing sort (L4), progress bar (L5).
3. Final validation: all three species, compare cluster-by-cluster.
4. Git commit: `feat: lace-fast 2.0 container, Python 3.12`.

### Phase 5 — Pipeline integration

1. Build `containers/lace/lace_fast_2.0.sif`.
2. Update `nextflow.config`: container path.
3. Update `modules/lace.nf` if CLI changes.
4. Update CLAUDE.md, CHANGELOG.md, README.md.

---

## Validation Protocol

Following the same practice as Corset v1.09 → v1.10:

1. **Sequence-level diff**: Compare `SuperDuper.fasta` header-by-header.
   Expect identical or near-identical sequences (minimap2 block boundaries
   may differ by a few bases at alignment edges).
2. **Length correlation**: Plot lace-fast vs baseline SuperTranscript
   lengths.  Expect Pearson r > 0.99.
3. **Downstream invariance**: Run the full pipeline (TD2 → SELECT_BEST_ORF
   → BUSCO).  BUSCO score must not decrease.
4. **Whirl statistics**: Compare whirl counts per cluster.  The new
   cycle-breaking algorithm may report different counts but should produce
   valid DAGs.

---

## Expected Final Performance

| Phase | BMAX est. | Speedup vs baseline |
|-------|-----------|---------------------|
| Baseline (v1.14.1) | 5 h 32 min | 1× |
| +Hygiene (Tier 0) | ~5 h 15 min | ~1.05× |
| +I/O elimination (Tier 1) | ~3 h | ~1.8× |
| +minimap2 + block graph (Tier 2) | **10–20 min** | **17–33×** |
| +Python 3.12 + polish (Tier 3–4) | **8–15 min** | **22–42×** |

The alignment step (currently ~60% of time) speeds up 30–100× with
minimap2, and the graph step (currently ~35% of time) speeds up 100–1000×
with block-level nodes.  The remaining 5% (I/O, pool overhead) is
addressed by Tiers 0–1.

Target: **Lace BMAX < 15 minutes on 16 cores** — comparable to Corset.

---

## Risks & Mitigations

| Risk | Likelihood | Mitigation |
|------|-----------|------------|
| minimap2 alignment boundaries differ from BLAT | Moderate | Validate with downstream BUSCO; tune `-k` and `-w` |
| Block-level graph misses intra-block variations | Low | Fall back to base-level for blocks with mismatches |
| `nx.simple_cycles` removal changes SuperTranscript topology | Low | Compare whirl stats; DFS cycle-breaking is strictly more conservative |
| Pickle overhead for large clusters in Pool | Low | Use `multiprocessing.shared_memory` for sequences > 1 MB |
| Python 3.12 breaks NetworkX | Very low | NetworkX 3.4 already supports 3.12 |

---

## Files to Create / Modify

| File | Action |
|------|--------|
| `lace-fast/` (new tree) | Fork of Lace source |
| `lace-fast/lace_fast/build_supertranscript.py` | Rewritten core (A1–A4) |
| `lace-fast/lace_fast/run.py` | Cleaned entry point (H1–H14, IO1, L3–L5) |
| `lace-fast/pyproject.toml` | Package metadata |
| `containers/lace/Lace_fast.def` | New container def (C1–C5) |
| `modules/lace.nf` | Update container + any CLI changes |
| `nextflow.config` | Container path |
| `CHANGELOG.md` | Document Lace upgrade |
| `CLAUDE.md` | Update Lace section |
---

## Implementation Log

### 2026-02-20 — Initial implementation (Tiers 0–4)

All tiers implemented in a single pass into `lace-fast/`:

| File | Tiers | Notes |
|------|-------|-------|
| `lace-fast/lace_fast/build_supertranscript.py` | H1–H14, IO1–IO3, A1–A3 | Full rewrite: minimap2 PAF parsing, block-level splice graph with union-find merging, DFS cycle breaking via `nx.find_cycle`, `str.maketrans` reverse complement, structured logging |
| `lace-fast/lace_fast/run.py` | H1–H14, IO1, L3–L5 | `ProcessPoolExecutor` + `as_completed`, largest-first scheduling, `tqdm` progress, f-string banner, `argparse` clean-up |
| `lace-fast/pyproject.toml` | C4 | PEP 621 metadata, `>=3.12` requirement, no matplotlib |
| `lace-fast/Lace_fast.def` | C1–C5 | Python 3.13, minimap2 2.28, pip install from source |
| `lace-fast/build.sh` | — | Container build helper |

#### BLAT vs minimap2 vs miniprot2

The plan specifies **minimap2**, not miniprot2. Key distinction:

- **minimap2**: nucleotide ↔ nucleotide aligner, `-X` all-vs-all mode.
  Direct replacement for BLAT in Lace's transcript alignment use case.
- **miniprot2**: protein → genome aligner. **Wrong tool** — Lace operates
  on nucleotide transcript sequences and needs nt-level block coordinates.

minimap2 flags used: `-c -X --eqx -k15 -w5 -N50 -p0.98`

#### Architecture decisions

1. **Block-level graph**: Uses alignment breakpoints to partition each
   transcript into intervals. Blocks across transcripts are merged via
   union-find when they span identical-length aligned regions. Produces
   ~10–30 nodes per cluster vs ~10 000 in the original.

2. **Cycle breaking**: Replaced `while len(whirls) > 0: nx.simple_cycles()`
   loop with iterative `nx.find_cycle()` + node duplication. Each iteration
   finds and breaks exactly one cycle. Much faster than enumerating all
   cycles (Johnson's algorithm is worst-case exponential).

3. **In-memory dispatch**: Sequences passed as `dict[str, str]` via pickle
   IPC. Workers write temp FASTA to `$TMPDIR` (node-local SSD) only for
   minimap2 invocation, then delete immediately.

4. **Python 3.13**: Container targets 3.13 for PEP 709 inlined comprehensions,
   specializing adaptive interpreter, and general 5–15% speedup over 3.10.

#### Remaining work

- [ ] Build container and validate on BMAX 10-cluster subset
- [ ] Full BMAX validation run, diff SuperDuper.fasta
- [ ] Update `nf-denovoslim/modules/lace.nf` container path
- [ ] Benchmark and record timings