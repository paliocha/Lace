# Algorithm

## Overview

Lace 2.0 constructs SuperTranscripts in three stages:

1. **Pairwise alignment** of all transcripts within a cluster
2. **Block-level directed graph** construction from alignment blocks
3. **Topological sort** to linearise the graph into a single SuperTranscript

```
Transcripts ──► minimap2 ──► Alignment blocks ──► Block graph ──► DAG ──► SuperTranscript
  (FASTA)        (PAF)        (AlignBlock)        (DiGraph)       ↓         (FASTA + GFF)
                                                                whirl
                                                              detection
```

## Stage 1 — minimap2 alignment (replaces BLAT)

Lace v1.14.1 used [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html) for pairwise transcript alignment.  Lace 2.0 replaces BLAT with [minimap2](https://github.com/lh3/minimap2) ([Li 2018](https://doi.org/10.1093/bioinformatics/bty191)), configured with flags that **mimic BLAT's alignment behaviour** for transcript-vs-transcript comparison:

```
minimap2 -c -X --eqx -k15 -w5 -N50 -p0.98 --no-long-join <fasta> <fasta>
```

| Flag | Purpose | BLAT equivalent |
|---|---|---|
| `-c` | Output base-level CIGAR in PAF | BLAT PSL with block coordinates |
| `-X` | All-vs-all mode (skip self-hits internally) | BLAT pairwise on same file |
| `--eqx` | Distinguish `=` (match) from `X` (mismatch) in CIGAR | Needed for accurate block extraction |
| `-k15` | 15-mer seed | Similar to BLAT's default 11-mer tileSize but tuned for transcripts |
| `-w5` | Window size 5 | Dense seeding comparable to BLAT's `stepSize=5` |
| `-N50` | Report up to 50 secondary alignments | Captures multi-exon overlaps like BLAT |
| `-p0.98` | Min 98% identity filter | Equivalent to BLAT's `-minIdentity=98` |
| `--no-long-join` | Disable long-gap joining | Prevents merging distinct exonic blocks |

This combination ensures minimap2 produces dense, high-identity pairwise alignments between transcripts — the same information BLAT provides — but runs **30–100× faster** and carries an MIT licence (no academic-only restriction).

### PAF parsing

Each alignment hit is decomposed into contiguous **alignment blocks** by walking the CIGAR string.  A block records: target name, query name, target start/end, query start/end, and strand.  Blocks shorter than 10 bp or with identity below the threshold are discarded.

## Stage 2 — Block-level directed graph

The block-level splice graph draws on the splice graph formalism introduced by [Heber *et al.* (2002)](https://doi.org/10.1093/bioinformatics/18.suppl_1.S181), adapted to operate on alignment blocks rather than individual bases.

Unlike v1.14.1's base-level coloured de Bruijn graph (one node per base, millions of nodes for large clusters), Lace 2.0 builds a **block-level** directed graph with typically 10–30 nodes per cluster:

1. Alignment blocks define shared regions between transcripts
2. Each unique genomic interval becomes a **node** containing the corresponding sequence
3. **Directed edges** preserve the ordering of blocks within each transcript — if block A precedes block B in any transcript, an edge A→B is created
4. Edges carry transcript membership for annotation

### Cycle handling ("whirls")

Transcripts with repeated or rearranged blocks can introduce cycles.  Lace detects these and breaks them by removing the minimum set of back-edges, converting the graph to a **DAG** (Directed Acyclic Graph).  Clusters with broken cycles are flagged with `Whirls:N` in the FASTA header (N = number of cycles removed).

## Stage 3 — Topological sort and sequence extraction

The DAG is topologically sorted using Kahn's algorithm.  The sorted node sequences are concatenated to produce the final SuperTranscript.  A GFF annotation records the original block coordinates on the SuperTranscript.

## Why block-level instead of base-level?

| Aspect | Base-level (v1.14.1) | Block-level (v2.0.0) |
|---|---|---|
| Nodes | ~millions per cluster | ~10–30 per cluster |
| Memory | High (1 node/base) | Low |
| Speed | Slow graph construction | Near-instant |
| Retained sequence | Lossy (k-mer collisions) | Lossless (full alignment blocks) |

The block-level approach retains **+54% more sequence** on the *Briza maxima* benchmark while being orders of magnitude faster.

## References

- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)
- Heber, S., Alekseyev, M., Sze, S.H., Tang, H. & Pevzner, P.A. (2002). Splicing graphs and EST assembly problem. *Bioinformatics*, 18(S1), S181–S188. [doi:10.1093/bioinformatics/18.suppl_1.S181](https://doi.org/10.1093/bioinformatics/18.suppl_1.S181)
- Davidson, N.M., Hawkins, A.D.K. & Oshlack, A. (2017). SuperTranscripts: a data driven reference for analysis and visualisation of transcriptomes. *Genome Biology*, 18, 148. [doi:10.1186/s13059-017-1284-1](https://doi.org/10.1186/s13059-017-1284-1)
