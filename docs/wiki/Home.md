# Lace 2.0 — SuperTranscript Construction

<p align="center">
<img src="https://github.com/Quarkins/SuperTranscript/raw/master/WikiFigs/logo.png" height="200" />
</p>

**Lace** constructs **SuperTranscripts** ([Davidson *et al.* 2017](https://doi.org/10.1186/s13059-017-1284-1)) — single, linear representations of the transcriptome that contain all the sequence from every transcript in a cluster while preserving base ordering.  Given a set of transcripts and a clustering (e.g. [Corset](https://github.com/Oshlack/Corset) output; [Davidson & Oshlack 2014](https://doi.org/10.1186/s13059-014-0410-6)), Lace produces one SuperTranscript per cluster together with a GFF annotation of the constituent blocks.

## What's new in v2.0.0

Lace 2.0 is a ground-up rewrite of the [original Lace](https://github.com/Oshlack/Lace) (v1.14.1, Hawkins *et al.* 2017).

| Feature | v1.14.1 | v2.0.0 |
|---|---|---|
| Aligner | BLAT (academic licence) | **[minimap2](https://github.com/lh3/minimap2)** ([Li 2018](https://doi.org/10.1093/bioinformatics/bty191), MIT licence) |
| Graph | Base-level coloured de Bruijn graph | **Block-level directed splice graph** ([Heber *et al.* 2002](https://doi.org/10.1093/bioinformatics/18.suppl_1.S181), ~10–30 nodes) |
| Parallelism | GNU Parallel + per-cluster disk I/O | **ProcessPoolExecutor** — in-memory dispatch |
| Python | 2.7 / 3.5 | **3.12+** (tested 3.12, 3.13, 3.14) |
| Code style | Mixed | PEP 8 snake_case, pylint 10/10 |

### Benchmark — *Briza maxima* (1.17 M transcripts, 16 cores)

| Metric | v1.14.1 | v2.0.0 | Change |
|---|---:|---:|---:|
| Wall-clock time | 19 938 s (5.5 h) | 94.5 s | **211× faster** |
| SuperTranscript count | 270 870 | 270 870 | identical |
| Total bases | 299 M bp | 461 M bp | **+54.2 %** |
| Mean SuperTranscript length | 1 104 bp | 1 703 bp | **+54.3 %** |
| Median SuperTranscript length | 634 bp | 643 bp | **+1.4 %** |

The increase in total bases reflects the block-level algorithm retaining more biological sequence from multi-transcript clusters compared to the de Bruijn graph approach.

## Summary output

Lace 2.0 logs a comprehensive before/after summary to stderr:

```
  Transcripts in FASTA:              1,171,166
    Total bases:                     793,367,210 bp
    Mean transcript length:          677 bp
    Median transcript length:        429 bp
    N50:                             1,071 bp
  Transcript clusters (total):     270,870
    Multi-transcript clusters:       116,458
    Singleton clusters:              154,412

  SuperTranscripts written:          270,870
    Assembled (multi-transcript):    116,458
    Pass-through (singletons):       154,412
    Total bases:                     461,233,884 bp
    Mean SuperTranscript length:     1,703 bp
    Median SuperTranscript length:   643 bp
    N50:                             2,911 bp
  Redundancy reduction:             76.9% (1,171,166 → 270,870 sequences)

BUILT SUPERTRANSCRIPTS ---- 94.5 seconds ----
```

## Wiki contents

- [Algorithm](Algorithm) — how minimap2 + block-level graph works
- [Installation](Installation) — container and source install
- [Usage](Usage) — CLI reference and example

## References

- Davidson, N.M., Hawkins, A.D.K. & Oshlack, A. (2017). SuperTranscripts: a data driven reference for analysis and visualisation of transcriptomes. *Genome Biology*, 18, 148. [doi:10.1186/s13059-017-1284-1](https://doi.org/10.1186/s13059-017-1284-1)
- Davidson, N.M. & Oshlack, A. (2014). Corset: enabling differential gene expression analysis for *de novo* assembled transcriptomes. *Genome Biology*, 15, 410. [doi:10.1186/s13059-014-0410-6](https://doi.org/10.1186/s13059-014-0410-6)
- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)
- Heber, S., Alekseyev, M., Sze, S.H., Tang, H. & Pevzner, P.A. (2002). Splicing graphs and EST assembly problem. *Bioinformatics*, 18(S1), S181–S188. [doi:10.1093/bioinformatics/18.suppl_1.S181](https://doi.org/10.1093/bioinformatics/18.suppl_1.S181)
