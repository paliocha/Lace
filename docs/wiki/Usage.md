# Usage

## Quick start

```bash
Lace transcripts.fasta clusters.txt --cores 16 -o output/
```

## CLI reference

```
usage: Lace [-h] [--cores CORES] [--maxTran MAXTRAN] [-o OUTPUTDIR]
            [-t] [-a] [-v]
            TranscriptsFile ClusterFile

Build SuperTranscripts from clustered transcript assemblies

positional arguments:
  TranscriptsFile       FASTA file containing all transcripts
  ClusterFile           Tab-delimited file mapping transcripts to clusters
                        (e.g. Corset output)

optional arguments:
  -h, --help            show this help message and exit
  --cores CORES         Number of parallel workers (default: 1)
  --maxTran MAXTRAN     Maximum transcripts per cluster (default: 50)
  -o, --outputDir DIR   Output directory (default: .)
  -t, --tidy            Remove intermediate files after running
  -a, --alternate       Create alternate annotations and metrics
  -v, --verbose         Enable debug logging
```

## Input files

### TranscriptsFile (FASTA)

A standard FASTA file with all transcripts.  Typically produced by a *de novo* assembler such as [Trinity](https://github.com/trinityrnaseq/trinityrnaseq):

```
>TRINITY_DN100_c0_g1_i1
ATGCGATCGATCG...
>TRINITY_DN100_c0_g1_i2
ATGCGATCGATCG...
```

### ClusterFile (tab-delimited)

Two columns: transcript ID and cluster name.  This is the standard output format of [Corset](https://github.com/Oshlack/Corset) ([Davidson & Oshlack 2014](https://doi.org/10.1186/s13059-014-0410-6)):

```
TRINITY_DN100_c0_g1_i1	Cluster-100.0
TRINITY_DN100_c0_g1_i2	Cluster-100.0
TRINITY_DN200_c0_g1_i1	Cluster-200.0
```

Clusters with a single transcript are passed through as-is (no alignment needed).

## Output files

| File | Description |
|---|---|
| `SuperDuper.fasta` | One SuperTranscript per cluster |
| `SuperDuper.gff` | Block annotation on the SuperTranscript coordinate system |

### FASTA header format

```
>Cluster-100.0 NoTrans:3,Whirls:0
ATGCGATCGATCG...
```

- **NoTrans** — number of transcripts used to build this SuperTranscript
- **Whirls** — number of graph cycles broken (0 = clean DAG)
- If both are `-1`, the cluster was too complex and the longest isoform was used as fallback

## Summary statistics

Lace 2.0 prints a before/after summary to stderr showing mean, median, and N50 lengths for both input transcripts and output SuperTranscripts:

```
  Transcripts in FASTA:              1,503,962
    Total bases:                     1,013,568,698 bp
    Mean transcript length:          674 bp
    Median transcript length:        401 bp
    N50:                             1,107 bp
  ...
  SuperTranscripts written:          463,540
    Total bases:                     654,732,311 bp
    Mean SuperTranscript length:     1,412 bp
    Median SuperTranscript length:   612 bp
    N50:                             2,404 bp
  Redundancy reduction:             69.2% (1,503,962 → 463,540 sequences)
```

## The `--maxTran` cap

Very large clusters (hundreds of transcripts) can be slow or produce overly complex graphs.  The `--maxTran` flag limits how many transcripts are used per cluster (default: 50).  Clusters exceeding this limit are logged with a warning:

```
WARNING: 329 cluster(s) hit the --maxTran=50 cap (only first 50 transcripts kept)
```

Increase `--maxTran` if you want to include more transcripts per cluster, or decrease it to speed up processing of noisy assemblies.

## Example

A 1% downsampled *Briza maxima* example is included in the repository:

```bash
git clone https://github.com/paliocha/Lace.git
cd Lace
Lace example/briza_maxima_1pct.fasta example/briza_maxima_1pct_clusters.txt \
    --cores 4 -o example_output/
```

This processes 4,583 transcripts across 2,708 clusters in a few seconds.
