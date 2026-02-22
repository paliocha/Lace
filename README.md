# Lace 2.0

Build **SuperTranscripts** — one consensus sequence per gene — from a
clustered *de novo* transcriptome assembly.

Lace takes a FASTA of transcripts and a cluster file (e.g. from
[Corset](https://github.com/paliocha/Corset), MMseqs2, or CD-HIT) and
produces a single SuperTranscript for each gene by aligning within-cluster
transcripts with **minimap2** and threading them through a **block-level
directed graph**.

## Quick start

```bash
# Container (recommended)
apptainer exec lace_2.0.sif Lace transcripts.fasta clusters.txt -o out/

# Source install
pip install .
Lace transcripts.fasta clusters.txt -o out/
```

## What's new in v2.0

| Feature | v1.14.1 | v2.0.0 |
|---|---|---|
| Aligner | BLAT (academic licence) | minimap2 (MIT) |
| Graph | Coloured de Bruijn | Block-level directed |
| Parallelism | GNU Parallel | ProcessPoolExecutor |
| Python | 2.7 / 3.x | 3.12+ |
| Performance (*Briza maxima*, 1.17 M tx) | 5 h 32 min | 94.5 s (**211× faster**) |

See [CHANGELOG.md](CHANGELOG.md) for full details.

## Installation

### Apptainer / Singularity container

Download from the [v2.0.0 release](https://github.com/paliocha/Lace/releases/tag/v2.0.0):

```bash
wget https://github.com/paliocha/Lace/releases/download/v2.0.0/lace_2.0.sif
apptainer exec lace_2.0.sif Lace --help
```

### From source

Requires Python >= 3.12 and minimap2 on `$PATH`.

```bash
pip install .
# or with conda:
conda env create -f environment.yml && conda activate lace
```

### Nextflow (nf-denovoslim)

```groovy
params.lace_container = "${projectDir}/containers/lace/lace_2.0.sif"
```

## Usage

```
Lace <FastaFile> <ClusterFile> [-o outdir] [--maxTran N] [-t] [--cores N]
```

| Argument | Description |
|---|---|
| `FastaFile` | FASTA of input transcripts |
| `ClusterFile` | Tab-separated cluster file (e.g. Corset, mmseqs2, CD-HIT) |
| `-o` | Output directory (default: `Lace_output`) |
| `--maxTran` | Cap transcripts per cluster (default: 50) |
| `-t` | Write individual cluster FASTAs |
| `--cores` | Worker processes (default: all CPUs) |

### Output

```
out/
├── SuperDuper.fasta    # SuperTranscript sequences
└── SuperDuper.gff      # Transcript-to-SuperTranscript annotation
```

### Example

```bash
cd example/
Lace briza_maxima_1pct.fasta briza_maxima_1pct_clusters.txt -o bmax_test
```

## Algorithm

1. **Align** — minimap2 with BLAT-equivalent settings (`-p0.98 --eqx -k15 -w5`)
2. **Build graph** — block-level directed graph from alignment coordinates
3. **Linearise** — topological sort, break cycles, emit consensus

See the [wiki](https://github.com/paliocha/Lace/wiki/Algorithm) for details.

## References

- Davidson et al. (2017) SuperTranscripts: a data driven reference for analysis and visualisation of transcriptomes. *Genome Biology* 18:148.
  [doi:10.1186/s13059-017-1284-1](https://doi.org/10.1186/s13059-017-1284-1)
- Davidson & Oshlack (2014) Corset: enabling differential gene expression analysis for *de novo* assembled transcriptomes. *Genome Biology* 15:410.
  [doi:10.1186/s13059-014-0410-6](https://doi.org/10.1186/s13059-014-0410-6)
- Li (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34(18):3094–3100.
  [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)
- Heber et al. (2002) Splicing graphs and EST assembly problem. *Bioinformatics* 18(S1):S181–S188.
  [doi:10.1093/bioinformatics/18.suppl_1.S181](https://doi.org/10.1093/bioinformatics/18.suppl_1.S181)
- Grabherr et al. (2011) Trinity: reconstructing a full-length transcriptome without a genome from RNA-Seq data. *Nature Biotechnology* 29(7):644–652.
  [doi:10.1038/nbt.1883](https://doi.org/10.1038/nbt.1883)

## Licence

GPL-3.0 — see [LICENSE](LICENSE).

## Credits

Originally developed by Nadia Davidson and colleagues at the
[Oshlack Lab](https://github.com/Oshlack/Lace).
v2.0 rewrite by Martin Paliocha ([NMBU](https://www.nmbu.no/)).
