# Installation

## Option 1 — Apptainer/Singularity container (recommended)

The easiest way to run Lace 2.0 is via the pre-built container image.  Download the SIF from the [GitHub Releases](https://github.com/paliocha/Lace/releases/tag/v2.0.0) page:

```bash
# Download (~281 MB)
wget https://github.com/paliocha/Lace/releases/download/v2.0.0/lace_2.0.sif

# Run
apptainer exec lace_2.0.sif Lace transcripts.fasta clusters.txt --cores 16 -o output/
```

The container bundles:
- Python 3.14 (via micromamba)
- minimap2 2.30
- networkx, tqdm, numpy, pandas

> **Note:** minimap2 is included inside the container — no separate installation needed.

## Option 2 — Install from source

### Prerequisites

- Python ≥ 3.12
- [minimap2](https://github.com/lh3/minimap2) on `$PATH`
- pip

### Steps

```bash
git clone https://github.com/paliocha/Lace.git
cd Lace
pip install -e .
```

### Python dependencies

```
networkx
tqdm
numpy
pandas
```

These are installed automatically by `pip install -e .`

### Verify

```bash
Lace --help
minimap2 --version
```

## Nextflow integration (nf-denovoslim)

When using Lace within the [nf-denovoslim](https://github.com/paliocha/nf-denovoslim) pipeline, the container is configured in `nextflow.config`:

```groovy
params {
    lace_container = "${projectDir}/containers/lace/lace_2.0.sif"
}
```

No manual installation is required — Nextflow manages the container automatically.
