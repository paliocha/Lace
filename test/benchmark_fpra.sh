#!/bin/bash
#SBATCH --job-name=lace_FPRA_v2
#SBATCH --partition=orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=test/lace_FPRA_v2_%j.out
#SBATCH --error=test/lace_FPRA_v2_%j.err
set -euo pipefail

# ── Benchmark Lace 2.0 on FPRA data ─────────────────────────────────
# Mirrors the nf-denovoslim LACE step but runs standalone.
# Output goes to test/ (gitignored), does NOT overwrite pipeline files.

# ── Paths ────────────────────────────────────────────────────────────
SIF="/mnt/users/martpali/AnnualPerennial/Lace/lace_2.0.sif"
TRINITY="/mnt/project/FjellheimLab/martpali/AnnualPerennial/assemblies/FPRA-Trinity1/FPRA-Trinity.fasta"
CLUSTERS="/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/FPRA/clustering/corset-clusters.txt"
OUTDIR="$TMPDIR"
CORES=${SLURM_CPUS_PER_TASK:-16}

# Validate inputs exist
for f in "$SIF" "$TRINITY" "$CLUSTERS"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing $f" >&2
        exit 1
    fi
done

echo "================================================================"
echo "Lace 2.0 FPRA benchmark"
echo "================================================================"
echo "SIF:       $SIF"
echo "Trinity:   $TRINITY"
echo "Clusters:  $CLUSTERS"
echo "Output:    $OUTDIR"
echo "Cores:     $CORES"
echo "TMPDIR:    $TMPDIR"
echo "Node:      $(hostname)"
echo "Date:      $(date -Iseconds)"
echo "================================================================"

# ── Transcript / cluster counts ──────────────────────────────────────
echo ""
echo "Transcripts: $(grep -c '^>' "$TRINITY")"
echo "Cluster map lines: $(wc -l < "$CLUSTERS")"
echo ""

# ── Run Lace 2.0 ────────────────────────────────────────────────────
echo "Starting Lace 2.0 ..."
START=$(date +%s)

apptainer exec \
    --bind "$TMPDIR:$TMPDIR" \
    --bind /mnt/project:/mnt/project:ro \
    "$SIF" \
    Lace \
        "$TRINITY" \
        "$CLUSTERS" \
        -t \
        --cores "$CORES" \
        -o "$OUTDIR"

END=$(date +%s)
ELAPSED=$((END - START))
echo ""
echo "================================================================"
echo "Lace 2.0 completed in ${ELAPSED}s ($(date -ud @${ELAPSED} +%H:%M:%S))"
echo "================================================================"

# ── Copy results back to test/ (on NFS, for comparison) ─────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULT_DIR="${SCRIPT_DIR}/lace_fpra_v2_results"
mkdir -p "$RESULT_DIR"
cp "$OUTDIR/SuperDuper.fasta" "$RESULT_DIR/"
cp "$OUTDIR/SuperDuper.gff"   "$RESULT_DIR/" 2>/dev/null || true

# ── Basic stats ──────────────────────────────────────────────────────
echo ""
echo "Output stats:"
echo "  SuperDuper.fasta: $(wc -l < "$RESULT_DIR/SuperDuper.fasta") lines, $(du -h "$RESULT_DIR/SuperDuper.fasta" | cut -f1)"
echo "  SuperTranscripts: $(grep -c '^>' "$RESULT_DIR/SuperDuper.fasta")"
echo ""

# ── Compare with v1.14.1 if available ───────────────────────────────
V1_FASTA="/mnt/project/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/FPRA/work/ad/ce86e83dd1ac7e8cdfb4c4068cf594/lace_out/SuperDuper.fasta"
if [[ -f "$V1_FASTA" ]]; then
    echo "Comparison with v1.14.1:"
    V1_COUNT=$(grep -c '^>' "$V1_FASTA")
    V2_COUNT=$(grep -c '^>' "$RESULT_DIR/SuperDuper.fasta")
    echo "  v1.14.1 SuperTranscripts: $V1_COUNT"
    echo "  v2.0.0  SuperTranscripts: $V2_COUNT"
    echo "  Difference: $((V2_COUNT - V1_COUNT))"
else
    echo "v1.14.1 output not yet available for comparison."
fi

echo ""
echo "Done. Results in: $RESULT_DIR"
