# Example — BMAX 1 % subsample

Quick-start example using 1 % of *Bromus maximus* Trinity transcripts
(2 708 clusters, 4 583 transcripts).

## Run

```bash
Lace bmax_1pct.fasta bmax_1pct_clusters.txt --cores 4 -o out
```

The full BMAX dataset (270 870 clusters, 1.17 M transcripts) completes
in ~95 seconds on 16 cores.

## Regenerate

The downsampling script requires access to the full BMAX Trinity assembly
and Corset clusters on the Orion NFS:

```bash
python3 downsample_bmax.py   # edit FRACTION for other sizes
```
