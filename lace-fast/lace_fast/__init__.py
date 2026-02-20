"""lace-fast: Fast SuperTranscript construction.

Rewrite of Lace 1.14.1 with:
- minimap2 replacing BLAT (30-100× faster alignment)
- Block-level splice graph instead of base-level (100-1000× fewer nodes)
- Single-pass DFS cycle breaking instead of iterative nx.simple_cycles
- In-memory cluster dispatch (eliminates ~376K NFS file operations)
- Python 3.13, modern typing, structured logging
"""

__version__ = "2.0.0"
