"""Lace: SuperTranscript construction from clustered transcript assemblies.

Version 2.0.0 — rewrite with:
- minimap2 replacing BLAT (30-100× faster alignment, MIT licence)
- Block-level splice graph instead of base-level (100-1000× fewer nodes)
- Single-pass DFS cycle breaking instead of iterative nx.simple_cycles
- In-memory cluster dispatch (eliminates ~376K NFS file operations)
- Python 3.12+, modern typing, structured logging
"""
# pylint: disable=invalid-name  # package name "Lace" follows upstream convention

__version__ = "2.0.0"
