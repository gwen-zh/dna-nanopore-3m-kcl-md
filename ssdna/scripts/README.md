# ssDNA scripts

This directory contains reusable analysis code only. Generated figures, tables
and representative structures belong in [`../results/`](../results/README.md).

- `analyze_equal_window.py`: current periodic-image-aware Python comparison of
  dA40 and repaired dT40. It creates all 2026-09-21 CSV, PDB, PNG, PDF and
  validation outputs.
- `analyze_conformation.tcl`: older VMD conformation analysis used by the
  archived 50 ns run workflow.
- `analyze_base_flipping.tcl`: older VMD base-orientation analysis used by the
  archived run and nanopore analysis workflows.

The Python script requires NumPy, pandas, SciPy and Matplotlib. Its input PSF
and DCD paths are explicit command-line arguments; its default output is
`ssdna/results/2026-09-21/`.

Simulation-specific NAMD inputs and exact job code snapshots remain beside the
corresponding simulation dataset because they are provenance records rather
than reusable analysis scripts.
