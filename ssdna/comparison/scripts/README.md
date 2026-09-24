# dA40/dT40 comparison scripts

This directory contains reusable analysis code only. Generated figures, tables
and representative structures belong in [`../results/`](../results/README.md).

- `analyze_equal_window.py`: periodic-image-aware Python comparison of dA40
  and dT40. It generates the CSV datasets, representative PDB structures,
  PNG/PDF figures and validation metadata in `../results/`.

The Python script requires NumPy, pandas, SciPy and Matplotlib. Its input PSF
and DCD paths are explicit command-line arguments.

Simulation-specific NAMD inputs and exact job code snapshots are kept in each
sequence's own `scripts/` directory.
