# dA40/dT40 comparison scripts

This directory contains reusable analysis code only. Generated figures, tables
and representative structures belong in [`../results/`](../results/README.md).

- `analyze_equal_window.py`: periodic-image-aware, matched-age comparison of
  the first 5 ns of dA40 and dT40 production. It generates the CSV datasets,
  representative PDB structures, PNG/PDF figures and validation metadata in
  `../results/`. The starting frame can be changed with `--da-start-frame` and
  `--dt-start-frame` when another matched window is needed.

The Python script requires NumPy, pandas, SciPy and Matplotlib. Its input PSF
and DCD paths are explicit command-line arguments.

Simulation-specific NAMD inputs and exact job code snapshots are kept in each
sequence's own `scripts/` directory.
