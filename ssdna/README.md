# ssDNA simulations

Each `<sequence>-3m-kcl/50ns/` directory contains the system, NAMD inputs,
stage logs, the latest checkpoint from each stage, analysis files and selected
structures from the 40–50 ns interval.

Production is a continuous 20 ns segment (`prod20`) plus a 30 ns continuation
(`prod20to50`), both at 2 fs. The preparation sequence is 50,000 minimization
steps, 0.5 ns heating, 2 ns restrained equilibration, 2.5 ns weak-restraint
equilibration and 5 ns unrestrained equilibration.

DCD filenames, sizes and hashes are in [the manifest](../data-manifest/ssdna_dcd.sha256).
Duplicate checkpoints and files ending in `.old` are not included.

## Analysis status

The archived CSVs, summaries, analysis scripts and representative PDB exports
need a periodic-image audit. They must not be treated as validated coiling or
base-exposure comparisons. No numerical results were recalculated in this
documentation update.

For nanopore construction, use the corrected `*_donor_intact.pdb` files in
[nanopore/common](../nanopore/common), not these older representative exports.
The [nanopore protocol](../nanopore/PROTOCOL.md) records their source frames.

Log normalization changes only the working-directory header to a repository
path. Energy, timing and other numerical records are unchanged.
