# Final formal ssDNA simulations in 3 M KCl

This repository intentionally contains only the final formal poly(dA)40 and
poly(dT)40 simulations and the formal nanopore workflows derived from them.
Preliminary ssDNA runs, abandoned pore builds, failed job logs, duplicate
checkpoints, and exploratory analyses are excluded.

## Final ssDNA protocol

- CHARMM36 nucleic-acid force field with explicit water and 3 M KCl
- 293 K and 1 atm
- 50,000 minimization steps
- 0.5 ns heating
- 2 ns strong/intermediate restraint equilibration
- 2.5 ns weak-restraint equilibration
- 5 ns unrestrained equilibration
- 50 ns production at a 2 fs timestep

The production trajectory was run continuously as an initial 20 ns segment
(`prod20`) followed from its checkpoint by a 30 ns extension (`prod20to50`).
Both files therefore belong to the same final 50 ns trajectory; they are not
separate trial simulations.

Final 10 ns means:

| Metric | poly(dA)40 | poly(dT)40 |
|---|---:|---:|
| Radius of gyration (A) | 25.8415 | 25.9275 |
| End-to-end distance (A) | 52.7580 | 41.9074 |
| Base SASA (A2) | 101.9342 | 114.5460 |
| Stacked adjacent pairs | 17.1500 | 2.1000 |
| Unstacked-base fraction | 0.3488 | 0.9200 |
| Syn-like fraction | 0.6300 | 0.3062 |

## Repository contents

- `final-ssdna/`: final inputs, stage logs, non-duplicate checkpoints,
  50 ns analyses, and last-10-ns representative structures
- `formal-nanopore/`: only the current four formal 7AHL/3B07 workflows,
  clean pore templates, and folded-DNA starting placements
- `data-manifest/`: sizes and SHA-256 checksums for the omitted formal DCD files
- `PORE_JOB_STATUS.md`: current formal nanopore Slurm job chains

The 14 formal DCD files total approximately 10.92 GiB and are not stored in
ordinary Git. Their exact filenames, sizes, and hashes are retained in the
manifest so archived copies can be verified.

## Software

- NAMD 3.0.2
- VMD 1.9.2
- CHARMM36/CHARMM36m
- Python 3 with NumPy and SciPy for pore-placement validation
- Slurm

These data contain one production trajectory per sequence. Statistical claims
requiring uncertainty estimates should use independent replicas.
