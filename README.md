# ssDNA and nanopore simulations in 3 M KCl

poly(dA)40 / poly(dT)40 DNA-only simulations and four DNA–7AHL/3B07
nanopore workflows. Earlier exploratory simulations are not included.

## dA trajectory results and movies — 2026-09-17

[Results, MP4 movies and downloadable inputs](nanopore/results/2026-09-17/README.md):
dA–7AHL completed its 5 ns production segment; dA–3B07 is published as an
explicitly labelled 0–4.825 ns in-progress preview. Neither exported trajectory
shows complete DNA passage. These are real-coordinate movies of pore-mouth
interactions and conformational changes, not illustrations of assumed passage.

The two old dT pore runs were stopped while the DNA-only periodic-image issue
is repaired. New dT pore jobs wait for successful sampling and validation.
The current dA trajectories use NAMD 2.14 CUDA; see the
[run-specific protocol](nanopore/results/2026-09-17/PROTOCOL.md).

## Protocols

- [Nanopore protocol / 过孔模拟流程](nanopore/PROTOCOL.md): system preparation,
  force fields, restraints, stage durations, voltage, Slurm execution and checks.
- [ssDNA dataset](ssdna/README.md): the continuous 50 ns DNA-only trajectories.
- [Job status](PORE_JOB_STATUS.md): timestamped cluster snapshot.

The ssDNA sequence is 50,000 minimization steps → 0.5 ns heating →
2 ns restrained equilibration → 2.5 ns weak-restraint equilibration →
5 ns unrestrained equilibration → 50 ns production at 2 fs, 293 K and 1 atm.
The production consists of a 20 ns segment followed by a 30 ns continuation.

The nanopore sequence is 50,000 minimization steps → 0.5 ns heating →
4 ns staged equilibration → 5 ns electric-field production at 1 fs.
The pore boxes are fixed; this is not the pressure-coupled DNA-only protocol.

## Repository contents

- `ssdna/<sequence>-3m-kcl/50ns/`: inputs, stage logs, checkpoints and analyses.
- `nanopore/<sequence>-3m-kcl-<pore>/run/`: DNA placement, system settings and
  ion-placement checks; `nanopore/common/` contains the shared workflow.
- `nanopore/7ahl-template/`, `nanopore/3b07-template/`: processed pore/membrane/water templates.
- `data-manifest/`: filenames, sizes and SHA-256 checksums for 14 omitted DCD
  files (approximately 10.92 GiB).

No trajectory or checkpoint was deleted from the running cluster project.
Repository naming is independent of that project's working-directory names.
In archived ssDNA logs, only the working-directory header has been normalized
to the repository path; numerical simulation records are unchanged.

## Data-quality note

A periodic-image problem was found in the earlier exported ssDNA structures.
The nanopore inputs now use intact DNA re-extracted from the original DCD,
joined through PSF bonds and validated. See the protocol for exact frame indices.

Existing ssDNA analysis CSVs, summaries and older representative exports are
retained as archived results, **pending an image-handling audit**. Do not use
them to establish differences in coiling or base exposure before that audit.
The original trajectories have not been modified. Each sequence currently has
one trajectory, not an independent-replica uncertainty estimate.

## Software

NAMD 2.14 CUDA for the current dA pore trajectories (NAMD 3.0.2 in earlier
workflow stages); VMD 1.9.2; CHARMM nucleic-acid, protein, lipid and water/ion
parameters as specified in the protocol; Python 3 with NumPy/SciPy; Slurm.
Force-field installation paths and the NAMD module name are cluster-specific.
