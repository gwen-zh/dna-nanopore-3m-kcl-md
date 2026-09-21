# ssDNA and nanopore simulations in 3 M KCl

poly(dA)40 / poly(dT)40 DNA-only simulations and four DNA–7AHL/3B07
nanopore workflows. Earlier exploratory simulations are not included.

## Validated ssDNA comparison — 2026-09-21

[Figures, CSV data, representative structures and reproducible analysis](ssdna/analysis-2026-09-21/README.md)
compare equal 5 ns windows after PSF-bond periodic reconstruction. In these
windows, dA40 is more compact and less rod-like than dT40, with higher syn and
outward base orientation. dT40 has substantially more loss of adjacent
stacking. The comparison uses dA 45.01–50.00 ns and the complete 0.01–5.00 ns
production from the [repaired dT branch](ssdna/polydt40-3m-kcl/repaired-5ns/README.md).

The repaired dT run completed normally and remained at least 85.886 Å from its
nearest periodic image. The older dT 50 ns trajectory is retained only as an
invalidated archive and is not used for current structural conclusions.

## dA trajectory results and movies — 2026-09-17

[Results, MP4 movies and downloadable inputs](nanopore/results/2026-09-17/README.md):
dA–7AHL completed its 5 ns production segment; dA–3B07 is published as an
explicitly labelled 0–4.825 ns in-progress preview. Neither exported trajectory
shows complete DNA passage. These are real-coordinate movies of pore-mouth
interactions and conformational changes, not illustrations of assumed passage.

The two old dT pore runs were stopped because their donor was affected by the
DNA-only periodic-image problem. The repaired dT DNA-only production is now
complete and validated; any new dT pore run must use that corrected branch.
The current dA trajectories use NAMD 2.14 CUDA; see the
[run-specific protocol](nanopore/results/2026-09-17/PROTOCOL.md).

## Protocols

- [Nanopore protocol / 过孔模拟流程](nanopore/PROTOCOL.md): system preparation,
  force fields, restraints, stage durations, voltage, Slurm execution and checks.
- [ssDNA dataset](ssdna/README.md): current validated sources and archived runs.
- [Job status](PORE_JOB_STATUS.md): timestamped cluster snapshot.

The DNA-only preparation is 50,000 minimization steps → 0.5 ns heating →
2 ns restrained equilibration → 2.5 ns weak-restraint equilibration →
5 ns unrestrained equilibration at 2 fs, 293 K and 1 atm. dA then has 50 ns
production; repaired dT currently has one validated 5 ns production block.

The nanopore sequence is 50,000 minimization steps → 0.5 ns heating →
4 ns staged equilibration → 5 ns electric-field production at 1 fs.
The pore boxes are fixed; this is not the pressure-coupled DNA-only protocol.

## Repository contents

- `ssdna/polyda40-3m-kcl/50ns/`: dA inputs, stage logs and checkpoints.
- `ssdna/polydt40-3m-kcl/repaired-5ns/`: repaired dT inputs and validation.
- `ssdna/analysis-2026-09-21/`: validated figures, tables, PDBs and source code.
- `nanopore/<sequence>-3m-kcl-<pore>/run/`: DNA placement, system settings and
  ion-placement checks; `nanopore/common/` contains the shared workflow.
- `nanopore/7ahl-template/`, `nanopore/3b07-template/`: processed pore/membrane/water templates.
- `data-manifest/`: filenames, sizes and SHA-256 checksums for omitted DCD and
  large solvated-system files.

No trajectory or checkpoint was deleted from the running cluster project.
Repository naming is independent of that project's working-directory names.
In archived ssDNA logs, only the working-directory header has been normalized
to the repository path; numerical simulation records are unchanged.

## Data-quality note

A periodic-image problem was confirmed in the older dT 50 ns trajectory. Its
files remain available for provenance but must not be used to establish dA/dT
differences. The repaired dT branch was rebuilt in a 180 Å water box and passed
all-frame bond and image-distance checks. The 2026-09-21 comparison explicitly
reconstructs both DNAs through PSF bonds before analysis.

Each sequence currently has one usable trajectory, not an independent-replica
uncertainty estimate. The equal windows are also at different trajectory ages:
dA uses its last 5 ns, whereas repaired dT uses its first 5 ns production block.

## Software

NAMD 2.14 CUDA for the current dA pore trajectories (NAMD 3.0.2 in earlier
workflow stages); VMD 1.9.2; CHARMM nucleic-acid, protein, lipid and water/ion
parameters as specified in the protocol; Python 3 with NumPy/SciPy; Slurm.
Force-field installation paths and the NAMD module name are cluster-specific.
