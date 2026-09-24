# ssDNA and nanopore simulations in 3 M KCl

poly(dA)40 / poly(dT)40 DNA-only simulations and four DNA–7AHL/3B07
nanopore workflows, with simulation inputs, analysis data and figures.

## ssDNA structural comparison

[Figures, CSV data, representative structures and reproducible analysis](ssdna/comparison/results/README.md)
compare equal 5 ns windows after PSF-bond periodic reconstruction. In these
windows, dA40 is more compact and less rod-like than dT40, with higher syn and
outward base orientation. dT40 has substantially more loss of adjacent
stacking. The comparison uses dA 45.01–50.00 ns and the complete 0.01–5.00 ns
production from the [dT simulation](ssdna/polydt40-3m-kcl/README.md).


## dA nanopore trajectories and movies

[Results, MP4 movies and downloadable inputs](nanopore/results/2026-09-17/README.md):
dA–7AHL completed its 5 ns production segment; dA–3B07 is published as an
explicitly labelled 0–4.825 ns in-progress preview. Neither exported trajectory
shows complete DNA passage. These are real-coordinate movies of pore-mouth
interactions and conformational changes, not illustrations of assumed passage.

The published dA trajectories use NAMD 2.14 CUDA; see the
[run-specific protocol](nanopore/results/2026-09-17/PROTOCOL.md).

## Protocols

- [Nanopore protocol / 过孔模拟流程](nanopore/PROTOCOL.md): system preparation,
  force fields, restraints, stage durations, voltage, Slurm execution and checks.
- [ssDNA dataset](ssdna/README.md): simulation scripts, results and analysis.
- [Job status](PORE_JOB_STATUS.md): nanopore workflow status.

The DNA-only preparation is 50,000 minimization steps → 0.5 ns heating →
2 ns restrained equilibration → 2.5 ns weak-restraint equilibration →
5 ns unrestrained equilibration at 2 fs, 293 K and 1 atm. dA then has 50 ns
production; dT currently has one validated 5 ns production block.

The nanopore sequence is 50,000 minimization steps → 0.5 ns heating →
4 ns staged equilibration → 5 ns electric-field production at 1 fs.
The pore boxes are fixed; this is not the pressure-coupled DNA-only protocol.

## Repository contents

- `ssdna/polyda40-3m-kcl/{scripts,results}/`: dA run inputs/code and outputs.
- `ssdna/polydt40-3m-kcl/{scripts,results}/`: dT code and outputs.
- `ssdna/comparison/{scripts,results}/`: cross-sequence analysis and figures.
- `nanopore/<sequence>-3m-kcl-<pore>/run/`: DNA placement, system settings and
  ion-placement checks; `nanopore/common/` contains the shared workflow.
- `nanopore/7ahl-template/`, `nanopore/3b07-template/`: processed pore/membrane/water templates.
- `data-manifest/`: filenames, sizes and SHA-256 checksums for omitted DCD and
  large solvated-system files.


## Data quality

Coordinates are reconstructed through PSF covalent bonds before analysis. The
dT trajectory remained at least 85.886 Å from its nearest periodic image and
passed all-frame bond and image-distance checks.

Each sequence currently has one usable trajectory, not an independent-replica
uncertainty estimate. The equal windows are also at different trajectory ages:
dA uses its last 5 ns, whereas dT uses its first 5 ns production block.

## Software

NAMD 2.14 CUDA; VMD 1.9.2; CHARMM nucleic-acid, protein, lipid and water/ion
parameters as specified in the protocols; Python 3 with NumPy/SciPy; Slurm.
Force-field installation paths and the NAMD module name are cluster-specific.
