# poly(dT)40 in 3 M KCl

This is the completed, image-safe poly(dT)40 DNA-only branch used by the
2026-09-21 structural analysis. Slurm job 57330935 ran on `dept_gpu` with one
non-L40 GPU and NAMD 2.14 CUDA, and ended normally at 2,500,000 production
steps on 2026-09-19 20:59 EDT.

The starting DNA was PSF-bond reconstructed from frame 1990 (19.91 ns) of the
earlier trajectory, centered, and resolvated in an explicit 180 Å cubic box.
It was not stretched or regenerated. The new system contains 526,789 atoms,
168,075 waters, 10,662 K⁺ and 10,623 Cl⁻ ions.

## Protocol

| Stage | Steps | Time | Backbone restraint scale |
|---|---:|---:|---:|
| Minimization | 50,000 | — | 1.0 |
| Heating, 50→293 K | 250,000 | 0.5 ns | 1.0 |
| NPT equilibration 1 | 1,000,000 | 2.0 ns | 0.5 |
| NPT equilibration 2 | 1,250,000 | 2.5 ns | 0.1 |
| NPT equilibration 3 | 2,500,000 | 5.0 ns | none |
| Production | 2,500,000 | 5.0 ns | none |

All dynamics use a 2 fs timestep, 293 K, 1.01325 bar, PME, a 12 Å cutoff,
10 Å switching distance, 14 Å pair list and `rigidBonds all`. Production saves
one DCD frame every 5,000 steps (10 ps).

## Validation

- 500/500 production frames were read successfully.
- Maximum reconstructed DNA bond: 1.703735 Å.
- Minimum DNA-to-periodic-image distance: 85.886434 Å.
- Endpoint Rg / end-to-end distance: 27.1185 / 73.4963 Å.
- Mean last-1-ns NAMD temperature: 292.246 K (501 energy records).
- Mean production box concentrations: 3.2778 M K⁺ and 3.2658 M Cl⁻. K⁺
  includes the counterions needed to neutralize the DNA.

Stage endpoint checks are in [`results/validation/`](results/validation/), and the full
production geometry trace is
[`results/validation/prod5.metrics.csv`](results/validation/prod5.metrics.csv).
The exact code snapshot used by the job is in [`scripts/`](scripts/README.md).
`scripts/box180/dna.pdb`, `dna.psf` and `donor.json` document the reconstructed
DNA donor.

[`results/validation/production_all_frame_summary.json`](results/validation/production_all_frame_summary.json)
collects the all-frame ranges, last-1-ns temperature and DCD hash in one file.

The large solvated PSF/PDB, binary checkpoints and DCD are not committed.
Their sizes and SHA-256 hashes are recorded in
[`data-manifest/polydt40_dcd.sha256`](../../data-manifest/polydt40_dcd.sha256).

The 5 ns production is retained because it is the valid dT source used by the
current structural comparison.
