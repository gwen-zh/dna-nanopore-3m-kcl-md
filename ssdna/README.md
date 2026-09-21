# ssDNA simulations

The ssDNA files are organized by sequence. Each retained sequence has its own
`scripts/` and `results/` directory; cross-sequence analysis is kept separately
under `comparison/`.

| Dataset | Scripts | Results |
|---|---|---|
| poly(dA)40 in 3 M KCl | [`polyda40-3m-kcl/scripts/`](polyda40-3m-kcl/scripts/README.md) | [`polyda40-3m-kcl/results/`](polyda40-3m-kcl/results/README.md) |
| poly(dT)40 in 3 M KCl | [`polydt40-3m-kcl/scripts/`](polydt40-3m-kcl/scripts/README.md) | [`polydt40-3m-kcl/results/`](polydt40-3m-kcl/results/README.md) |
| dA40/dT40 comparison | [`comparison/scripts/`](comparison/scripts/README.md) | [`comparison/results/`](comparison/results/README.md) |

The common preparation is 50,000 minimization steps, 0.5 ns heating, 2 ns
restrained equilibration, 2.5 ns weak-restraint equilibration and 5 ns
unrestrained equilibration, with a 2 fs timestep at 293 K. dA has 50 ns
production. The retained, image-safe dT dataset currently has one validated
5 ns production block after the same preparation schedule.

The dT 5 ns block is retained because it is the dT trajectory used by the
validated comparison. The older dT 50 ns export had periodic-image contacts and
has been removed from the current repository tree.

Large omitted trajectories are recorded in
[`polyda40_dcd.sha256`](../data-manifest/polyda40_dcd.sha256) and
[`polydt40_dcd.sha256`](../data-manifest/polydt40_dcd.sha256).

The validated comparison uses 500 frames per sequence after PSF-bond
reconstruction. It covers Rg, end-to-end and contour lengths, shape, stacking,
χ/syn occupancy, outward base orientation, contacts, hydrogen bonds, RMSF,
ion association and hydration. There is one usable trajectory per sequence,
so frame variability is descriptive rather than independent-replica uncertainty.
