# ssDNA simulations

The current structural comparison, all plots, per-frame/per-residue data and
representative PDBs are in
[`analysis-2026-09-21/`](analysis-2026-09-21/README.md).

| Dataset | Current use |
|---|---|
| `polyda40-3m-kcl/50ns/` | Valid dA source; comparison uses 45.01–50.00 ns. |
| `polydt40-3m-kcl/repaired-5ns/` | Valid image-safe dT source; comparison uses 0.01–5.00 ns. |
| `polydt40-3m-kcl/50ns/` | Archived invalid source; do not use for structural conclusions. |

The common preparation is 50,000 minimization steps, 0.5 ns heating, 2 ns
restrained equilibration, 2.5 ns weak-restraint equilibration and 5 ns
unrestrained equilibration, with a 2 fs timestep. dA production is a continuous
20 ns segment plus 30 ns continuation. Repaired dT currently contains one
validated 5 ns production block after the same preparation schedule.

DCD filenames, sizes and hashes are in
[the original manifest](../data-manifest/ssdna_dcd.sha256) and
[the repaired-dT manifest](../data-manifest/repaired_dt_5ns.sha256).

## Analysis status

The 2026-09-21 results are recalculated from 500 frames per sequence after
PSF-bond reconstruction. They include Rg, end-to-end/contour length, shape,
stacking, χ/syn occupancy, outward base orientation, nonlocal contacts,
hydrogen bonds, RMSF, ion association and hydration.

The older CSVs and representative structures inside `polydt40-3m-kcl/50ns/`
are retained only for provenance. They are superseded and must not be mixed
with the repaired results.

For nanopore construction, use the validated donor files recorded by the
[nanopore protocol](../nanopore/PROTOCOL.md), not an older representative export.

Log normalization changes only the working-directory header to a repository
path. Energy, timing and other numerical records are unchanged.
