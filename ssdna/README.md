# ssDNA simulations

| Dataset | Scripts | Results |
|---|---|---|
| poly(dA)40 in 3 M KCl | [`polyda40-3m-kcl/scripts/`](polyda40-3m-kcl/scripts/README.md) | [`polyda40-3m-kcl/results/`](polyda40-3m-kcl/results/README.md) |
| poly(dT)40 in 3 M KCl | [`polydt40-3m-kcl/scripts/`](polydt40-3m-kcl/scripts/README.md) | [`polydt40-3m-kcl/results/`](polydt40-3m-kcl/results/README.md) |
| dA40/dT40 comparison | [`comparison/scripts/`](comparison/scripts/README.md) | [`comparison/results/`](comparison/results/README.md) |

Both sequences contain 50 ns of production after the same minimization,
heating and 9.5 ns equilibration schedule. The structural comparison uses the
matched 45.01–50.00 ns window with 500 frames per sequence.

![Global structural metrics](comparison/results/01_global_structure.png)

![Per-residue structural profiles](comparison/results/02_residue_profiles.png)

The analysis includes size and shape, stacking, χ/syn state, outward-base
orientation, nonlocal contacts, hydrogen bonds, ion association, hydration,
representative structures and periodic-image validation.
