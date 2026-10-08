# poly(dT)40 in 3 M KCl

This dataset contains 50 ns of unrestrained production after minimization,
heating and three equilibration stages. The simulation used a 2 fs timestep,
293 K, 1.01325 bar, PME, a 12 Å cutoff, 10 Å switching distance and a 14 Å
pair list. One non-L40 `dept_gpu` GPU ran NAMD 2.14 CUDA.

| Stage | Steps | Time |
|---|---:|---:|
| Minimization | 50,000 | — |
| Heating, 50→293 K | 250,000 | 0.5 ns |
| Restrained equilibration | 1,000,000 | 2.0 ns |
| Weak-restraint equilibration | 1,250,000 | 2.5 ns |
| Unrestrained equilibration | 2,500,000 | 5.0 ns |
| Unrestrained production | 25,000,000 | 50.0 ns |

All 5,000 production frames passed bonded-geometry and periodic-image checks.
The maximum reconstructed DNA bond is 1.737 Å
and the minimum periodic-image separation is
85.493 Å.

- [Simulation and analysis results](results/README.md)
- [Executed production scripts](scripts/README.md)
- [Matched dA40/dT40 structural comparison](../comparison/results/README.md)
