# poly(dA)40 in 3 M KCl

- [`scripts/`](scripts/README.md): NAMD inputs, Slurm launchers, force-field
  supplement, starting system and legacy VMD analysis scripts.
- [`results/`](results/README.md): stage logs, restart checkpoints and archived
  sequence-specific analysis outputs.

The production trajectory is 50 ns at 293 K and 1.01325 bar with a 2 fs
timestep. The validated dA40/dT40 comparison uses its 45.01–50.00 ns window and
is stored under [`../comparison/results/`](../comparison/results/README.md).
