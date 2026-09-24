# poly(dA)40 in 3 M KCl

- [`scripts/`](scripts/README.md): NAMD inputs, Slurm launchers, force-field
  supplement and starting system.
- [`results/`](results/README.md): stage logs and restart checkpoints for the
  50 ns simulation.

The production trajectory is 50 ns at 293 K and 1.01325 bar with a 2 fs
timestep. The validated dA40/dT40 comparison uses its 45.01–50.00 ns window and
is stored under [`../comparison/results/`](../comparison/results/README.md).
