# Shared nanopore workflow

See [the protocol](../PROTOCOL.md) for parameters and execution order.

- `extract_intact_donor.tcl`: extract one raw DCD frame, join DNA through PSF bonds and check bond lengths.
- `polyda40_donor_intact.pdb`, `polydt40_donor_intact.pdb`: validated donor coordinates.
- `place_dna.py`, `build_pore.tcl`: placement, DNA topology, displaced-water removal and solvation.
- `fast_autoionize.tcl`, `fast_ion_positions.py`: retain stock ion counts/topology
  while using periodic neighbor queries for site selection.
- `validate_bonds.py`: pre-GPU coordinate gate.
- `make_restraints.tcl`: reference coordinates and per-atom restraint coefficients.
- `min.namd`, `heat.namd`, `eq1.namd`, `eq2.namd`, `eq3.namd`, `prod5.namd`: sequential MD stages.
- `build_job.sbatch`, `prep_job.sbatch`, `prod_job.sbatch`, `analyze_job.sbatch`: Slurm jobs.
- `finish_ionization.tcl`, `finish_build_job.sbatch`: finish an already-solvated build.
- `analyze_pore.tcl`: trajectory-analysis implementation pending the image/time-axis audit described in the protocol.

All NAMD stages enter the exported `SYSTEM_DIR` before reading inputs.
Preparation and production request one GPU and 16 GB RAM on `dept_gpu`.

These scripts reflect the running workflow, with repository-relative paths.
No new cluster jobs are launched by cloning or reading this repository.
