# Scripts

- `ordered_translocation.py`: advances successive three-nucleotide fronts and
  validates contiguous 5′-first passage.
- `run.py`: shared stage-generation and geometry helpers captured with the job.
- `stage.namd`: NAMD configuration used for each stage.
- `submit.sbatch`: one-GPU `dept_gpu` submission script.
- `analyze_translocation.py`: reconstructs periodic molecules and computes the
  structural and interaction observables documented in the protocol.
- `trajectory.py`: PSF, DCD and coordinate I/O helpers.
- `make_movie.py` and `render.tcl`: VMD/FFmpeg trajectory rendering.
- `requirements.txt`: Python analysis dependencies.

The simulation scripts are the code snapshot associated with this pathway and
retain cluster-workspace imports and paths. They are provided for method review;
those paths and system/checkpoint locations must be adapted before reuse in a
different workspace. The analysis scripts operate on the segmented trajectory
and require a strict 40/40 workflow-validation record.
