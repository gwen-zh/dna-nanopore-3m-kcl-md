# Pore trajectory analysis and actual-coordinate movies

These scripts produced the [2026-09-17 dA release](../results/2026-09-17/README.md).
They never submit or alter a simulation job. Existing result directories are
not overwritten by the numerical-analysis entry point.

`trajectory.py` reads PSF/NAMD DCD data, reconstructs bonded fragments and writes
display-only subsets. `analyze.py` takes a cluster project root, pore and fresh
output directory; the two source job IDs are deliberately explicit for this dataset.

```bash
python analyze.py --workspace /path/to/cluster/project --pore 7ahl --out /path/to/new/results
VMD_BIN=/path/to/vmd python movie.py /path/to/new/results
python -m unittest discover -s . -p test_results.py -v
```

To recreate an already published movie, decompress `visualization.dcd.gz` beside
the provided visualization PSF and run `movie.py` on that result directory.
VMD uses two CPU threads and no GPU. The output is checked for blank/frozen
molecular frames, compared against independent VMD geometry/contact measurements,
and completely decoded after H.264 encoding. No coordinates are interpolated.

The seven unit tests cover covalent PBC reconstruction, disconnected-fragment
rejection, proper rigid rotations, partial DCD tails and time axes, display DCD
cell flags, fixed-payload hashing and PDB column alignment. All seven passed.
See [methods](../results/2026-09-17/METHODS.md) for scientific definitions and limitations.
