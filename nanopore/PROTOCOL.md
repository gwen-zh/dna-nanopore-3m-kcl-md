# Nanopore translocation protocol

## Systems

- ssDNA: poly(dA)40 or poly(dT)40 after 50 ns DNA-only sampling in nominal
  3 M KCl. The selected input frames are 45.11 ns for dA and 47.71 ns for dT.
- Pores: 3B07 and 7AHL in explicit membrane, water and ions.
- Temperature: 293 K.
- Integration timestep: 1 fs.
- Direction: cis to trans, with the 5′ end entering first.

## Driving protocol

The electric-field parameter is `eField 0.0 0.0 0.1695`. With the actual box
length this corresponds to approximately 1.91 V for 3B07 and 1.84 V for 7AHL.
This field is deliberately much larger than the experimental −120 mV bias so
that a complete atomistic pathway can be sampled within accessible compute time.

Ordered steered MD advances the next contiguous three-nucleotide front toward
the trans chamber at 50 Å/ns. Once a front has crossed, residues 1–3 are held in
the trans chamber to prevent reversal while the remaining chain is advanced.
Each stage is accepted only when the crossed residues form a contiguous 5′-first
prefix. Completion requires all 40 C1′ atoms below the trans-side reference and
no out-of-order crossing.

## Analysis

Every saved frame is reconstructed through the PSF covalent-bond graph before
geometry is measured. The analysis reports:

- nucleotide progress relative to the cis entrance, constriction and trans exit;
- DNA–protein atom contacts and residue-pair contact occupancy;
- constriction contacts and per-nucleotide contact occupancy;
- geometric hydrogen-bond candidates;
- potassium-mediated DNA–acidic-residue bridges;
- DNA–constriction water bridges and DNA hydration;
- base orientation and base/backbone contact classes;
- instantaneous SMD spring-force traces;
- minimum DNA–protein separation and maximum DNA bond length;
- start and endpoint structures, plots and an annotated trajectory movie.

Hydrogen bonds, ion bridges and water bridges are geometric assignments. Contact
occupancy is the fraction of saved frames in this one nonequilibrium pathway; it
is not an equilibrium probability or a free-energy estimate.

## Reproducibility and validation

The `script/` directory in each system contains the job script, NAMD stage
configuration, ordered-translocation driver, analysis code and movie renderer.
The simulation driver writes stage-level order and geometry checks. Analysis is
released only after the workflow records strict 40/40 passage with no
out-of-order event. Raw DCD segments and restart checkpoints remain on the
cluster because of their size.
