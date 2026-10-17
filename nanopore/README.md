# ssDNA translocation through 3B07 and 7AHL

This directory is organized as four sequence–pore combinations. Each combination
contains the simulation and analysis code in `script/` and its validated outputs
in `result/`.

| ssDNA | Pore | Result |
|---|---|---|
| poly(dA)40 | 3B07 | complete 5′-first cis-to-trans pathway; analysis and movie included |
| poly(dA)40 | 7AHL | complete 5′-first cis-to-trans pathway; analysis and movie included |
| poly(dT)40 | 3B07 | complete 5′-first cis-to-trans pathway; analysis and movie included |
| poly(dT)40 | 7AHL | complete 5′-first cis-to-trans pathway; analysis and movie included |

- [poly(dA)40–3B07](polyda40-3b07/README.md)
- [poly(dA)40–7AHL](polyda40-7ahl/README.md)
- [poly(dT)40–3B07](polydt40-3b07/README.md)
- [poly(dT)40–7AHL](polydt40-7ahl/README.md)
- [Simulation and analysis protocol](PROTOCOL.md)

## Completed-pathway overview

| Pathway | Saved structures | Path time (ns) | Mean contacts at 3.5 Å | Mean H-bond candidates | Mean K+ bridges | Mean constriction water bridges |
|---|---:|---:|---:|---:|---:|---:|
| poly(dA)40–3B07 | 386 | 9.625 | 6.05 | 1.77 | 0.26 | 2.67 |
| poly(dA)40–7AHL | 430 | 10.725 | 13.77 | 3.78 | 0.44 | 7.93 |
| poly(dT)40–3B07 | 376 | 9.375 | 9.74 | 1.72 | 1.05 | 2.71 |
| poly(dT)40–7AHL | 433 | 10.800 | 12.06 | 2.02 | 1.09 | 6.48 |

All four sequence–pore combinations are now represented, which permits matched
within-pore descriptive comparisons. Each row is still one driven pathway, so
differences should be treated as interaction hypotheses rather than converged
sequence-specific ensemble averages.

## Interpretation boundary

These are high-field, ordered steered-MD pathways designed to obtain complete
translocation geometries and identify contacts along the pore. They support
residue-level interaction hypotheses for the experimental current signal, but
they do not directly predict the experimental voltage dwell time, unbiased
translocation rate, binding free energy or blockade current. Current prediction
requires a separate ion-current calculation or conductance workflow.

Large solvated trajectories and restart checkpoints are retained in the cluster
workspace and are not duplicated in Git. Each completed `result/` contains the
full rendered pathway, structures, plots, interaction tables and validation
metadata needed to inspect the published analysis.
