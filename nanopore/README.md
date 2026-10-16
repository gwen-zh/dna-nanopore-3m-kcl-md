# ssDNA translocation through 3B07 and 7AHL

This directory is organized as four sequence–pore combinations. Each combination
contains the simulation and analysis code in `script/` and its validated outputs
in `result/`.

| ssDNA | Pore | Result |
|---|---|---|
| poly(dA)40 | 3B07 | complete 5′-first cis-to-trans pathway; analysis and movie included |
| poly(dA)40 | 7AHL | simulation running on `dept_gpu`; result directory reserved |
| poly(dT)40 | 3B07 | simulation running on `dept_gpu`; result directory reserved |
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
| poly(dT)40–7AHL | 433 | 10.800 | 12.06 | 2.02 | 1.09 | 6.48 |

The two rows involve different pores as well as different sequences, so their
numerical difference is not a controlled dA-versus-dT sequence effect. The two
remaining pathways will supply the matched within-pore comparisons.

## Interpretation boundary

These are high-field, ordered steered-MD pathways designed to obtain complete
translocation geometries and identify contacts along the pore. They support
residue-level interaction hypotheses for the experimental current signal, but
they do not directly predict the experimental −120 mV dwell time, unbiased
translocation rate, binding free energy or blockade current. Current prediction
requires a separate ion-current calculation or conductance workflow.

Large solvated trajectories and restart checkpoints are retained in the cluster
workspace and are not duplicated in Git. Each completed `result/` contains the
full rendered pathway, structures, plots, interaction tables and validation
metadata needed to inspect the published analysis.
