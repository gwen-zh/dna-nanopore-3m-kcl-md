# ssDNA and nanopore simulations in 3 M KCl

This repository contains poly(dA)40 and poly(dT)40 DNA-only simulations plus
DNA–7AHL and DNA–3B07 nanopore workflows, with executable inputs, validation
data, structural analysis and figures.

## ssDNA datasets

Both dA40 and dT40 contain 50 ns of unrestrained production after the same
minimization, heating and 9.5 ns equilibration schedule at 293 K. The matched
structural comparison uses 45.01–50.00 ns from each sequence.

- [ssDNA datasets and figure previews](ssdna/README.md)
- [Complete dA40/dT40 analysis](ssdna/comparison/results/README.md)
- [dT40 50 ns results](ssdna/polydt40-3m-kcl/results/README.md)

## Nanopore workflows

Four 3 M KCl systems cover dA40 and dT40 with 7AHL and 3B07. Each combination
has its own `script/` and `result/` directories under
[`nanopore/`](nanopore/README.md). The [nanopore protocol](nanopore/PROTOCOL.md)
describes the applied field, ordered 5′-first pathway, interaction analysis and
validation checks.

Large raw trajectories and solvated restart files are omitted from Git. Compact
analysis tables, structures, figures and trajectory movies are included for
completed pathways.
