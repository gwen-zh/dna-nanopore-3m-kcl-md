# Nanopore simulations

Four systems in nominal 3 M KCl:

- poly(dA)40–7AHL
- poly(dT)40–7AHL
- poly(dA)40–3B07
- poly(dT)40–3B07

Read [PROTOCOL.md](PROTOCOL.md) for the complete method, stage table, force
fields, electric-field conversion and submission example.

`common/` contains shared inputs and build/validation helpers.
`<sequence>-3m-kcl-<pore>/run/` contains system-specific settings and the
corrected DNA placement. The two `*-template/` directories contain processed
protein/membrane/water templates without ssDNA or ions.

DNA donors were extracted from original frames in the 40–50 ns portion of
the corresponding ssDNA trajectory, joined through their bond connectivity,
and placed independently for each pore. These are not guaranteed median-Rg
frames of a corrected analysis.

Preparation is followed by 5 ns production at 1 fs. All NAMD branches use
one `dept_gpu` GPU each, with at most four concurrently. See the
[timestamped status](../PORE_JOB_STATUS.md); a configured production stage is
not a claim that DNA has completely translocated.
