# Formal nanopore workflows

This directory contains only the four current formal translocation systems:

- poly(dA)40 through 7AHL in 3 M KCl
- poly(dT)40 through 7AHL in 3 M KCl
- poly(dA)40 through 3B07 in 3 M KCl
- poly(dT)40 through 3B07 in 3 M KCl

The placed DNA structures come from the final 10 ns of the corresponding
50 ns ssDNA trajectory. `7ahl-template` and `3b07-template` contain clean
protein/membrane/water templates; no preliminary translocation trajectory is
included.

The formal protocol is 50,000 minimization steps, 0.5 ns heating, 4 ns staged
equilibration, then 5 ns electric-field production at a 1 fs timestep and a
1.0 V box-spanning potential. Each NAMD branch requests one `dept_gpu` GPU.
