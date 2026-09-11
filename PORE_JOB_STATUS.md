# Nanopore job status

Snapshot: 2026-09-11 09:37 EDT (America/New_York).
All four preparation jobs are RUNNING on `dept_gpu`, one GPU per system.

| System | Build (completed) | Preparation | Node | Observed minimization step | Production (pending) | Analysis (pending) |
|---|---:|---:|---|---:|---:|---:|
| poly(dA)40–7AHL | 57318990 | 57318979 | g007 | 43,848 / 50,000 | 57311449 | 57311450 |
| poly(dT)40–7AHL | 57318991 | 57318981 | g012 | 41,020 / 50,000 | 57311453 | 57311454 |
| poly(dA)40–3B07 | 57318992 | 57318983 | g009 | 23,633 / 50,000 | 57311457 | 57311458 |
| poly(dT)40–3B07 | 57318993 | 57318985 | g010 | 28,679 / 50,000 | 57311461 | 57311462 |

The sequence is minimization, 0.5 ns heating, 4 ns staged equilibration,
and 5 ns production at 1 fs. Production waits for successful preparation;
it is not running yet at this snapshot. Preparation/production request
one GPU, 8 CPU cores and 16 GB RAM, at normal priority. Maximum simultaneous
GPU demand for these four chains is four. Other user jobs were not modified.

This is a timestamped record, not live monitoring or a completion forecast.
The analysis jobs are queued, but their output will need the periodic-image,
time-axis and geometry audit described in [the protocol](nanopore/PROTOCOL.md).

The published inputs include the repaired intact DNA donors, protein-parameter
fallback and ion-placement validation used for these builds. Repository path
renaming was done in a separate clone; running cluster directories were not
renamed and no new simulation jobs were submitted for this repository update.
