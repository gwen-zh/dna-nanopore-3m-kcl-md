# Pore workflow status

Snapshot: 2026-09-09 (America/New_York)

| System | Build | GPU preparation | GPU production | Analysis |
|---|---:|---:|---:|---:|
| poly(dA)40–7AHL | 57311447 | 57311448 | 57311449 | 57311450 |
| poly(dT)40–7AHL | 57311451 | 57311452 | 57311453 | 57311454 |
| poly(dA)40–3B07 | 57311455 | 57311456 | 57311457 | 57311458 |
| poly(dT)40–3B07 | 57311459 | 57311460 | 57311461 | 57311462 |

At the 10:20 EDT snapshot, both 7AHL builds have completed and passed final
PSF/PDB, ion-content, clash, and box-clearance validation. Their GPU preparation
jobs are pending for `dept_gpu` priority. Both larger 3B07 builds are still
actively placing 3 M KCl ions. Every NAMD preparation and production job is
pending on its corresponding successful predecessor and requests one GPU from
`dept_gpu`; the four GPU branches may run concurrently. Analysis jobs run after
successful production. The active job limits are 2 days for preparation and 3
days for production, based on measured performance with additional margin.

The formal build validates final PSF/PDB and KCl content and performs a
deterministic protein-clash and periodic-box-clearance check for the folded
high-salt DNA.
