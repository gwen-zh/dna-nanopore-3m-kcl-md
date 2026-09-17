# Nanopore job status

Snapshot: 2026-09-17 16:32 EDT (America/New_York).

| System / task | Job | State at snapshot | Published result |
|---|---:|---|---|
| poly(dA)40–7AHL | 57323201 | 5 ns production completed; endpoint geometry passed | [5.000 ns result and movie](nanopore/results/2026-09-17/polyda40-3m-kcl-7ahl/) |
| poly(dA)40–3B07 | 57323197 | RUNNING on dept_gpu | [0–4.825 ns in-progress preview](nanopore/results/2026-09-17/polyda40-3m-kcl-3b07-preview/) |
| New dT DNA-only, larger box | 57330935 | RUNNING on dept_gpu | Not yet an approved new pore donor |
| New poly(dT)40–7AHL | 57334241 | PENDING, afterok:57330935 | No new pore trajectory yet |
| New poly(dT)40–3B07 | 57334242 | PENDING, afterok:57330935 | No new pore trajectory yet |

The dA–7AHL log ended normally at step 5000000 on September 17 at 09:41 EDT.
Completion means the scheduled 5 ns segment, not complete DNA translocation.
Neither exported dA trajectory shows full passage; see the [results](nanopore/results/2026-09-17/README.md).
The 3B07 movie is a fixed, labelled snapshot, not a continuously updated feed.

Two GPUs were allocated to this research at the snapshot (dA–3B07 and new dT
DNA-only); waiting dependencies allocate none. These MD jobs use dept_gpu,
one non-L40 GPU and eight CPU workers per job. Other user jobs are untouched.
The actual dA production engine is NAMD 2.14 CUDA, not the earlier NAMD 3.0.2 route.

The two previous dT pore jobs were stopped to replace their periodic-image-biased
DNA donor. New dT pore starts require successful bulk sampling, geometry/image
and drift checks, safe placement and membrane checks; they are not automatic
resumes of the old dT trajectories. No dT result is included in this dA release.

Video rendering and analysis use CPUs only and do not interrupt simulations.
