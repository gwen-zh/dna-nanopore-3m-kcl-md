# Equal-window dA40/dT40 structural analysis in 3 M KCl

This directory contains the validated, periodic-image-aware comparison generated
on 2026-09-21. Coordinates are reconstructed through PSF covalent bonds before
any molecular geometry is measured.

## Compared windows

| Sequence | Trajectory used | Frames | Source time | Window length |
|---|---|---:|---:|---:|
| poly(dA)40 | mature 50 ns run, final segment | 2500–2999 | 45.01–50.00 ns | 5.00 ns |
| poly(dT)40 | repaired 180 Å-box production | 0–499 | 0.01–5.00 ns | 5.00 ns |

Both DCDs contain one frame every 10 ps. The dT trajectory followed 0.5 ns
heating and 9.5 ns equilibration. The older dT 50 ns trajectory is not used
because it contains genuine periodic-image contacts. Exact input hashes and
selected frame indices are in [`validation.json`](validation.json).

## Main result

Values are mean ± frame-to-frame standard deviation over 500 frames. Frames are
autocorrelated, so these standard deviations are descriptive and are not
independent-replica uncertainty estimates.

| Metric | dA40 | dT40 |
|---|---:|---:|
| Radius of gyration, Å | 25.13 ± 0.80 | 27.32 ± 1.05 |
| 5′–3′ C1′ distance, Å | 53.57 ± 2.66 | 70.83 ± 4.38 |
| C1′ contour length, Å | 264.90 ± 2.49 | 293.72 ± 3.12 |
| End-to-end / contour | 0.202 ± 0.010 | 0.241 ± 0.014 |
| Shape anisotropy | 0.537 ± 0.034 | 0.709 ± 0.039 |
| Adjacent stacked pairs | 11.84 ± 1.83 | 4.01 ± 0.75 |
| Unstacked-from-neighbors fraction | 0.509 ± 0.064 | 0.841 ± 0.030 |
| syn fraction | 0.648 ± 0.061 | 0.278 ± 0.062 |
| Mean outward base projection | 0.221 ± 0.039 | 0.004 ± 0.034 |
| Nonlocal residue-contact pairs | 30.19 ± 1.80 | 29.06 ± 2.00 |
| Base–base hydrogen bonds | 3.89 ± 1.34 | 0.78 ± 0.43 |
| All inter-residue hydrogen bonds | 7.86 ± 1.73 | 4.49 ± 1.25 |
| K⁺ within 3.5 Å of DNA | 36.35 ± 4.57 | 29.67 ± 3.76 |
| Water O within 3.5 Å of DNA | 473.73 ± 11.29 | 465.06 ± 10.55 |

The global measures consistently show that dA40 is more compact and less
rod-like than dT40 in these windows: its Rg is about 8% lower, end-to-end
distance about 24% lower, and shape anisotropy about 24% lower.
Across all ten 0.5 ns blocks, dA40 has lower end-to-end distance and anisotropy,
higher syn fraction and outward projection, and lower unstacked-neighbor
fraction than dT40; these directions are therefore not caused by one brief
frame excursion.

“Base flipping” is not one unique observable. Two orientation measures support
more prominent base reorientation in dA40: its syn occupancy is 64.8% versus
27.8% for dT40, and its mean radial projection is more positive. dA residues
1, 2, 7–9, 12, 25, 28, 31, 33 and 40 have syn occupancy above 0.90; residues
8–10, 13, 18, 23, 25, 35 and 40 have mean outward projection above 0.70.

At the same time, dT40 has much more loss of adjacent stacking, while dA40
retains more adjacent stacks and more inter-base hydrogen bonds. Thus
“outward/syn” and “unstacked” should not be treated as synonyms: the present
data support a compact, substantially stacked dA ensemble with more syn/outward
bases, and a more extended dT ensemble with greater local stacking disorder.

![Global structure time series](01_global_structure.png)

![Per-residue profiles](02_residue_profiles.png)

![Base-state heat maps](03_base_state_heatmaps.png)

In the state heat maps, purple is 0 and yellow is 1. In the top row, 1 means
syn; in the bottom row, 1 means no strict stack with either sequence neighbor.

## Operational definitions

- Rg is mass-weighted over DNA heavy atoms. Shape anisotropy ranges from 0 for
  an isotropic object to 1 for an ideal line.
- End-to-end distance joins residue 1 and 40 C1′ atoms; contour length sums all
  adjacent C1′ distances.
- A strict base stack requires centroid distance ≤5.5 Å, plane-normal alignment
  within 45°, normal separation 2.0–4.5 Å and lateral displacement ≤3.5 Å.
- χ uses O4′–C1′–N9–C4 for adenine and O4′–C1′–N1–C2 for thymine. `|χ| < 90°`
  is classified as syn.
- Outward projection is the dot product of the unit C1′→base-centroid vector
  and the unit molecular-center→C1′ vector. Positive values point outward.
- A nonlocal contact is any heavy-atom distance ≤4.5 Å between residues
  separated by at least three positions.
- Hydrogen bonds use donor–acceptor distance ≤3.5 Å and D–H–A angle ≥150°.
- Ion and water contacts are counted within 3.5 Å of any DNA heavy atom every
  100 ps using periodic nearest-neighbor distances.

The largest reconstructed DNA bond over all analyzed frames is 1.707 Å for
dA40 and 1.704 Å for dT40. The repaired dT production has a minimum DNA-to-image
distance of 85.886 Å, so no short-range periodic-image contact occurs.

## Files

- [`summary.csv`](summary.csv): one-row-per-sequence global summary.
- [`timeseries.csv`](timeseries.csv): 1,000 frame records for global metrics.
- [`block_means_0p5ns.csv`](block_means_0p5ns.csv): 0.5 ns block means.
- [`per_residue_summary.csv`](per_residue_summary.csv): 40 residue summaries per sequence.
- [`per_residue_timeseries.csv`](per_residue_timeseries.csv): all 40,000 residue/frame records.
- [`environment_timeseries.csv`](environment_timeseries.csv): ion, water and box concentration data.
- `dA40_*_occupancy.csv`, `dT40_*_occupancy.csv`: 40 × 40 contact, stack and hydrogen-bond occupancies.
- `dA40_representative.pdb`, `dT40_representative.pdb`: real sampled representative frames.
- `01`–`07` `.png` and `.pdf`: raster and vector versions of every plot.
- [`../../scripts/analyze_equal_window.py`](../../scripts/analyze_equal_window.py):
  standalone Python analysis and plotting source.

## Limits

This is an equal-duration comparison, but not an equal-trajectory-age
comparison: dA uses a late 5 ns window whereas dT currently has only
its first 5 ns production block. There is one trajectory per sequence, no
independent replica, and no claim of ensemble convergence. More repaired dT
sampling and matched replicas are needed before assigning population-level
free-energy differences.
