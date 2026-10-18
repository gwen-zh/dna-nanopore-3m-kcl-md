# Four-pathway translocation comparison

| Mean pathway metric | [dA40–3B07](polyda40-3b07/result/README.md) | [dT40–3B07](polydt40-3b07/result/README.md) | [dA40–7AHL](polyda40-7ahl/result/README.md) | [dT40–7AHL](polydt40-7ahl/result/README.md) |
|---|---:|---:|---:|---:|
| DNA–protein atom contacts at 3.5 Å | 6.05 | 9.74 | 13.77 | 12.06 |
| Constriction residue-pair contacts | 0.60 | 0.66 | 2.30 | 2.18 |
| Hydrogen-bond candidates | 1.77 | 1.72 | 3.78 | 2.02 |
| K+ bridges | 0.26 | 1.05 | 0.44 | 1.09 |
| Constriction water bridges | 2.67 | 2.71 | 7.93 | 6.48 |
| Adjacent base stacks | 12.58 | 4.98 | 14.56 | 4.99 |
| DNA radius of gyration (Å) | 41.59 | 39.65 | 43.12 | 40.25 |
| SMD force magnitude (pN) | 638 | 670 | 1055 | 1061 |

## Main conclusions

### Sequence structure

- dA40 retains 2.5–2.9 times more adjacent-base stacking than dT40 in both
  pores and has the larger radius of gyration.
- dA40 therefore passes through the pores as the more extended and ordered
  chain, whereas dT40 is more flexible and adapts more readily to the pore wall.

### Sequence discrimination in 3B07

- dT40 forms 61% more DNA–protein atom contacts, contacts 40% more nucleotides
  per frame and forms approximately four times as many K+-mediated bridges as
  dA40.
- dA40 and dT40 have similar hydrogen-bond, constriction-contact and water-bridge
  counts at the narrowest region. The main 3B07 sequence difference therefore
  comes from broader wall contact and ion-mediated anchoring by dT40.
- GLN112 and LYS108 are the main shared contact sites. dT40 additionally makes
  broader contacts involving GLN118, ASN144 and other lumen residues.

### Sequence discrimination in 7AHL

- dA40 forms 14% more DNA–protein atom contacts, 1.87 times as many candidate
  hydrogen bonds and 22% more constriction water bridges than dT40.
- dT40 forms 2.45 times as many K+-mediated bridges as dA40.
- The two sequences have nearly identical mean SMD force magnitudes in 7AHL,
  but their interaction chemistry differs: dA40 is dominated by direct
  hydrogen-bond and water-bridge networks, whereas dT40 has the stronger
  ion-mediated component.
- LYS147 is the dominant shared interaction site. GLU111 and THR9 contribute
  more strongly to dA40 contacts, while ASN17 and K+-mediated contacts are more
  prominent for dT40.

### Pore-dependent behavior

- 7AHL produces more DNA–protein contacts, approximately three to four times
  more constriction contacts, 2.4–3.0 times more constriction water bridges and
  a larger SMD force magnitude than 3B07.
- The strongest dA40/dT40 contrast in 3B07 is the broader wall contact and K+
  bridging of dT40.
- The strongest dA40/dT40 contrast in 7AHL is the direct hydrogen-bond and water
  network of dA40 versus the ion-mediated interaction pattern of dT40.

### Experimental signal interpretation

- In 3B07, the broader and more ion-mediated dT40 contact pattern supports a
  stronger and more heterogeneous current-blockade response than dA40.
- In 7AHL, dA40 and dT40 have similar mechanical resistance but distinct local
  interaction chemistries, supporting differences in blockade amplitude,
  waveform and current fluctuation pattern.
- The more flexible dT40 chain can access more wall-bound geometries, whereas
  the more strongly stacked dA40 chain retains a more ordered translocation
  conformation.

## Trajectory movies

- [dA40 through 3B07](polyda40-3b07/result/trajectory.mp4)
- [dA40 through 7AHL](polyda40-7ahl/result/trajectory.mp4)
- [dT40 through 3B07](polydt40-3b07/result/trajectory.mp4)
- [dT40 through 7AHL](polydt40-7ahl/result/trajectory.mp4)
