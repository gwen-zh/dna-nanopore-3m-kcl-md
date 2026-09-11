#!/usr/bin/env python3
"""Choose distinct water sites with periodic solute/ion exclusion distances."""
import json
import re
import sys
import zlib
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree


def choose_sites(candidates, solute, box, count, from_distance, between, seed):
    xyz = np.mod(candidates, box)
    nearest = cKDTree(np.mod(solute, box), boxsize=box).query(xyz)[0]
    eligible = nearest >= from_distance
    tree = cKDTree(xyz, boxsize=box)
    selected = []
    for index in np.random.default_rng(seed).permutation(len(xyz)):
        if not eligible[index]:
            continue
        selected.append(index)
        eligible[tree.query_ball_point(xyz[index], between)] = False
        if len(selected) == count:
            break
    if len(selected) != count:
        raise ValueError(f'Could place only {len(selected)} of {count} ions')
    selected = np.asarray(selected)
    ion_xyz = xyz[selected]
    min_ion = cKDTree(ion_xyz, boxsize=box).query(ion_xyz, k=2)[0][:, 1].min()
    min_solute = nearest[selected].min()
    if min_ion < between - 1e-6 or min_solute < from_distance - 1e-6:
        raise ValueError('Ion distance validation failed')
    return selected, float(min_ion), float(min_solute)


def main():
    count, from_distance, between = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    label = sys.argv[4]
    settings = Path('system_settings.namd').read_text()
    box = np.array([float(re.search(rf'^set cell_{a} (\S+)', settings, re.M)[1])
                    for a in 'xyz'])
    candidates = np.loadtxt('ion_candidates.tsv')
    solute = np.loadtxt('ion_solute.tsv')
    seed = zlib.crc32(label.encode()) ^ 20260911
    sites, min_ion, min_solute = choose_sites(candidates[:, 1:], solute, box,
                                             count, from_distance, between, seed)
    Path('ion_placement_validation.json').write_text(json.dumps({
        'count': count, 'seed': seed, 'box_A': box.tolist(),
        'minimum_ion_ion_A': min_ion, 'minimum_ion_solute_A': min_solute,
        'requested_from_A': from_distance, 'requested_between_A': between,
        'periodic_distances': True,
    }, indent=2) + '\n')
    print(' '.join(str(int(i)) for i in candidates[sites, 0]))


if __name__ == '__main__':
    main()
