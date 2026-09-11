#!/usr/bin/env python3
"""Reject broken coordinates before spending GPU time; orthorhombic cells."""
import re
from pathlib import Path
import numpy as np


def main():
    settings = Path('system_settings.namd').read_text()
    box = np.array([float(re.search(rf'^set cell_{a} (\S+)', settings, re.M)[1])
                    for a in 'xyz'])
    with Path('system.psf').open() as stream:
        for line in stream:
            if '!NATOM' in line:
                natoms = int(line.split()[0])
                break
        atoms = [next(stream).split() for _ in range(natoms)]
        for line in stream:
            if '!NBOND' in line:
                nbonds = int(line.split()[0])
                break
        values = []
        while len(values) < 2 * nbonds:
            values.extend(map(int, next(stream).split()))
    bonds = np.array(values).reshape(-1, 2) - 1
    xyz = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])]
                    for line in Path('system.pdb').open()
                    if line.startswith(('ATOM  ', 'HETATM'))])
    if xyz.shape != (natoms, 3) or not np.isfinite(xyz).all():
        raise SystemExit('INVALID coordinate count or non-finite coordinates')
    delta = xyz[bonds[:, 0]] - xyz[bonds[:, 1]]
    delta -= box * np.rint(delta / box)
    length = np.linalg.norm(delta, axis=1)
    dna = np.array([a[1] == 'AN1' for a in atoms])
    dna_bonds = dna[bonds].all(axis=1)
    # The inherited membrane has strained bonds up to 3.32 A, to be minimized.
    # DNA must already be chemically intact before minimization.
    bad = (length < 0.5) | (length > 3.5) | (dna_bonds & (length > 2.2))
    if dna.sum() != 1279 or dna_bonds.sum() < 1278:
        raise SystemExit('INVALID DNA topology/atom count')
    if bad.any():
        for index in np.flatnonzero(bad)[:20]:
            i, j = bonds[index]
            print('INVALID_BOND', atoms[i][:6], atoms[j][:6], length[index])
        raise SystemExit(f'INVALID {bad.sum()} bonds')
    print(f'BONDS_VALIDATED atoms={natoms} DNA_atoms={dna.sum()} '
          f'DNA_max_bond_A={length[dna_bonds].max():.6f} '
          f'all_max_bond_A={length.max():.6f}')


if __name__ == '__main__':
    main()
