"""DNA-only input/restart QA and read-only auditing of a growing NAMD DCD."""
import itertools
import json
import struct
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2] / 'nanopore' / 'analysis'))
from analyze import read_psf, traversal, join, minimum_image

IMAGE_STOP_A = 24.0  # 12 A cutoff plus 12 A precautionary margin, not a force.

def nearest_nonzero_image(x, box):
    x, box = np.asarray(x), np.asarray(box)
    if not np.isfinite(x).all() or not np.isfinite(box).all() or np.any(box <= 0):
        raise ValueError('Non-finite coordinates or invalid box')
    if np.any(np.ptp(x, axis=0) >= box):
        raise ValueError('DNA spans an entire cell dimension; stop and reassess box')
    tree = cKDTree(x)
    return min(float(tree.query(x + np.array(shift) * box)[0].min())
               for shift in itertools.product([-1, 0, 1], repeat=3) if shift != (0, 0, 0))

def record(stream):
    size = stream.read(4)
    if len(size) != 4:
        raise EOFError('Incomplete DCD record')
    n = struct.unpack('<i', size)[0]
    if n < 0 or n > 100000000:
        raise ValueError('Invalid DCD record size')
    data = stream.read(n)
    if len(data) != n or stream.read(4) != size:
        raise EOFError('Incomplete DCD record or mismatched footer')
    return data

def live_dna_frames(path, natoms, ndna, start=0, strict=False):
    """Yield complete frames only; do not trust growing-file NSET header.

    DNA is explicitly validated to be the first ndna atoms. Coordinate-record
    sizes and footers are checked while seeking past the solvent, avoiding
    repeatedly reading megabytes of solvent to monitor a 1279-atom DNA chain.
    """
    path = Path(path)
    if not path.exists():
        if strict:
            raise ValueError('Missing stage DCD')
        return
    with path.open('rb') as stream:
        header = record(stream)
        if len(header) != 84 or header[:4] != b'CORD':
            raise ValueError('Unsupported DCD header')
        control = np.frombuffer(header[4:], dtype='<i4')
        first, stride = map(int, control[1:3])
        if control[8] != 0 or control[10] != 1 or stride <= 0:
            raise ValueError('Fixed atoms or missing DCD cell unsupported')
        record(stream)
        count = struct.unpack('<i', record(stream))[0]
        if count != natoms or not (0 < ndna <= count):
            raise ValueError('DCD/PSF atom count mismatch')
        offset = stream.tell()
        framebytes = 56 + 3 * (8 + 4 * count)
        payload = path.stat().st_size - offset
        available, trailing = divmod(payload, framebytes)
        if strict and trailing:
            raise ValueError('Completed stage has a partial trailing DCD frame')
        for index in range(start, available):
            stream.seek(offset + index * framebytes)
            cell = np.frombuffer(record(stream), dtype='<f8')
            if len(cell) != 6 or not np.all((abs(cell[[1, 3, 4]]) < 1e-6) |
                                          (abs(cell[[1, 3, 4]] - 90) < 1e-6)):
                raise ValueError('Only orthorhombic boxes supported')
            axes = []
            for _ in range(3):
                size = stream.read(4)
                if len(size) != 4 or struct.unpack('<i', size)[0] != 4 * count:
                    raise EOFError('Incomplete coordinate record')
                xyz = stream.read(4 * ndna)
                if len(xyz) != 4 * ndna:
                    raise EOFError('Incomplete DNA coordinates')
                stream.seek(4 * (count - ndna), 1)
                if stream.read(4) != size:
                    raise EOFError('Coordinate record is still being written')
                axes.append(np.frombuffer(xyz, dtype='<f4').astype(float))
            yield index, first + index * stride, np.column_stack(axes), cell[[0, 2, 5]].copy()

class Checker:
    def __init__(self, system):
        self.system = Path(system)
        self.atoms, self.bonds = read_psf(self.system / 'system.psf')
        self.n = len(self.atoms)
        ids = self.atoms.index[self.atoms.segid == 'DNA'].to_numpy()
        if not np.array_equal(ids, np.arange(1279)):
            raise ValueError('DNA atom order is not the expected leading 1279 atoms')
        dna, db = read_psf(self.system / 'dna.psf')
        columns = ['segid', 'resid', 'resname', 'name', 'type', 'charge', 'mass']
        if not dna[columns].equals(self.atoms.iloc[:1279][columns]):
            raise ValueError('DNA topology, charges or atom names changed')
        self.dbonds = self.bonds[(self.bonds < 1279).all(1)]
        if not np.array_equal(self.dbonds, db):
            raise ValueError('DNA covalent topology changed')
        self.edges = traversal(1279, self.dbonds)
        self.heavy = dna.mass.to_numpy() > 2
        self.mass = dna.mass.to_numpy()[self.heavy]
        self.c1 = dna.index[dna.name == "C1'"].to_numpy()
        if len(self.c1) != 40 or set(dna.resname) != {'THY'}:
            raise ValueError('Expected poly(dT)40')
        self.nk = int((self.atoms.name == 'POT').sum())
        self.ncl = int((self.atoms.name == 'CLA').sum())
        self.nwater = int(((self.atoms.resname == 'TIP3') & (self.atoms.name == 'OH2')).sum())
        if self.nk - self.ncl != 39 or abs(self.atoms.charge.sum()) > 1e-4:
            raise ValueError('Ion count/neutrality mismatch')

    def dna_metrics(self, xyz, box):
        whole = join(np.asarray(xyz[:1279]), box, self.edges)
        length = np.linalg.norm(whole[self.dbonds[:, 0]] - whole[self.dbonds[:, 1]], axis=1)
        if not np.isfinite(whole).all() or length.min() <= .5 or length.max() >= 2.2:
            raise ValueError('Invalid complete DNA bonded geometry')
        x = whole[self.heavy]
        gap = nearest_nonzero_image(x, box)
        center = np.average(x, axis=0, weights=self.mass)
        result = {
            'rg_A': float(np.sqrt(np.average(np.sum((x - center)**2, axis=1), weights=self.mass))),
            'end_to_end_A': float(np.linalg.norm(whole[self.c1[-1]] - whole[self.c1[0]])),
            'dna_max_bond_A': float(length.max()), 'nearest_periodic_image_A': gap,
            'box_x_A': float(box[0]), 'box_y_A': float(box[1]), 'box_z_A': float(box[2]),
            'K_box_molar': self.nk / (float(np.prod(box)) * 6.02214076e-4),
            'Cl_box_molar': self.ncl / (float(np.prod(box)) * 6.02214076e-4),
        }
        if gap <= IMAGE_STOP_A:
            raise ValueError(f'Precautionary periodic-image stop: {json.dumps(result)}')
        return result

    def binary(self, path):
        with Path(path).open('rb') as stream:
            if struct.unpack('<i', stream.read(4))[0] != self.n:
                raise ValueError('Binary atom count mismatch')
            xyz = np.fromfile(stream, dtype='<f8')
        if xyz.size != 3 * self.n or not np.isfinite(xyz).all():
            raise ValueError('Incomplete or non-finite binary data')
        return xyz.reshape(self.n, 3)

    def validate_coords(self, xyz, box):
        delta = minimum_image(xyz[self.bonds[:, 0]] - xyz[self.bonds[:, 1]], box)
        length = np.linalg.norm(delta, axis=1)
        if not np.isfinite(xyz).all() or length.min() <= .5 or length.max() >= 2.5:
            raise ValueError('Invalid whole-system covalent geometry')
        result = self.dna_metrics(xyz, box)
        result.update(atoms=self.n, dna_atoms=1279, n_K=self.nk, n_Cl=self.ncl,
                      n_waters=self.nwater, all_min_bond_A=float(length.min()),
                      all_max_bond_A=float(length.max()))
        return result

    def initial(self):
        xyz = np.array([[float(s[30:38]), float(s[38:46]), float(s[46:54])]
                        for s in (self.system / 'system.pdb').open()
                        if s.startswith(('ATOM  ', 'HETATM'))])
        if xyz.shape != (self.n, 3):
            raise ValueError('PDB atom count mismatch')
        donor = np.load(self.system / 'donor_coordinates_A.npy')
        if np.max(abs(xyz[:1279] - donor)) > .0011:
            raise ValueError('Solvation changed DNA conformation')
        box = np.repeat(180., 3)
        result = self.validate_coords(xyz, box)
        heavy = xyz[self.atoms.mass.to_numpy() > 2] % box
        closest = float(cKDTree(heavy, boxsize=box).query(heavy, k=2)[0][:, 1].min())
        if closest <= .75:
            raise ValueError(f'Near-coincident periodic heavy atoms: {closest}')
        result['minimum_any_heavy_atom_pair_A'] = closest
        return result

    def restart(self, prefix):
        prefix = str(prefix)
        xyz = self.binary(prefix + '.coor')
        self.binary(prefix + '.vel')
        xsc = np.loadtxt(prefix + '.xsc', comments='#')
        if np.any(abs(xsc[[2, 3, 4, 6, 7, 8]]) > 1e-6):
            raise ValueError('Nonorthorhombic XSC')
        result = self.validate_coords(xyz, xsc[[1, 5, 9]])
        result['stage_step'] = int(xsc[0])
        return result
