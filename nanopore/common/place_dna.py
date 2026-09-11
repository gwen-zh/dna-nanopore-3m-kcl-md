#!/usr/bin/env python3
import math
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}


def read_atoms(path):
    records = []
    for line in Path(path).read_text().splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        records.append(
            {
                "line": line,
                "name": line[12:16].strip(),
                "resname": line[17:21].strip(),
                "resid": int(line[22:26]),
                "element": (line[76:78].strip() or line[12:16].strip()[0]).upper(),
                "xyz": np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                ),
            }
        )
    return records


def rotation_between(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    cross = np.cross(a, b)
    cosine = float(np.clip(np.dot(a, b), -1.0, 1.0))
    sine = np.linalg.norm(cross)
    if sine < 1.0e-12:
        if cosine > 0:
            return np.eye(3)
        axis = np.cross(a, np.array([1.0, 0.0, 0.0]))
        if np.linalg.norm(axis) < 1.0e-12:
            axis = np.cross(a, np.array([0.0, 1.0, 0.0]))
        axis /= np.linalg.norm(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    kx, ky, kz = cross
    k = np.array([[0.0, -kz, ky], [kz, 0.0, -kx], [-ky, kx, 0.0]])
    return np.eye(3) + k + k @ k * ((1.0 - cosine) / (sine * sine))


def count_pairs(tree, xyz, cutoff):
    return sum(len(neighbors) for neighbors in tree.query_ball_point(xyz, cutoff))


def main():
    if len(sys.argv) != 13:
        raise SystemExit(
            "usage: pore.pdb donor.pdb output.pdb target_z minx miny minz "
            "maxx maxy maxz label margin"
        )
    pore_path, donor_path, output_path = sys.argv[1:4]
    target_z = float(sys.argv[4])
    bounds = np.array([float(v) for v in sys.argv[5:11]]).reshape(2, 3)
    label = sys.argv[11]
    margin = float(sys.argv[12])

    donor = read_atoms(donor_path)
    pore = read_atoms(pore_path)
    coords = np.array([record["xyz"] for record in donor])
    lead = [r["xyz"] for r in donor if r["resid"] == 1 and r["name"] == "C1'"]
    tail = [r["xyz"] for r in donor if r["resid"] == 40 and r["name"] == "C1'"]
    if len(lead) != 1 or len(tail) != 1:
        raise SystemExit("donor must contain one C1' atom in residues 1 and 40")

    centered = coords - lead[0]
    align = rotation_between(tail[0] - lead[0], np.array([0.0, 0.0, 1.0]))
    oriented = centered @ align.T
    protein_xyz = np.array(
        [r["xyz"] for r in pore if r["resname"] in AMINO_ACIDS and r["element"] != "H"]
    )
    tree = cKDTree(protein_xyz)
    heavy_mask = np.array([r["element"] != "H" for r in donor])

    best = None
    for dz in range(0, 13, 2):
        found_at_z = []
        for phi in range(0, 360, 15):
            angle = math.radians(phi)
            rz = np.array(
                [
                    [math.cos(angle), -math.sin(angle), 0.0],
                    [math.sin(angle), math.cos(angle), 0.0],
                    [0.0, 0.0, 1.0],
                ]
            )
            candidate = oriented @ rz.T + np.array([0.0, 0.0, target_z + dz])
            low = candidate.min(axis=0)
            high = candidate.max(axis=0)
            inside = bool(np.all(low >= bounds[0] + margin) and np.all(high <= bounds[1] - margin))
            heavy = candidate[heavy_mask]
            hard = count_pairs(tree, heavy, 1.2)
            soft = count_pairs(tree, heavy, 2.2)
            if hard == 0 and inside:
                clearance = float(np.min(np.concatenate((low - bounds[0], bounds[1] - high))))
                found_at_z.append((soft, -clearance, phi, candidate, low, high))
        if found_at_z:
            best = min(found_at_z, key=lambda item: (item[0], item[1], item[2]))
            best_z = target_z + dz
            break

    if best is None:
        raise SystemExit(f"no clash-free placement with {margin:g} A box clearance")

    soft, neg_clearance, phi, placed, low, high = best
    lines = []
    atom_index = 0
    for line in Path(donor_path).read_text().splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            x, y, z = placed[atom_index]
            newline = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            lines.append(newline)
            atom_index += 1
        else:
            lines.append(line)
    Path(output_path).write_text("".join(lines))
    print(
        f"PLACEMENT label={label} target_z={best_z:g} rotation_z={phi} "
        f"hard_contacts=0 soft_contacts={soft} clearance={-neg_clearance:.3f} "
        f"bounds={low.tolist()} {high.tolist()}"
    )


if __name__ == "__main__":
    main()
