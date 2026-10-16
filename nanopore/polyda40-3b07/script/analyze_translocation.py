#!/usr/bin/env python3
"""Analyze a complete segmented 5'-first nanopore translocation trajectory."""

import argparse
import collections
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from trajectory import DCD, binary, edges, psf, whole, write_dcd, write_pdb, write_psf


AMINO = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}
BACKBONE = {
    "P", "O1P", "O2P", "O5'", "C5'", "H5'", "H5''", "C4'", "H4'",
    "O4'", "C1'", "H1'", "C2'", "H2'", "H2''", "C3'", "H3'", "O3'",
}


def rotation(mobile, target):
    center = mobile.mean(0)
    reference = target.mean(0)
    u, _, vt = np.linalg.svd((mobile - center).T @ (target - reference))
    matrix = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
    return center, matrix, reference


def minimum_image(delta, box):
    return delta - box * np.rint(delta / box)


def atom_element(name, mass):
    if mass < 2:
        return "H"
    text = str(name).lstrip("0123456789")
    return text[:1].upper()


def connected_tree(global_ids, bonds, total_atoms):
    inverse = np.full(total_atoms, -1, dtype=int)
    inverse[global_ids] = np.arange(len(global_ids))
    local = inverse[bonds[np.isin(bonds, global_ids).all(1)]]
    return local, edges(len(global_ids), local)


def read_source(workspace, plan, atom_count):
    prefix = workspace / plan["source"]
    xyz = binary(str(prefix) + ".coor", atom_count)
    box = np.loadtxt(str(prefix) + ".xsc", comments="#")[[1, 5, 9]]
    return prefix, xyz, box


def constriction_mask(atoms, pore, protein):
    segid = atoms.segid.to_numpy()
    resid = atoms.resid.to_numpy()
    if pore == "7ahl":
        return protein & np.isin(resid, [111, 113, 147])
    return protein & (
        (np.isin(segid, list("BDFH")) & np.isin(resid, [108, 132]))
        | (np.isin(segid, list("ACEG")) & np.isin(resid, [114, 140]))
    )


def donor_hydrogen_pairs(group_mask, atoms, bonds):
    mass = atoms.mass.to_numpy()
    name = atoms.name.to_numpy()
    allowed = np.array([
        atom_element(n, m) in {"N", "O", "S"} for n, m in zip(name, mass)
    ])
    result = []
    for i, j in bonds:
        if mass[i] < 2 <= mass[j] and group_mask[j] and allowed[j]:
            result.append((int(j), int(i)))
        elif mass[j] < 2 <= mass[i] and group_mask[i] and allowed[i]:
            result.append((int(i), int(j)))
    return result


def parse_colvars(run, stages):
    rows = []
    offset_ns = 0.0
    conversion_pn = 69.479
    for stage in stages:
        label = stage["stage"]
        path = run / f"{label}.colvars.traj"
        config = run / f"{label}.colvars"
        if not path.is_file() or not config.is_file():
            continue
        force_constant = 5.0
        lines = config.read_text().splitlines()
        in_move = False
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("name move_current_front_to_trans"):
                in_move = True
            elif in_move and stripped.startswith("forceConstant"):
                force_constant = float(stripped.split()[1])
                break
        header = None
        stage_last_step = 0
        for line in path.read_text().splitlines():
            if line.startswith("# step"):
                header = line[1:].split()
                continue
            if not line.strip() or line.startswith("#") or header is None:
                continue
            values = list(map(float, line.split()))
            record = dict(zip(header, values))
            stage_last_step = max(stage_last_step, int(record["step"]))
            displacement = record.get("front_z", np.nan) - record.get("x0_front_z", np.nan)
            rows.append({
                "stage": label,
                "stage_step": int(record["step"]),
                "path_time_ns": offset_ns + record["step"] * 1e-6,
                "front_z_A": record.get("front_z"),
                "target_z_A": record.get("x0_front_z"),
                "front_radius_A": record.get("front_radius"),
                "harmonic_force_kcal_mol_A": -force_constant * displacement,
                "harmonic_force_pN": -force_constant * displacement * conversion_pn,
            })
        if stage_last_step <= 0:
            raise RuntimeError(f"{label}: colvars trajectory contains no data")
        offset_ns += stage_last_step * 1e-6
    return pd.DataFrame(rows)


def run_analysis(args):
    workspace = args.workspace.resolve()
    run = args.run.resolve()
    system = args.system.resolve()
    output = args.out.resolve()
    output.mkdir(parents=True, exist_ok=False)

    plan = json.loads((run / "run_plan.json").read_text())
    workflow_path = run / "workflow_validation.json"
    if not workflow_path.is_file():
        raise RuntimeError("The trajectory is not complete: workflow_validation.json is absent")
    workflow = json.loads(workflow_path.read_text())
    final_order = workflow["final_order"]
    if not (
        workflow.get("complete_passage_observed")
        and final_order.get("complete_passage")
        and final_order.get("passed_prefix_count") == 40
        and not final_order.get("out_of_order_passage")
    ):
        raise RuntimeError("Strict ordered 40/40 completion gate failed")
    stages = workflow["ordered_stages"]

    atoms, bonds = psf(system / "system.psf")
    n_atoms = len(atoms)
    segid = atoms.segid.to_numpy()
    resid = atoms.resid.to_numpy()
    resname = atoms.resname.to_numpy()
    name = atoms.name.to_numpy()
    mass = atoms.mass.to_numpy()
    charge = atoms.charge.to_numpy()
    element = np.array([atom_element(n, m) for n, m in zip(name, mass)])

    dna_mask = segid == "AN1"
    protein_mask = np.isin(resname, list(AMINO))
    dna = np.flatnonzero(dna_mask)
    protein = np.flatnonzero(protein_mask)
    protein_heavy = np.flatnonzero(protein_mask & (mass > 2))
    protein_ca = np.flatnonzero(protein_mask & (name == "CA"))
    dna_heavy = np.flatnonzero(dna_mask & (mass > 2))
    dna_polar = np.flatnonzero(dna_mask & np.isin(element, ["N", "O", "S"]) & (mass > 2))
    dna_oxygen = np.flatnonzero(dna_mask & (element == "O") & (mass > 2))
    constriction = constriction_mask(atoms, args.pore, protein_mask)
    constriction_heavy = np.flatnonzero(constriction & (mass > 2))
    constriction_polar = np.flatnonzero(constriction & np.isin(element, ["N", "O", "S"]) & (mass > 2))
    acidic_oxygen = np.flatnonzero(
        protein_mask & np.isin(resname, ["ASP", "GLU"]) & (element == "O") & (mass > 2)
    )
    potassium = np.flatnonzero((name == "POT") | (resname == "POT"))
    water_oxygen = np.flatnonzero(
        np.isin(resname, ["TIP3", "HOH", "WAT", "SOL"]) & np.isin(name, ["OH2", "O"])
    )
    lipid_p = np.flatnonzero((resname == "POPC") & (name == "P"))
    c1 = np.flatnonzero(dna_mask & (name == "C1'"))
    if len(c1) != 40:
        raise ValueError(f"Expected 40 C1' atoms, found {len(c1)}")

    dna_bonds, dna_tree = connected_tree(dna, bonds, n_atoms)
    dna_inverse = np.full(n_atoms, -1, dtype=int)
    dna_inverse[dna] = np.arange(len(dna))
    chain_trees = []
    for chain in atoms.iloc[protein].segid.unique():
        ids = np.flatnonzero(protein_mask & (segid == chain))
        local_bonds, tree = connected_tree(ids, bonds, n_atoms)
        chain_trees.append((ids, local_bonds, tree))

    source_prefix, reference, reference_box = read_source(workspace, plan, n_atoms)
    for ids, _, tree in chain_trees:
        reference[ids] = whole(reference[ids], reference_box, tree)
    ca_reference = reference[protein_ca].copy()
    reference_center = ca_reference.mean(0)
    pore_xy = reference_center[:2]

    protein_acceptor = protein_mask & np.isin(element, ["N", "O", "S"]) & (charge <= -0.10) & (mass > 2)
    dna_acceptor = dna_mask & np.isin(element, ["N", "O", "S"]) & (charge <= -0.10) & (mass > 2)
    protein_acceptors = np.flatnonzero(protein_acceptor)
    dna_acceptors = np.flatnonzero(dna_acceptor)
    protein_donors = donor_hydrogen_pairs(protein_mask, atoms, bonds)
    dna_donors = donor_hydrogen_pairs(dna_mask, atoms, bonds)

    base_atoms = {}
    normal_atoms = {}
    for nt in range(1, 41):
        ids = np.flatnonzero(dna_mask & (resid == nt) & ~np.isin(name, list(BACKBONE)) & (mass > 2))
        base_atoms[nt] = ids
        wanted = ["N9", "C4", "C8"] if args.sequence == "da" else ["N1", "C2", "C6"]
        triplet = [np.flatnonzero(dna_mask & (resid == nt) & (name == atom_name)) for atom_name in wanted]
        if any(len(v) != 1 for v in triplet):
            raise ValueError(f"Missing ring-normal atoms for nucleotide {nt}")
        normal_atoms[nt] = np.array([int(v[0]) for v in triplet])

    selected = np.sort(np.r_[protein, dna, lipid_p])
    selected_inverse = np.full(n_atoms, -1, dtype=int)
    selected_inverse[selected] = np.arange(len(selected))
    visual_bonds = selected_inverse[bonds[np.isin(bonds, selected).all(1)]]

    protein_contact_frames = collections.Counter()
    nucleotide_contact_frames = collections.Counter()
    residue_pair_frames = collections.Counter()
    constriction_pair_frames = collections.Counter()
    hbond_frames = collections.Counter()
    potassium_bridge_frames = collections.Counter()
    base_angle_values = collections.defaultdict(list)
    base_contact_frames = collections.Counter()
    rows = []
    z_rows = []
    movie = []
    last_dna_center = reference[dna].mean(0)
    cumulative_ns = 0.0
    frame_index = 0
    dcd_header = None
    max_bond = 0.0
    minimum_gap = float("inf")

    def global_xyz(atom_id, raw, dna_whole):
        local = dna_inverse[atom_id]
        return dna_whole[local] if local >= 0 else raw[atom_id]

    def candidate_hbonds(donors, acceptors, raw, dna_whole, box, direction):
        if not donors or len(acceptors) == 0:
            return set()
        acceptor_xyz = np.mod(np.array([global_xyz(i, raw, dna_whole) for i in acceptors]), box)
        tree = cKDTree(acceptor_xyz, boxsize=box)
        found = set()
        for donor, hydrogen in donors:
            d = global_xyz(donor, raw, dna_whole)
            h = global_xyz(hydrogen, raw, dna_whole)
            near = tree.query_ball_point(np.mod(d, box), 3.5)
            for local_acceptor in near:
                acceptor = int(acceptors[local_acceptor])
                a = global_xyz(acceptor, raw, dna_whole)
                hd = minimum_image(d - h, box)
                ha = minimum_image(a - h, box)
                norm = np.linalg.norm(hd) * np.linalg.norm(ha)
                if norm == 0:
                    continue
                angle = math.degrees(math.acos(np.clip(np.dot(hd, ha) / norm, -1, 1)))
                if angle >= 135:
                    if direction == "dna_to_protein":
                        dna_atom, protein_atom = donor, acceptor
                    else:
                        dna_atom, protein_atom = acceptor, donor
                    found.add((
                        int(resid[dna_atom]), str(name[dna_atom]), str(segid[protein_atom]),
                        int(resid[protein_atom]), str(resname[protein_atom]), str(name[protein_atom]), direction,
                    ))
        return found

    def consume(stage, stage_time_ns, path_time_ns, raw, box):
        nonlocal last_dna_center, frame_index, max_bond, minimum_gap
        dna_xyz = whole(raw[dna], box, dna_tree)
        dna_xyz += box * np.rint((last_dna_center - dna_xyz.mean(0)) / box)
        last_dna_center = dna_xyz.mean(0)
        lengths = np.linalg.norm(dna_xyz[dna_bonds[:, 0]] - dna_xyz[dna_bonds[:, 1]], axis=1)
        if lengths.min() <= 0.5 or lengths.max() >= 2.2:
            raise ValueError(f"DNA covalent geometry failed in frame {frame_index}")
        max_bond = max(max_bond, float(lengths.max()))

        protein_xyz = raw.copy()
        for ids, _, tree in chain_trees:
            part = whole(raw[ids], box, tree)
            part += box * np.rint((reference[ids].mean(0) - part.mean(0)) / box)
            protein_xyz[ids] = part
        center, matrix, aligned_center = rotation(protein_xyz[protein_ca], ca_reference)
        dna_aligned = (dna_xyz - center) @ matrix + aligned_center
        protein_aligned = (protein_xyz - center) @ matrix + aligned_center

        p_tree = cKDTree(np.mod(raw[protein_heavy], box), boxsize=box)
        d_heavy_local = dna_inverse[dna_heavy]
        neighbours = p_tree.query_ball_point(np.mod(dna_xyz[d_heavy_local], box), 3.5)
        nearest = p_tree.query(np.mod(dna_xyz[d_heavy_local], box))[0]
        minimum_gap = min(minimum_gap, float(nearest.min()))
        frame_protein = set()
        frame_nt = set()
        frame_pairs = set()
        frame_constriction = set()
        frame_base_nt = set()
        contact_pair_count = 0
        constriction_set = set(map(int, constriction_heavy))
        for dna_atom, near in zip(dna_heavy, neighbours):
            nt = int(resid[dna_atom])
            region = "backbone" if name[dna_atom] in BACKBONE else "base"
            for local_protein in near:
                protein_atom = int(protein_heavy[local_protein])
                pkey = (str(segid[protein_atom]), int(resid[protein_atom]), str(resname[protein_atom]))
                frame_protein.add(pkey)
                frame_nt.add(nt)
                frame_pairs.add((nt, *pkey, region))
                if protein_atom in constriction_set:
                    frame_constriction.add((nt, *pkey, region))
                if region == "base":
                    frame_base_nt.add(nt)
                contact_pair_count += 1
        protein_contact_frames.update(frame_protein)
        nucleotide_contact_frames.update(frame_nt)
        residue_pair_frames.update(frame_pairs)
        constriction_pair_frames.update(frame_constriction)
        base_contact_frames.update(frame_base_nt)

        hbonds = candidate_hbonds(dna_donors, protein_acceptors, raw, dna_xyz, box, "dna_to_protein")
        hbonds |= candidate_hbonds(protein_donors, dna_acceptors, raw, dna_xyz, box, "protein_to_dna")
        hbond_frames.update(hbonds)

        dna_o_xyz = np.mod(dna_xyz[dna_inverse[dna_oxygen]], box)
        potassium_xyz = np.mod(raw[potassium], box)
        acid_xyz = np.mod(raw[acidic_oxygen], box)
        dna_o_tree = cKDTree(dna_o_xyz, boxsize=box)
        acid_tree = cKDTree(acid_xyz, boxsize=box)
        near_dna = dna_o_tree.query_ball_point(potassium_xyz, 3.5)
        near_acid = acid_tree.query_ball_point(potassium_xyz, 3.5)
        frame_k_bridges = set()
        for dlist, alist in zip(near_dna, near_acid):
            for di in dlist:
                dna_atom = int(dna_oxygen[di])
                for ai in alist:
                    protein_atom = int(acidic_oxygen[ai])
                    frame_k_bridges.add((
                        int(resid[dna_atom]), str(name[dna_atom]), str(segid[protein_atom]),
                        int(resid[protein_atom]), str(resname[protein_atom]), str(name[protein_atom]),
                    ))
        potassium_bridge_frames.update(frame_k_bridges)

        water_mod = np.mod(raw[water_oxygen], box)
        water_tree = cKDTree(water_mod, boxsize=box)
        dna_water = set()
        for values in water_tree.query_ball_point(np.mod(dna_xyz[dna_inverse[dna_polar]], box), 3.5):
            dna_water.update(values)
        constriction_water = set()
        for values in water_tree.query_ball_point(np.mod(raw[constriction_polar], box), 3.5):
            constriction_water.update(values)
        water_bridges = dna_water & constriction_water

        centers = []
        normals = []
        angles = []
        for nt in range(1, 41):
            base_local = dna_inverse[base_atoms[nt]]
            triplet = dna_aligned[dna_inverse[normal_atoms[nt]]]
            normal = np.cross(triplet[1] - triplet[0], triplet[2] - triplet[0])
            normal /= np.linalg.norm(normal)
            angle = math.degrees(math.acos(np.clip(abs(normal[2]), 0, 1)))
            centers.append(dna_aligned[base_local].mean(0))
            normals.append(normal)
            angles.append(angle)
            base_angle_values[nt].append(angle)
        centers = np.asarray(centers)
        normals = np.asarray(normals)
        delta = centers[1:] - centers[:-1]
        distance = np.linalg.norm(delta, axis=1)
        normal_alignment = np.abs(np.sum(normals[1:] * normals[:-1], axis=1))
        stacked = (distance <= 5.5) & (normal_alignment >= math.cos(math.pi / 4))

        c1_z = dna_aligned[dna_inverse[c1], 2]
        ca_z = protein_aligned[protein_ca, 2]
        constr_z = float(protein_aligned[np.flatnonzero(constriction & (name == "CA")), 2].mean())
        if args.pore == "7ahl":
            cis_ca = np.flatnonzero(protein_mask & (name == "CA") & np.isin(resid, [17, 19]))
        else:
            cis_ca = np.flatnonzero(protein_mask & (name == "CA") & (reference[:, 2] >= 65) & (reference[:, 2] <= 75))
        cis_z = float(protein_aligned[cis_ca, 2].mean())
        trans_z = float(ca_z.min())
        heavy_aligned = dna_aligned[dna_inverse[dna_heavy]]
        weights = mass[dna_heavy]
        com = np.average(heavy_aligned, axis=0, weights=weights)
        rg = float(np.sqrt(np.average(np.sum((heavy_aligned - com) ** 2, axis=1), weights=weights)))

        rows.append({
            "frame": frame_index,
            "stage": stage,
            "stage_time_ns": stage_time_ns,
            "path_time_ns": path_time_ns,
            "n_C1_below_cis": int((c1_z < cis_z).sum()),
            "n_C1_below_constriction": int((c1_z < constr_z).sum()),
            "n_C1_below_trans": int((c1_z < trans_z).sum()),
            "five_prime_C1_z_A": float(c1_z[0]),
            "three_prime_C1_z_A": float(c1_z[-1]),
            "cis_reference_z_A": cis_z,
            "constriction_reference_z_A": constr_z,
            "trans_reference_z_A": trans_z,
            "Rg_A": rg,
            "end_to_end_A": float(np.linalg.norm(dna_aligned[dna_inverse[c1[0]]] - dna_aligned[dna_inverse[c1[-1]]])),
            "DNA_protein_atom_contact_pairs_3p5A": contact_pair_count,
            "contacting_DNA_nucleotides_3p5A": len(frame_nt),
            "contacting_protein_residues_3p5A": len(frame_protein),
            "constriction_residue_pairs_3p5A": len(frame_constriction),
            "candidate_hydrogen_bonds": len(hbonds),
            "potassium_mediated_DNA_acidic_residue_bridges": len(frame_k_bridges),
            "DNA_hydration_waters_3p5A": len(dna_water),
            "DNA_constriction_water_bridges_3p5A": len(water_bridges),
            "adjacent_base_stack_count": int(stacked.sum()),
            "mean_base_normal_tilt_deg": float(np.mean(angles)),
            "minimum_DNA_protein_heavy_A": float(nearest.min()),
            "DNA_max_bond_A": float(lengths.max()),
        })
        z_rows.append({"frame": frame_index, "path_time_ns": path_time_ns, **{
            f"nt_{i}_C1_z_A": float(value) for i, value in enumerate(c1_z, 1)
        }})

        display = protein_xyz[selected].copy()
        lp = selected_inverse[lipid_p]
        display[lp] = raw[lipid_p] + box * np.rint((aligned_center - raw[lipid_p]) / box)
        display = (display - center) @ matrix + aligned_center
        display[selected_inverse[dna]] = dna_aligned
        movie.append(display.astype(np.float32))
        frame_index += 1
        if frame_index % 25 == 0:
            print(f"{args.sequence}-{args.pore}: processed {frame_index} frames", flush=True)

    consume("source", 0.0, 0.0, reference.copy(), reference_box)
    for stage_record in stages:
        label = stage_record["stage"]
        dcd = DCD(run / f"{label}.dcd", n_atoms, True)
        if dcd_header is None:
            dcd_header = dcd.header
        stage_duration_ns = 0.0
        for step, raw, box in dcd:
            stage_duration_ns = step * 1e-6
            consume(label, stage_duration_ns, cumulative_ns + stage_duration_ns, raw, box)
        if stage_duration_ns <= 0:
            raise RuntimeError(f"{label}: DCD contains no complete frames")
        cumulative_ns += stage_duration_ns

    table = pd.DataFrame(rows)
    table.to_csv(output / "timeseries.csv", index=False)
    pd.DataFrame(z_rows).to_csv(output / "nucleotide_z.csv", index=False)
    denominator = len(table)

    def counter_table(counter, columns):
        records = []
        for key, value in counter.items():
            values = (key,) if len(columns) == 1 else key
            records.append(dict(
                zip(columns, values),
                contact_frames=value,
                frame_fraction=value / denominator,
            ))
        return pd.DataFrame(records, columns=columns + ["contact_frames", "frame_fraction"])

    protein_table = counter_table(protein_contact_frames, ["chain", "protein_resid", "protein_resname"])
    if len(protein_table):
        protein_table = protein_table.sort_values(["frame_fraction", "chain", "protein_resid"], ascending=[False, True, True])
    protein_table.to_csv(output / "protein_contact_occupancy.csv", index=False)
    nt_table = counter_table(nucleotide_contact_frames, ["nucleotide"])
    if len(nt_table):
        nt_table = nt_table.sort_values("nucleotide")
    nt_table.to_csv(output / "nucleotide_contact_occupancy.csv", index=False)
    pair_table = counter_table(
        residue_pair_frames,
        ["nucleotide", "chain", "protein_resid", "protein_resname", "DNA_region"],
    )
    if len(pair_table):
        pair_table = pair_table.sort_values("frame_fraction", ascending=False)
    pair_table.to_csv(output / "DNA_protein_residue_pair_occupancy.csv", index=False)
    constriction_table = counter_table(
        constriction_pair_frames,
        ["nucleotide", "chain", "protein_resid", "protein_resname", "DNA_region"],
    )
    if len(constriction_table):
        constriction_table = constriction_table.sort_values("frame_fraction", ascending=False)
    constriction_table.to_csv(output / "constriction_contact_occupancy.csv", index=False)
    hbond_table = counter_table(
        hbond_frames,
        ["nucleotide", "DNA_atom", "chain", "protein_resid", "protein_resname", "protein_atom", "direction"],
    )
    if len(hbond_table):
        hbond_table = hbond_table.sort_values("frame_fraction", ascending=False)
    hbond_table.to_csv(output / "candidate_hydrogen_bond_occupancy.csv", index=False)
    bridge_table = counter_table(
        potassium_bridge_frames,
        ["nucleotide", "DNA_oxygen", "chain", "protein_resid", "protein_resname", "protein_oxygen"],
    )
    if len(bridge_table):
        bridge_table = bridge_table.sort_values("frame_fraction", ascending=False)
    bridge_table.to_csv(output / "potassium_bridge_occupancy.csv", index=False)

    base_rows = []
    for nt in range(1, 41):
        values = np.asarray(base_angle_values[nt])
        base_rows.append({
            "nucleotide": nt,
            "mean_normal_tilt_deg": float(values.mean()),
            "median_normal_tilt_deg": float(np.median(values)),
            "standard_deviation_deg": float(values.std()),
            "base_contact_frame_fraction": base_contact_frames[nt] / denominator,
        })
    base_table = pd.DataFrame(base_rows)
    base_table.to_csv(output / "base_orientation_and_contact.csv", index=False)

    force_table = parse_colvars(run, stages)
    force_table.to_csv(output / "pulling_force.csv", index=False)

    visual_atoms = atoms.iloc[selected].reset_index(drop=True)
    visual_xyz = np.asarray(movie, dtype=np.float32)
    write_psf(output / "visualization.psf", visual_atoms, visual_bonds)
    write_pdb(output / "start.pdb", visual_atoms, visual_xyz[0])
    write_pdb(output / "endpoint.pdb", visual_atoms, visual_xyz[-1])
    write_dcd(output / "visualization.dcd", visual_xyz, dcd_header)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
    axes[0].plot(table.path_time_ns, table.n_C1_below_cis, label="below cis")
    axes[0].plot(table.path_time_ns, table.n_C1_below_constriction, label="below constriction")
    axes[0].plot(table.path_time_ns, table.n_C1_below_trans, label="below trans")
    axes[0].set_ylabel("C1' count"); axes[0].legend(ncol=3)
    axes[1].plot(table.path_time_ns, table.DNA_protein_atom_contact_pairs_3p5A, label="DNA–protein atom contacts")
    axes[1].plot(table.path_time_ns, table.constriction_residue_pairs_3p5A, label="constriction residue pairs")
    axes[1].set_ylabel("contacts"); axes[1].legend()
    axes[2].plot(table.path_time_ns, table.candidate_hydrogen_bonds, label="candidate H bonds")
    axes[2].plot(table.path_time_ns, table.potassium_mediated_DNA_acidic_residue_bridges, label="K+ bridges")
    axes[2].plot(table.path_time_ns, table.DNA_constriction_water_bridges_3p5A, label="water bridges")
    axes[2].set_ylabel("interaction count"); axes[2].legend(ncol=3)
    axes[3].plot(table.path_time_ns, table.Rg_A, label="Rg")
    axes[3].plot(table.path_time_ns, table.end_to_end_A, label="end-to-end")
    axes[3].set_ylabel("Å"); axes[3].set_xlabel("cumulative simulated path time (ns)"); axes[3].legend()
    fig.suptitle(f"poly(d{args.sequence[-1].upper()})40 through {args.pore.upper()}: ordered nonequilibrium path")
    fig.tight_layout(); fig.savefig(output / "overview.png", dpi=180); plt.close(fig)

    top_keys = list(protein_contact_frames.most_common(24))
    top_residues = [key for key, _ in top_keys]
    matrix = np.zeros((40, len(top_residues)))
    for j, pkey in enumerate(top_residues):
        for nt in range(1, 41):
            matrix[nt - 1, j] = sum(
                value for key, value in residue_pair_frames.items()
                if key[0] == nt and tuple(key[1:4]) == pkey
            ) / denominator
    fig, ax = plt.subplots(figsize=(max(10, len(top_residues) * 0.45), 7))
    image = ax.imshow(matrix, aspect="auto", origin="lower", cmap="magma", vmin=0)
    ax.set_yticks(np.arange(0, 40, 4)); ax.set_yticklabels(np.arange(1, 41, 4))
    ax.set_xticks(np.arange(len(top_residues)))
    ax.set_xticklabels([f"{c}:{r}{n}" for c, r, n in top_residues], rotation=60, ha="right")
    ax.set_ylabel("DNA nucleotide"); ax.set_xlabel("most-contacted protein residues")
    fig.colorbar(image, ax=ax, label="trajectory frame fraction")
    fig.tight_layout(); fig.savefig(output / "contact_heatmap.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(base_table.nucleotide, base_table.mean_normal_tilt_deg, marker="o", ms=3)
    axes[0].fill_between(
        base_table.nucleotide,
        base_table.mean_normal_tilt_deg - base_table.standard_deviation_deg,
        base_table.mean_normal_tilt_deg + base_table.standard_deviation_deg,
        alpha=0.2,
    )
    axes[0].set_ylabel("base-normal tilt (deg)")
    axes[1].bar(base_table.nucleotide, base_table.base_contact_frame_fraction)
    axes[1].set_ylabel("base-contact fraction"); axes[1].set_xlabel("nucleotide")
    fig.tight_layout(); fig.savefig(output / "base_orientation.png", dpi=180); plt.close(fig)

    if len(force_table):
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(force_table.path_time_ns, force_table.harmonic_force_pN, lw=0.7)
        ax.axhline(0, color="black", lw=0.5)
        ax.set_xlabel("cumulative simulated path time (ns)")
        ax.set_ylabel("instantaneous SMD spring force (pN)")
        ax.set_title("Biased pulling-force trace; not an equilibrium free-energy profile")
        fig.tight_layout(); fig.savefig(output / "pulling_force.png", dpi=180); plt.close(fig)

    summary = {
        "sequence": plan["sequence"],
        "pore": plan["pore"],
        "KCl_M_nominal": plan["KCl_M"],
        "input_ssDNA_frame_ns": plan.get("source_frame_ns"),
        "field_z": plan["field_z"],
        "nominal_box_voltage_V": plan["nominal_voltage_V"],
        "pull_rate_A_per_ns": plan["pull_rate_A_per_ns"],
        "direction": plan["direction"],
        "complete_ordered_passage": True,
        "frames_including_source": denominator,
        "cumulative_simulated_path_ns": float(cumulative_ns),
        "first_complete_frame": int(table.index[table.n_C1_below_trans == 40][0]),
        "first_complete_path_time_ns": float(table.loc[table.n_C1_below_trans == 40, "path_time_ns"].iloc[0]),
        "minimum_DNA_protein_heavy_A": minimum_gap,
        "maximum_DNA_bond_A": max_bond,
        "mean_DNA_protein_atom_contacts_3p5A": float(table.DNA_protein_atom_contact_pairs_3p5A.mean()),
        "maximum_DNA_protein_atom_contacts_3p5A": int(table.DNA_protein_atom_contact_pairs_3p5A.max()),
        "mean_candidate_hydrogen_bonds": float(table.candidate_hydrogen_bonds.mean()),
        "mean_potassium_bridges": float(table.potassium_mediated_DNA_acidic_residue_bridges.mean()),
        "mean_DNA_constriction_water_bridges": float(table.DNA_constriction_water_bridges_3p5A.mean()),
        "visualization_atom_count": int(len(selected)),
        "interpretation": (
            "Single high-field ordered-SMD pathway. Contacts, candidate hydrogen bonds, ion/water bridges, "
            "base orientations and passage order are geometric nonequilibrium observables; they are not "
            "binding free energies, experimental -120 mV dwell times, or an unbiased translocation rate."
        ),
        "candidate_hydrogen_bond_definition": "donor-acceptor <=3.5 A and D-H-A >=135 degrees",
        "potassium_bridge_definition": "one K+ within 3.5 A of a DNA oxygen and an ASP/GLU oxygen",
        "water_bridge_definition": "one water oxygen within 3.5 A of DNA polar atoms and constriction polar atoms",
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    zmin = float(visual_xyz[:, :, 2].min())
    zmax = float(visual_xyz[:, :, 2].max())
    (output / "render.json").write_text(json.dumps({
        "pore": args.pore,
        "sequence": args.sequence,
        "center": [float(pore_xy[0]), float(pore_xy[1]), (zmin + zmax) / 2],
        "z_bounds": [zmin, zmax],
        "frames": denominator,
        "status": "complete_40_of_40",
    }, indent=2) + "\n")
    (output / "analysis_validation.json").write_text(json.dumps({
        "workflow_complete": True,
        "strict_passed_prefix_count": 40,
        "out_of_order_passage": False,
        "frames_processed": denominator,
        "endpoint_C1_below_trans": int(table.n_C1_below_trans.iloc[-1]),
        "minimum_DNA_protein_heavy_A": minimum_gap,
        "maximum_DNA_bond_A": max_bond,
    }, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--system", type=Path, required=True)
    parser.add_argument("--sequence", choices=["da", "dt"], required=True)
    parser.add_argument("--pore", choices=["7ahl", "3b07"], required=True)
    parser.add_argument("--out", type=Path, required=True)
    run_analysis(parser.parse_args())


if __name__ == "__main__":
    main()
