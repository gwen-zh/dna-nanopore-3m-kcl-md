#!/usr/bin/env python3
"""Strict contiguous 5'-first SMD pathway for dA40/dT40 and 7AHL/3B07."""
import importlib.util
import json
import math
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
WORKSPACE = ROOT.parent

SPEC = importlib.util.spec_from_file_location(
    "accelerated", WORKSPACE / "accelerated-translocation-20261003/run.py"
)
accel = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(accel)

SEQUENCE = os.environ.get("SEQUENCE", "").lower()
PORE = os.environ.get("PORE", "").lower()
SYSTEM = Path(os.environ.get("SYSTEM_DIR", ""))
SOURCE = Path(os.environ.get("SOURCE_PREFIX", ""))
SOURCE_FRAME_NS = os.environ.get("SOURCE_FRAME_NS", "")
SOURCE_DONOR_SHA256 = os.environ.get("SOURCE_DONOR_SHA256", "")
FIELD = 0.1695
GROUP_SIZE = 3
PULL_RATE_A_PER_NS = 50.0
MIN_STAGE_STEPS = 250_000
MAX_STAGE_STEPS = 2_500_000
FRAME_QUANTUM = 25_000
TARGET_BELOW_TRANS_A = 10.0
FIRST_TARGET_BELOW_TRANS_A = 18.0
MAX_ATTEMPTS_PER_GROUP = 3


def require_configuration():
    if SEQUENCE not in {"da", "dt"}:
        raise ValueError("SEQUENCE must be da or dt")
    if PORE not in {"7ahl", "3b07"}:
        raise ValueError("PORE must be 7ahl or 3b07")
    for path in [SYSTEM / "system.psf", SYSTEM / "system.pdb",
                 SYSTEM / "system_settings.namd",
                 SYSTEM / "restraints_prod.pdb"]:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    for suffix in ["coor", "vel", "xsc"]:
        path = Path(str(SOURCE) + "." + suffix)
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)


def pore_from_system(system):
    name = str(system).lower()
    if PORE not in name:
        raise RuntimeError(f"Configured pore {PORE} is not present in {system}")
    return PORE


class PoreAwareAuditor:
    """Adapter around the validated dA and dT monitors with pore-aware rings."""
    def __init__(self, sequence, system, reference):
        self.sequence = sequence
        self.system = Path(system)
        self.reference = Path(reference)
        pore = pore_from_system(system)
        if sequence == "da":
            sys.path.insert(0, str(WORKSPACE / "da-entry-20260919"))
            from monitor import Monitor, binary
            self.monitor = Monitor(self.system / "system.psf", reference)
            self.binary = binary
            if pore == "3b07":
                atoms = self.monitor.atoms
                segid = np.asarray([a[1] for a in atoms])
                resid = np.asarray([int(a[2]) for a in atoms])
                resname = np.asarray([a[3] for a in atoms])
                name = np.asarray([a[4] for a in atoms])
                mass = np.asarray([float(a[7]) for a in atoms])
                protein = np.isin(resname, list(self.monitor.__class__.__module__ and {
                    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
                    "HIS", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "MET",
                    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                }))
                ring = protein & (
                    (np.isin(segid, list("BDFH")) & np.isin(resid, [108, 132])) |
                    (np.isin(segid, list("ACEG")) & np.isin(resid, [114, 140]))
                )
                self.monitor.constriction_ca = np.flatnonzero(ring & (name == "CA"))
                self.monitor.constriction_heavy = np.flatnonzero(ring & (mass > 2))
                reference_xyz = self.monitor.reference
                self.monitor.cis_ca = np.flatnonzero(
                    protein & (name == "CA") &
                    (reference_xyz[:, 2] >= 65) & (reference_xyz[:, 2] <= 75)
                )
                if (len(self.monitor.constriction_ca) < 8 or
                        len(self.monitor.cis_ca) < 8):
                    raise RuntimeError("Incomplete 3B07 landmarks for dA monitor")
            self.n = self.monitor.n
            self.first_c1 = int(self.monitor.c1[0])
            self.lead = np.asarray(self.monitor.lead, dtype=int)
            self.protein_ca = np.asarray(self.monitor.ca, dtype=int)
            self.constriction_ca = np.asarray(
                self.monitor.constriction_ca, dtype=int
            )
        else:
            sys.path.insert(0, str(WORKSPACE / "dt-pore-auto-20260915"))
            from validate import PoreChecker
            from checks import live_dna_frames
            self.checker = PoreChecker(self.system, pore)
            self.live_dna_frames = live_dna_frames
            self.n = self.checker.n
            self.first_c1 = int(self.checker.ids[int(self.checker.c1[0])])
            self.lead = np.asarray(
                self.checker.ids[self.checker.lead_local], dtype=int
            )
            self.protein_ca = np.asarray(self.checker.protein_ca, dtype=int)
            self.constriction_ca = np.asarray(
                self.checker.constriction_ca, dtype=int
            )

    def restart_metrics(self, prefix):
        prefix = Path(prefix)
        if self.sequence == "da":
            xsc = np.loadtxt(str(prefix) + ".xsc", comments="#")
            return self.monitor.measure(
                int(xsc[0]), self.binary(str(prefix) + ".coor", self.n),
                xsc[[1, 5, 9]],
            )
        self.checker.restart(prefix)
        return self.checker.entry_metrics(
            self.checker.binary(str(prefix) + ".coor"), self.checker.box
        )

    def frames(self, path, start, strict=False):
        if self.sequence == "da":
            yield from self.monitor.new_frames(path, start)
            return
        for index, step, raw, box in self.live_dna_frames(
            path, self.checker.n, int(self.checker.ids[-1]) + 1,
            start, strict,
        ):
            expected = (index + 1) * FRAME_QUANTUM
            if step != expected:
                raise ValueError(f"Unexpected dT DCD step {step}; expected {expected}")
            yield index, dict(
                step=step, time_ns=step * 1e-6,
                **self.checker.entry_metrics(raw, box),
            )


# run_stage constructs an Auditor internally. Replace it only in this process.
accel.base.Auditor = PoreAwareAuditor


def cell_z_A():
    text = (SYSTEM / "system_settings.namd").read_text()
    match = re.search(r"^set cell_z (\S+)", text, re.M)
    if not match:
        raise RuntimeError("cell_z missing from system_settings.namd")
    return float(match.group(1))


accel.base.voltage = lambda field: field * cell_z_A() / 23.060548867


def residue_group(auditor, first, last):
    if SEQUENCE == "da":
        resid = auditor.monitor.resid
        ids = auditor.monitor.dna[
            (resid[auditor.monitor.dna] >= first) &
            (resid[auditor.monitor.dna] <= last)
        ]
    else:
        atoms = auditor.checker.atoms
        mask = (
            (atoms.segid.to_numpy() == "AN1") &
            (atoms.resid.to_numpy() >= first) &
            (atoms.resid.to_numpy() <= last)
        )
        ids = np.flatnonzero(mask)
    if ids.size == 0:
        raise RuntimeError(f"No atoms found for DNA residues {first}-{last}")
    return np.asarray(ids, dtype=int)


def passage_order(auditor, prefix):
    prefix = Path(prefix)
    if SEQUENCE == "da":
        raw = auditor.binary(str(prefix) + ".coor", auditor.n)
        xsc = np.loadtxt(str(prefix) + ".xsc", comments="#")
        dna = auditor.monitor.whole_dna(raw, xsc[[1, 5, 9]])
        local = {int(g): i for i, g in enumerate(auditor.monitor.dna)}
        c1_z = np.asarray([
            dna[local[int(g)], 2] for g in auditor.monitor.c1
        ])
    else:
        from checks import join
        raw = auditor.checker.binary(str(prefix) + ".coor")
        dna = join(raw[auditor.checker.ids], auditor.checker.box,
                   auditor.checker.edges)
        c1_z = dna[auditor.checker.c1, 2]
    trans_z = float(raw[auditor.protein_ca, 2].min())
    passed = c1_z < trans_z
    prefix_count = 0
    for value in passed:
        if not value:
            break
        prefix_count += 1
    return dict(
        passed_prefix_count=prefix_count,
        total_C1_below_trans=int(passed.sum()),
        out_of_order_passage=bool(passed[prefix_count:].any()),
        complete_passage=bool(passed.all()),
        trans_reference_z_A=trans_z,
        C1_z_A=[float(v) for v in c1_z],
    )


def ordered_pull_config(auditor, moving, start, target, steps,
                        anchor=None, anchor_z=None):
    text = (
        "units real\ncolvarsTrajFrequency 1000\n"
        "colvarsRestartFrequency 25000\n"
        + accel.coordinate_block("front_z", moving, auditor.protein_ca,
                                 "distanceZ")
        + accel.coordinate_block("front_radius", moving, auditor.protein_ca,
                                 "distanceXY")
    )
    text += f"""harmonic {{
 name move_current_front_to_trans
 colvars front_z
 centers {start['relative_z_A']:.8f}
 targetCenters {target:.8f}
 targetNumSteps {steps}
 forceConstant 5.0
 outputCenters on
}}
harmonic {{
 name keep_current_front_on_axis
 colvars front_radius
 centers {start['radius_A']:.8f}
 targetCenters 0.0
 targetNumSteps {steps}
 forceConstant 0.5
 outputCenters on
}}
"""
    if anchor is not None:
        text += accel.coordinate_block(
            "five_prime_anchor_z", anchor, auditor.protein_ca, "distanceZ"
        )
        text += f"""harmonic {{
 name keep_five_prime_in_trans_chamber
 colvars five_prime_anchor_z
 centers {anchor_z:.8f}
 forceConstant 2.0
}}
"""
    return text


def stage_steps(displacement_A):
    raw = abs(displacement_A) / PULL_RATE_A_PER_NS * 1_000_000
    quantized = int(math.ceil(raw / FRAME_QUANTUM) * FRAME_QUANTUM)
    return min(MAX_STAGE_STEPS, max(MIN_STAGE_STEPS, quantized))


def validate_start(config):
    require_configuration()
    auditor = PoreAwareAuditor(SEQUENCE, SYSTEM, SOURCE)
    metrics = auditor.restart_metrics(SOURCE)
    accel.base.metrics_safe(metrics)
    order = passage_order(auditor, SOURCE)
    if order["out_of_order_passage"]:
        raise RuntimeError("Starting checkpoint already has out-of-order passage")
    return metrics, order


def run(models):
    config = {
        "system": SYSTEM,
        "seed": 202610080 + (0 if SEQUENCE == "da" else 10) +
                (0 if PORE == "7ahl" else 1),
    }
    start_metrics, start_order = validate_start(config)
    job_id = os.environ["SLURM_JOB_ID"]
    work = ROOT / f"{SEQUENCE}-{PORE}-run-{job_id}"
    work.mkdir(exist_ok=False)
    snapshot = work / "code-snapshot"
    snapshot.mkdir()
    for path in [Path(__file__), ROOT / "ordered_translocation.sbatch",
                 WORKSPACE / "accelerated-translocation-20261003/run.py",
                 WORKSPACE / "accelerated-translocation-20261003/stage.namd"]:
        shutil.copyfile(path, snapshot / path.name)

    plan = dict(
        sequence=f"poly(d{SEQUENCE[-1].upper()})40", pore=PORE.upper(),
        KCl_M=3, source=str(SOURCE.relative_to(WORKSPACE)),
        source_frame_ns=SOURCE_FRAME_NS,
        source_donor_sha256=SOURCE_DONOR_SHA256,
        system=str(SYSTEM.relative_to(WORKSPACE)),
        source_metrics=start_metrics, source_order=start_order,
        direction="5-prime first, cis to trans",
        strategy=("advance only the next contiguous three-residue front; "
                  "the 3-prime tail remains unbiased until it is the front"),
        field_z=FIELD, nominal_voltage_V=accel.base.voltage(FIELD),
        pull_rate_A_per_ns=PULL_RATE_A_PER_NS,
        group_size_residues=GROUP_SIZE,
        complete_passage_definition=(
            "all 40 C1 atoms below the trans reference with no out-of-order crossing"
        ),
        interpretation=("nonequilibrium ordered SMD pathway; not a voltage-only "
                        "translocation rate trajectory"),
    )
    accel.base.write_json(work / "run_plan.json", plan)
    accel.base.event(work, "job_start", job_id=job_id,
                     node=os.environ.get("SLURMD_NODENAME"),
                     models=models, plan=plan)

    current = SOURCE
    records = []
    while True:
        auditor = PoreAwareAuditor(SEQUENCE, SYSTEM, current)
        order = passage_order(auditor, current)
        if order["out_of_order_passage"]:
            raise RuntimeError("Detected non-5'-first passage; aborting")
        if order["complete_passage"]:
            break
        first = order["passed_prefix_count"] + 1
        last = min(40, first + GROUP_SIZE - 1)
        moving = residue_group(auditor, first, last)
        anchor = residue_group(auditor, 1, 3) if first > 3 else None
        succeeded = False
        for attempt in range(1, MAX_ATTEMPTS_PER_GROUP + 1):
            auditor = PoreAwareAuditor(SEQUENCE, SYSTEM, current)
            geometry = accel.group_geometry(auditor, current, moving)
            depth = (FIRST_TARGET_BELOW_TRANS_A if first == 1 else
                     TARGET_BELOW_TRANS_A + 5.0 * (attempt - 1))
            target = geometry["trans_relative_z_A"] - depth
            steps = stage_steps(geometry["relative_z_A"] - target)
            anchor_z = None if anchor is None else accel.group_geometry(
                auditor, current, anchor
            )["relative_z_A"]
            label = f"pull_r{first:02d}_{last:02d}_a{attempt}"
            colvars = work / f"{label}.colvars"
            colvars.write_text(ordered_pull_config(
                auditor, moving, geometry, target, steps, anchor, anchor_z
            ))
            accel.base.event(
                work, "ordered_front_plan", stage=label,
                residue_first=first, residue_last=last,
                start_geometry=geometry, target_relative_z_A=target,
                steps=steps, duration_ns=steps * 1e-6,
            )
            current, result = accel.base.run_stage(
                work, SEQUENCE, config, label, current, steps, FIELD, colvars
            )
            endpoint = PoreAwareAuditor(SEQUENCE, SYSTEM, current)
            endpoint_order = passage_order(endpoint, current)
            if endpoint_order["out_of_order_passage"]:
                raise RuntimeError(f"{label}: noncontiguous passage detected")
            record = dict(stage=label, residue_first=first,
                          residue_last=last, result=result,
                          passage_order=endpoint_order)
            records.append(record)
            accel.base.write_json(work / f"{label}.order.json", endpoint_order)
            if endpoint_order["passed_prefix_count"] >= last:
                succeeded = True
                break
        if not succeeded:
            raise RuntimeError(f"Residues {first}-{last} failed after three attempts")

    final_auditor = PoreAwareAuditor(SEQUENCE, SYSTEM, current)
    final_metrics = final_auditor.restart_metrics(current)
    final_order = passage_order(final_auditor, current)
    if not final_order["complete_passage"] or final_order["out_of_order_passage"]:
        raise RuntimeError("Strict 40/40 ordered completion gate failed")
    summary = dict(**plan, complete_passage_observed=True,
                   final_prefix=str(current.relative_to(WORKSPACE)),
                   final_metrics=final_metrics, final_order=final_order,
                   ordered_stages=records)
    accel.base.write_json(work / "workflow_validation.json", summary)
    accel.base.event(work, "workflow_complete", **summary)


def main():
    config = {"system": SYSTEM}
    metrics, order = validate_start(config)
    if os.environ.get("CHECK_ONLY") == "yes":
        print(json.dumps({"sequence": SEQUENCE, "pore": PORE,
                          "system": str(SYSTEM), "source": str(SOURCE),
                          "metrics": metrics, "passage_order": order}, indent=2))
        return
    run(accel.base.validate_allocation())


if __name__ == "__main__":
    main()
