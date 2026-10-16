#!/usr/bin/env python3
"""Deep insertion followed by a voltage trial and accelerated tail SMD."""
import argparse
import csv
import importlib.util
import json
import os
import shutil
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
WORKSPACE = ROOT.parent
BASE_PATH = (
    WORKSPACE / "constriction-start-translocation-20260929/run.py"
)
SPEC = importlib.util.spec_from_file_location("constriction_base", BASE_PATH)
base = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(base)
base.ROOT = ROOT

FIELD = 0.1695
DEEP_LEAD_TARGET_A = -24.3846618711014
HOLD_STEPS = 250_000
TRIAL_STEPS = 500_000
REINSERT_STEPS = 500_000
CENTER_STEPS = 500_000
PULL_STEPS = 4_000_000
TAIL_TARGET_BELOW_TRANS_A = 12.0

RUNS = {
    "da": {
        "label": "poly(dA)40",
        "system": WORKSPACE / "polyda40-3m-kcl-7ahl/formal-50ns",
        "source": (
            WORKSPACE /
            "da-entry-20260919/run-57400896/deeper_beta"
        ),
        "guide_steps": 0,
        "seed": 202610031,
    },
    "dt": {
        "label": "poly(dT)40",
        "system": (
            WORKSPACE /
            "dt-pore-auto-20260915/7ahl/run-57334241/system"
        ),
        "source": (
            WORKSPACE /
            "dt-pore-auto-20260915/7ahl/run-57334241/entry15.restart"
        ),
        "guide_steps": 1_000_000,
        "seed": 202610032,
    },
}


def dna_group(auditor, which):
    if auditor.sequence == "da":
        return np.asarray(
            getattr(auditor.monitor, which), dtype=int
        )
    local = getattr(auditor.checker, f"{which}_local")
    return np.asarray(auditor.checker.ids[local], dtype=int)


def group_geometry(auditor, prefix, group):
    prefix = Path(prefix)
    if auditor.sequence == "da":
        raw = auditor.binary(str(prefix) + ".coor", auditor.n)
        xsc = np.loadtxt(str(prefix) + ".xsc", comments="#")
        dna = auditor.monitor.whole_dna(raw, xsc[[1, 5, 9]])
        local = {g: i for i, g in enumerate(auditor.monitor.dna)}
        ids = np.asarray([local[int(i)] for i in group], dtype=int)
        xyz = dna[ids]
        weights = auditor.monitor.mass[group]
        protein_weights = auditor.monitor.mass[auditor.protein_ca]
    else:
        from checks import join
        raw = auditor.checker.binary(str(prefix) + ".coor")
        dna = join(
            raw[auditor.checker.ids],
            auditor.checker.box,
            auditor.checker.edges,
        )
        ids = group - int(auditor.checker.ids[0])
        xyz = dna[ids]
        weights = auditor.checker.dna_mass[ids]
        protein_weights = auditor.checker.atoms.mass.to_numpy()[
            auditor.protein_ca
        ]
    center = np.average(xyz, axis=0, weights=weights)
    reference = np.average(
        raw[auditor.protein_ca], axis=0, weights=protein_weights
    )
    trans_z = float(raw[auditor.protein_ca, 2].min())
    delta = center - reference
    return dict(
        relative_z_A=float(delta[2]),
        radius_A=float(np.linalg.norm(delta[:2])),
        trans_relative_z_A=float(trans_z - reference[2]),
    )


def numbers(ids):
    return " ".join(str(int(i) + 1) for i in ids)


def coordinate_block(name, group, reference, component):
    return f"""colvar {{
 name {name}
 width 1.0
 {component} {{
  main {{ atomNumbers {{ {numbers(group)} }} }}
  ref {{ atomNumbers {{ {numbers(reference)} }} }}
  axis (0.0, 0.0, 1.0)
  forceNoPBC yes
 }}
}}
"""


def lead_config(auditor, initial, target, steps):
    lead = dna_group(auditor, "lead")
    header = (
        "units real\ncolvarsTrajFrequency 1000\n"
        "colvarsRestartFrequency 25000\n"
    )
    return header + coordinate_block(
        "lead_z", lead, auditor.protein_ca, "distanceZ"
    ) + coordinate_block(
        "lead_radius", lead, auditor.protein_ca, "distanceXY"
    ) + f"""harmonic {{
 name axial_lead
 colvars lead_z
 centers {initial:.8f}
 targetCenters {target:.8f}
 targetNumSteps {steps}
 forceConstant 5.0
 outputCenters on
}}
harmonic {{
 name radial_lead
 colvars lead_radius
 centers 0.0
 forceConstant 0.5
}}
"""


def lead_hold_config(auditor, center):
    lead = dna_group(auditor, "lead")
    header = (
        "units real\ncolvarsTrajFrequency 1000\n"
        "colvarsRestartFrequency 25000\n"
    )
    return header + coordinate_block(
        "lead_z", lead, auditor.protein_ca, "distanceZ"
    ) + coordinate_block(
        "lead_radius", lead, auditor.protein_ca, "distanceXY"
    ) + f"""harmonic {{
 name axial_lead_hold
 colvars lead_z
 centers {center:.8f}
 forceConstant 2.0
}}
harmonic {{
 name radial_lead_hold
 colvars lead_radius
 centers 0.0
 forceConstant 0.5
}}
"""


def tail_center_config(auditor, lead_z, tail):
    lead = dna_group(auditor, "lead")
    tail_ids = dna_group(auditor, "tail")
    header = (
        "units real\ncolvarsTrajFrequency 1000\n"
        "colvarsRestartFrequency 25000\n"
    )
    return header + coordinate_block(
        "lead_z", lead, auditor.protein_ca, "distanceZ"
    ) + coordinate_block(
        "tail_z", tail_ids, auditor.protein_ca, "distanceZ"
    ) + coordinate_block(
        "tail_radius", tail_ids, auditor.protein_ca, "distanceXY"
    ) + f"""harmonic {{
 name hold_threaded_lead
 colvars lead_z
 centers {lead_z:.8f}
 forceConstant 2.0
}}
harmonic {{
 name hold_tail_height
 colvars tail_z
 centers {tail['relative_z_A']:.8f}
 forceConstant 2.0
}}
harmonic {{
 name center_tail
 colvars tail_radius
 centers {tail['radius_A']:.8f}
 targetCenters 0.0
 targetNumSteps {CENTER_STEPS}
 forceConstant 0.5
 outputCenters on
}}
"""


def tail_pull_config(auditor, tail, target):
    tail_ids = dna_group(auditor, "tail")
    header = (
        "units real\ncolvarsTrajFrequency 1000\n"
        "colvarsRestartFrequency 25000\n"
    )
    return header + coordinate_block(
        "tail_z", tail_ids, auditor.protein_ca, "distanceZ"
    ) + coordinate_block(
        "tail_radius", tail_ids, auditor.protein_ca, "distanceXY"
    ) + f"""harmonic {{
 name pull_tail_through_pore
 colvars tail_z
 centers {tail['relative_z_A']:.8f}
 targetCenters {target:.8f}
 targetNumSteps {PULL_STEPS}
 forceConstant 5.0
 outputCenters on
}}
harmonic {{
 name keep_tail_on_pore_axis
 colvars tail_radius
 centers 0.0
 forceConstant 0.5
}}
"""


def full_passage_from_csv(path):
    with Path(path).open() as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        if int(row["n_C1_below_trans"]) == 40:
            return dict(
                step=int(row["step"]),
                time_ns=float(row["time_ns"]),
                n_C1_below_trans=int(row["n_C1_below_trans"]),
                three_prime_C1_z_A=float(row["three_prime_C1_z_A"]),
                trans_reference_z_A=float(row["trans_reference_z_A"]),
            )
    return None


def require_threaded(result, label):
    endpoint = result["endpoint"]
    if endpoint["n_C1_below_constriction"] < 3:
        raise RuntimeError(
            f"{label}: fewer than three nucleotide C1 atoms are threaded"
        )
    if endpoint["minimum_DNA_protein_heavy_A"] <= 0.8:
        raise RuntimeError(f"{label}: DNA/protein overlap")


def run(sequence, config, models):
    job_id = os.environ["SLURM_JOB_ID"]
    work = ROOT / f"{sequence}-run-{job_id}"
    work.mkdir(exist_ok=False)
    snapshot = work / "code-snapshot"
    snapshot.mkdir()
    for name in [
        "README.md", "run.py", "stage.namd",
        f"run_{sequence}.sbatch",
    ]:
        shutil.copyfile(ROOT / name, snapshot / name)

    source = config["source"]
    auditor = base.Auditor(sequence, config["system"], source)
    source_metrics = auditor.restart_metrics(source)
    base.metrics_safe(source_metrics)
    plan = dict(
        sequence=config["label"], pore="7AHL", KCl_M=3,
        source=str(source.relative_to(WORKSPACE)),
        source_metrics=source_metrics,
        deep_start_definition="at least 3 nucleotide C1 atoms below constriction",
        high_field_z=FIELD,
        nominal_voltage_V=base.voltage(FIELD),
        field_only_trial_ns=TRIAL_STEPS * 1e-6,
        centered_tail_setup_ns=CENTER_STEPS * 1e-6,
        tail_SMD_ns=PULL_STEPS * 1e-6,
        tail_target_below_trans_A=TAIL_TARGET_BELOW_TRANS_A,
        interpretation=(
            "field-only trial followed, if needed, by nonequilibrium "
            "accelerated SMD; SMD is not voltage-only kinetics"
        ),
    )
    base.write_json(work / "run_plan.json", plan)
    base.event(
        work, "job_start", job_id=job_id,
        node=os.environ.get("SLURMD_NODENAME"), models=models, plan=plan,
    )

    current = source
    if config["guide_steps"]:
        geometry = group_geometry(auditor, current, dna_group(auditor, "lead"))
        guide = work / "position_lead.colvars"
        guide.write_text(lead_config(
            auditor, geometry["relative_z_A"],
            DEEP_LEAD_TARGET_A, config["guide_steps"],
        ))
        current, result = base.run_stage(
            work, sequence, config, "position_lead", current,
            config["guide_steps"], 0.0, guide,
        )
        require_threaded(result, "position_lead")

    auditor = base.Auditor(sequence, config["system"], current)
    geometry = group_geometry(auditor, current, dna_group(auditor, "lead"))
    hold = work / "hold_threaded.colvars"
    hold.write_text(lead_hold_config(auditor, geometry["relative_z_A"]))
    current, hold_result = base.run_stage(
        work, sequence, config, "hold_threaded", current,
        HOLD_STEPS, 0.0, hold,
    )
    require_threaded(hold_result, "hold_threaded")

    current, trial_result = base.run_stage(
        work, sequence, config, "field_trial", current,
        TRIAL_STEPS, FIELD, None,
    )
    trial_passage = full_passage_from_csv(
        work / "field_trial.metrics.csv"
    )
    if trial_passage:
        summary = dict(
            **plan, complete_passage_observed=True,
            passage_mode="field_only", first_passage=trial_passage,
            final_prefix=str(current.relative_to(WORKSPACE)),
        )
        base.write_json(work / "workflow_validation.json", summary)
        base.event(work, "workflow_complete", **summary)
        return

    if trial_result["endpoint"]["n_C1_below_constriction"] < 3:
        auditor = base.Auditor(sequence, config["system"], current)
        geometry = group_geometry(
            auditor, current, dna_group(auditor, "lead")
        )
        reinsert = work / "reinsert_lead.colvars"
        reinsert.write_text(lead_config(
            auditor, geometry["relative_z_A"],
            DEEP_LEAD_TARGET_A, REINSERT_STEPS,
        ))
        current, reinsert_result = base.run_stage(
            work, sequence, config, "reinsert_lead", current,
            REINSERT_STEPS, 0.0, reinsert,
        )
        require_threaded(reinsert_result, "reinsert_lead")

    auditor = base.Auditor(sequence, config["system"], current)
    lead = group_geometry(auditor, current, dna_group(auditor, "lead"))
    tail = group_geometry(auditor, current, dna_group(auditor, "tail"))
    center_config = work / "center_tail.colvars"
    center_config.write_text(tail_center_config(
        auditor, lead["relative_z_A"], tail
    ))
    current, center_result = base.run_stage(
        work, sequence, config, "center_tail", current,
        CENTER_STEPS, 0.0, center_config,
    )
    require_threaded(center_result, "center_tail")

    auditor = base.Auditor(sequence, config["system"], current)
    tail = group_geometry(auditor, current, dna_group(auditor, "tail"))
    target = (
        tail["trans_relative_z_A"] - TAIL_TARGET_BELOW_TRANS_A
    )
    pull_config = work / "pull_tail.colvars"
    pull_config.write_text(tail_pull_config(auditor, tail, target))
    base.event(
        work, "accelerated_SMD_plan",
        tail_start=tail, tail_target_relative_z_A=target,
        axial_shift_A=tail["relative_z_A"] - target,
        steps=PULL_STEPS, duration_ns=PULL_STEPS * 1e-6,
    )
    current, pull_result = base.run_stage(
        work, sequence, config, "pull_tail", current,
        PULL_STEPS, FIELD, pull_config,
    )
    first_passage = full_passage_from_csv(
        work / "pull_tail.metrics.csv"
    )
    summary = dict(
        **plan,
        complete_passage_observed=first_passage is not None,
        passage_mode="accelerated_tail_SMD",
        first_passage=first_passage,
        final_prefix=str(current.relative_to(WORKSPACE)),
        final_endpoint=pull_result["endpoint"],
        minimum_DNA_protein_heavy_A=(
            pull_result["minimum_DNA_protein_heavy_A"]
        ),
    )
    base.write_json(work / "workflow_validation.json", summary)
    base.event(work, "workflow_complete", **summary)


def check(sequence, config):
    auditor = base.Auditor(sequence, config["system"], config["source"])
    metrics = auditor.restart_metrics(config["source"])
    lead = group_geometry(
        auditor, config["source"], dna_group(auditor, "lead")
    )
    tail = group_geometry(
        auditor, config["source"], dna_group(auditor, "tail")
    )
    return dict(
        sequence=config["label"],
        source=str(config["source"].relative_to(WORKSPACE)),
        source_metrics=metrics, lead_geometry=lead, tail_geometry=tail,
        guide_steps=config["guide_steps"],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence", choices=sorted(RUNS), required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    config = RUNS[args.sequence]
    for suffix in ["coor", "vel", "xsc"]:
        path = Path(str(config["source"]) + "." + suffix)
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    if args.check_only:
        print(json.dumps(check(args.sequence, config), indent=2))
        return
    models = base.validate_allocation()
    run(args.sequence, config, models)


if __name__ == "__main__":
    main()
