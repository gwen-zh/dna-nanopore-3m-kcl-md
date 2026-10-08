#!/usr/bin/env python3
"""Extend the validated 5 ns dT40 3 M KCl production to 50 ns."""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
WORKSPACE = ROOT.parent
DT_ROOT = WORKSPACE / "dt-image-repair-20260915"
BOX = DT_ROOT / "box180"
SOURCE = DT_ROOT / "run-57330935/prod5"
NAMD = "/opt/NAMD/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2"
DCD_STRIDE = 5000

sys.path.insert(0, str(DT_ROOT))
from checks import Checker, live_dna_frames  # noqa: E402


CHUNKS = [
    dict(label="prod20", start_step=2_500_000, run_steps=7_500_000,
         final_step=10_000_000, final_time_ns=20.0, seed="202609273"),
    dict(label="prod35", start_step=10_000_000, run_steps=7_500_000,
         final_step=17_500_000, final_time_ns=35.0, seed="202609274"),
    dict(label="prod50", start_step=17_500_000, run_steps=7_500_000,
         final_step=25_000_000, final_time_ns=50.0, seed="202609275"),
]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1_048_576), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def record(work, label, **fields):
    value = dict(utc=datetime.now(timezone.utc).isoformat(), event=label, **fields)
    with (work / "events.jsonl").open("a") as stream:
        stream.write(json.dumps(value) + "\n")
    print(json.dumps(value), flush=True)


def terminate_own_process(process):
    process.terminate()
    try:
        process.wait(timeout=30)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()


def validate_allocation():
    if os.environ.get("SLURM_JOB_PARTITION") != "dept_gpu":
        raise RuntimeError("Only dept_gpu is authorized")
    if int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) != 8:
        raise RuntimeError("Exactly 8 CPU workers are required")
    models = subprocess.check_output(
        ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"], text=True
    ).splitlines()
    if not models or any(re.search(r"\bL40", model, re.I) for model in models):
        raise RuntimeError("A visible L40 GPU is not authorized")
    return models


def run_chunk(work, checker, source, chunk):
    label = chunk["label"]
    source_validation = checker.restart(source)
    if source_validation["stage_step"] != chunk["start_step"]:
        raise RuntimeError(
            f"{label}: source step {source_validation['stage_step']} does not match "
            f"expected {chunk['start_step']}"
        )
    output = work / label
    env = dict(
        os.environ,
        DT_BOX=str(BOX),
        CODE_DIR=str(ROOT),
        INPUT_PREFIX=str(source),
        OUTPUT_PREFIX=str(output),
        START_STEP=str(chunk["start_step"]),
        RUN_STEPS=str(chunk["run_steps"]),
        MD_SEED=chunk["seed"],
        OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
    )
    record(
        work,
        "stage_start",
        stage=label,
        source=str(source),
        source_validation=source_validation,
        start_step=chunk["start_step"],
        run_steps=chunk["run_steps"],
        final_step=chunk["final_step"],
        final_time_ns=chunk["final_time_ns"],
        timestep_fs=2,
        temperature_K=293,
        pressure_bar=1.01325,
        field=False,
        restraints=False,
    )

    seen = 0
    minimum_image_gap = float("inf")
    first = None
    last = None
    with (work / f"{label}.metrics.csv").open("x", newline="") as metric_stream, \
            (work / f"{label}.log").open("x") as log:
        writer = None

        def audit(strict=False):
            nonlocal seen, minimum_image_gap, first, last, writer
            try:
                for index, step, xyz, box in live_dna_frames(
                    str(output) + ".dcd", checker.n, 1279, seen, strict
                ):
                    expected_step = chunk["start_step"] + (index + 1) * DCD_STRIDE
                    if step != expected_step:
                        raise ValueError(
                            f"{label}: DCD step {step} does not match {expected_step}"
                        )
                    row = dict(
                        production_step=step,
                        production_time_ns=step * 2e-6,
                        **checker.dna_metrics(xyz, box),
                    )
                    if writer is None:
                        writer = csv.DictWriter(metric_stream, fieldnames=list(row))
                        writer.writeheader()
                    writer.writerow(row)
                    metric_stream.flush()
                    seen = index + 1
                    first = first or row
                    last = row
                    minimum_image_gap = min(
                        minimum_image_gap, row["nearest_periodic_image_A"]
                    )
                    if seen % 50 == 0:
                        record(
                            work,
                            "periodic_guard_pass",
                            stage=label,
                            frame_count=seen,
                            metrics=row,
                        )
            except EOFError:
                if strict:
                    raise

        process = subprocess.Popen(
            [NAMD, "+p8", "+devices", "0", "+idlepoll", str(ROOT / "stage.namd")],
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        try:
            while process.poll() is None:
                audit()
                time.sleep(5)
            audit(strict=True)
        except BaseException as exc:
            if process.poll() is None:
                terminate_own_process(process)
            record(work, "stage_stopped_by_guard", stage=label, reason=str(exc))
            raise

    text = (work / f"{label}.log").read_text(errors="replace")
    if (process.returncode != 0 or "End of program" not in text or
            re.search(r"FATAL ERROR|\bnan\b", text, re.I)):
        raise RuntimeError(f"NAMD failed during {label}")
    if not re.search(rf"^ENERGY:\s+{chunk['final_step']}\s", text, re.M):
        raise RuntimeError(f"{label}: final ENERGY record is missing")
    expected_frames = chunk["run_steps"] // DCD_STRIDE
    if seen != expected_frames or first is None or last is None:
        raise RuntimeError(
            f"{label}: expected {expected_frames} DCD frames, found {seen}"
        )
    restart_validation = checker.restart(output)
    if restart_validation["stage_step"] != chunk["final_step"]:
        raise RuntimeError(f"{label}: final restart step mismatch")
    result = dict(
        stage=label,
        source=str(source),
        output_prefix=str(output),
        frames=seen,
        first=first,
        endpoint=last,
        minimum_periodic_image_A=minimum_image_gap,
        restart_validation=restart_validation,
        final_production_time_ns=chunk["final_time_ns"],
    )
    write_json(work / f"{label}.validation.json", result)
    record(work, "stage_complete", **result)
    return output, result


def preflight(checker):
    source_validation = checker.restart(SOURCE)
    if source_validation["stage_step"] != 2_500_000:
        raise RuntimeError("The retained dT checkpoint is not the 5 ns endpoint")
    if source_validation["nearest_periodic_image_A"] <= 24:
        raise RuntimeError("The retained dT checkpoint fails the periodic-image gate")
    return dict(
        system="poly(dT)40 in 3 M KCl, explicit 180 A periodic box",
        source=str(SOURCE.relative_to(WORKSPACE)),
        source_validation=source_validation,
        existing_production_ns=5,
        additional_steps=22_500_000,
        additional_ns=45,
        final_production_ns=50,
        timestep_fs=2,
        chunks=CHUNKS,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    checker = Checker(BOX)
    plan = preflight(checker)
    if args.check_only:
        print(json.dumps(plan, indent=2))
        return

    models = validate_allocation()
    work = ROOT / ("run-" + os.environ["SLURM_JOB_ID"])
    work.mkdir(exist_ok=False)
    snapshot = work / "code-snapshot"
    snapshot.mkdir()
    for name in ["README.md", "common.namd", "stage.namd", "run.py", "run.sbatch"]:
        shutil.copyfile(ROOT / name, snapshot / name)

    manifest_files = [
        ROOT / "common.namd",
        ROOT / "stage.namd",
        ROOT / "run.py",
        ROOT / "run.sbatch",
        DT_ROOT / "checks.py",
        BOX / "system.psf",
        BOX / "system.pdb",
        BOX / "dna.psf",
        BOX / "restraints.pdb",
        BOX / "water_ions_namd.prm",
        Path(str(SOURCE) + ".coor"),
        Path(str(SOURCE) + ".vel"),
        Path(str(SOURCE) + ".xsc"),
    ]
    write_json(
        work / "source_manifest.json",
        [
            dict(
                path=str(path.relative_to(WORKSPACE)),
                bytes=path.stat().st_size,
                sha256=sha256(path),
            )
            for path in manifest_files
        ],
    )
    write_json(work / "run_plan.json", plan)
    record(
        work,
        "job_start",
        job_id=os.environ["SLURM_JOB_ID"],
        node=os.environ.get("SLURMD_NODENAME"),
        partition=os.environ.get("SLURM_JOB_PARTITION"),
        models=models,
        threads=8,
        plan=plan,
    )

    source = SOURCE
    results = []
    try:
        for chunk in CHUNKS:
            source, result = run_chunk(work, checker, source, chunk)
            results.append(result)
        summary = dict(
            production_prefix=str(source),
            total_production_ns=50,
            added_production_ns=45,
            completed_chunks=[item["stage"] for item in results],
            result=results[-1],
        )
        write_json(work / "workflow_validation.json", summary)
        record(work, "workflow_complete", **summary)
    except BaseException as exc:
        record(
            work,
            "workflow_stopped",
            reason=str(exc),
            completed_chunks=[item["stage"] for item in results],
            automatic_retry=False,
        )
        raise


if __name__ == "__main__":
    main()
