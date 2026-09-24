#!/usr/bin/env python3
"""One GPU, preserved DNA topology, original preparation then a gated 5 ns block."""
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

from checks import Checker, live_dna_frames

ROOT = Path(__file__).resolve().parent
BOX = ROOT / 'box180'
NAMD = '/opt/NAMD/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2'
assert os.environ['SLURM_JOB_PARTITION'] == 'dept_gpu'
assert int(os.environ['SLURM_CPUS_PER_TASK']) == 8
models = subprocess.check_output(['nvidia-smi', '--query-gpu=name', '--format=csv,noheader'], text=True)
assert not re.search(r'\bL40', models, re.I)
WORK = ROOT / ('run-' + os.environ['SLURM_JOB_ID'])
WORK.mkdir(exist_ok=False)
SNAP = WORK / 'code-snapshot'
SNAP.mkdir()
hashes = {}
for name in ['run.py', 'checks.py', 'common.namd', 'stage.namd', 'build.tcl', 'run.sbatch']:
    shutil.copyfile(ROOT / name, SNAP / name)
    hashes[name] = hashlib.sha256((SNAP / name).read_bytes()).hexdigest()
(WORK / 'code_sha256.json').write_text(json.dumps(hashes, indent=2) + '\n')

def record(event, **fields):
    data = {'utc': datetime.now(timezone.utc).isoformat(), 'event': event, **fields}
    print(json.dumps(data), flush=True)
    with (WORK / 'events.jsonl').open('a') as stream:
        stream.write(json.dumps(data) + '\n')

base_env = dict(os.environ, DT_BUILD_DIR=str(BOX), DT_SCRIPT_DIR=str(SNAP))
record('job_start', job_id=os.environ['SLURM_JOB_ID'], node=os.environ.get('SLURMD_NODENAME'),
       model_names=models.splitlines(), cpus=8, gpus=1,
       purpose='dT DNA-only production in an enlarged periodic box')
with (WORK / 'build.log').open('x') as stream:
    rc = subprocess.run(['/opt/bin/vmd', '-dispdev', 'text', '-e', str(SNAP / 'build.tcl')],
                        env=base_env, stdout=stream, stderr=subprocess.STDOUT).returncode
build_log = (WORK / 'build.log').read_text(errors='replace')
if rc != 0 or 'BUILD_COMPLETE' not in build_log or 'BUILD_FAILED' in build_log:
    raise RuntimeError('Solvent build failed; outputs retained, MD not started')
checker = Checker(BOX)
initial = checker.initial()
(WORK / 'build_validation.json').write_text(json.dumps(initial, indent=2) + '\n')
record('build_validated', metrics=initial)

def terminate_own_process(proc):
    proc.terminate()
    try:
        proc.wait(timeout=30)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()

def stage(label, source, mode, steps, scale):
    prefix = WORK / label
    env = dict(base_env, OUTPUT_PREFIX=str(prefix), INPUT_PREFIX=str(source or ''),
               STAGE_MODE=mode, STAGE_STEPS=str(steps), RESTRAINT_SCALE=str(scale))
    record('stage_start', stage=label, mode=mode, steps=steps, timestep_fs=2,
           restraint_scale=scale, input=str(source), output=str(prefix), field=False, threads=8)
    seen = 0
    with (WORK / f'{label}.metrics.csv').open('x', newline='') as metric_stream, \
         (WORK / f'{label}.log').open('x') as log:
        writer = None

        def audit(strict=False):
            nonlocal seen, writer
            try:
                for index, step, xyz, box in live_dna_frames(str(prefix) + '.dcd', checker.n, 1279, seen, strict):
                    row = {'stage_step': step, 'stage_time_ns': step * 2e-6,
                           **checker.dna_metrics(xyz, box)}
                    if writer is None:
                        writer = csv.DictWriter(metric_stream, fieldnames=list(row))
                        writer.writeheader()
                    writer.writerow(row)
                    metric_stream.flush()
                    seen = index + 1
                    if row['nearest_periodic_image_A'] < 36 or seen % 50 == 0:
                        record('periodic_guard_pass', stage=label, frame_count=seen, metrics=row)
            except EOFError:
                if strict:
                    raise

        proc = subprocess.Popen([NAMD, '+p8', '+devices', '0', '+idlepoll', str(SNAP / 'stage.namd')],
                                env=env, stdout=log, stderr=subprocess.STDOUT)
        try:
            while proc.poll() is None:
                if mode != 'min':
                    audit()
                time.sleep(5)
            if mode != 'min':
                audit(strict=True)
                if seen != steps // 5000:
                    raise ValueError(f'Missing DCD frames: {seen}, expected {steps // 5000}')
        except BaseException as exc:
            if proc.poll() is None:
                terminate_own_process(proc)
            record('stage_stopped_by_guard', stage=label, reason=str(exc))
            raise
    text = (WORK / f'{label}.log').read_text(errors='replace')
    if proc.returncode != 0 or 'End of program' not in text or re.search(r'FATAL ERROR|\bnan\b', text, re.I):
        record('stage_failed', stage=label, returncode=proc.returncode)
        raise RuntimeError('NAMD failed; no single-thread fallback or automatic unsafe continuation')
    if not re.search(rf'^ENERGY:\s+{steps}\s', text, re.M):
        raise RuntimeError('Expected final ENERGY step missing')
    metrics = checker.restart(prefix)
    if metrics['stage_step'] != steps:
        raise RuntimeError('Final restart step mismatch')
    (WORK / f'{label}.validation.json').write_text(json.dumps(metrics, indent=2) + '\n')
    record('stage_complete', stage=label, audited_frames=seen, metrics=metrics)
    return prefix

# Match the existing preparation: 50k min, 0.5 ns heat, then 2+2.5+5 ns equil.
source = None
for label, mode, steps, scale in [
    ('min', 'min', 50000, 1.0), ('heat', 'heat', 250000, 1.0),
    ('eq1', 'md', 1000000, .5), ('eq2', 'md', 1250000, .1),
    ('eq3', 'md', 2500000, 0), ('prod5', 'md', 2500000, 0),
]:
    source = stage(label, source, mode, steps, scale)
record('first_sampling_block_complete', production_ns=5, production_prefix=str(source),
       donor_ready=False, next_action='Review image clearance, salt density and conformational drift before extending sampling or selecting new pore donors')
