#!/usr/bin/env python3
"""Render every saved translocation frame and encode a validated H.264 movie."""

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

import imageio_ffmpeg
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
BOLD = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"


def compose(picture, row, summary):
    width, height = 1280, 960
    canvas = Image.new("RGB", (width, height), "white")
    picture = picture.convert("RGB")
    picture.thumbnail((950, 830), Image.Resampling.LANCZOS)
    canvas.paste(picture, ((950 - picture.width) // 2, 90 + (830 - picture.height) // 2))
    draw = ImageDraw.Draw(canvas)

    def label(x, y, text, size=22, color="#24384a", bold=False):
        draw.text((x, y), text, font=ImageFont.truetype(BOLD if bold else FONT, size), fill=color)

    draw.rectangle((0, 0, width, 80), fill="#153348")
    label(28, 14, f"{summary['sequence']} through {summary['pore']} | nominal 3 M KCl", 27, "white", True)
    label(28, 49, "High-field ordered-SMD pathway | complete 40/40 passage", 17, "#cae2ee")
    draw.line((959, 110, 959, 875), fill="#dfe6ea", width=2)
    label(985, 116, "PATH TIME", 16, "#647a89", True)
    label(985, 146, f"{row['path_time_ns']:.3f} ns", 31, bold=True)
    label(985, 205, "Passage progress", 19, bold=True)
    label(985, 240, f"cis: {int(row['n_C1_below_cis'])} / 40", 19)
    label(985, 271, f"constriction: {int(row['n_C1_below_constriction'])} / 40", 18)
    label(985, 302, f"trans: {int(row['n_C1_below_trans'])} / 40", 19, bold=True)
    label(985, 359, "DNA conformation", 19, bold=True)
    label(985, 394, f"Rg: {row['Rg_A']:.1f} Å", 19)
    label(985, 425, f"End-to-end: {row['end_to_end_A']:.1f} Å", 18)
    label(985, 482, "Interactions", 19, bold=True)
    label(985, 517, f"Atom contacts: {int(row['DNA_protein_atom_contact_pairs_3p5A'])}", 18)
    label(985, 548, f"Candidate H bonds: {int(row['candidate_hydrogen_bonds'])}", 17)
    label(985, 579, f"K+ bridges: {int(row['potassium_mediated_DNA_acidic_residue_bridges'])}", 17)
    for y, color, text in [
        (650, "#f28f1f", "DNA"),
        (687, "#518cb5", "Pore protein"),
        (724, "#d62e29", "5′ nucleotide"),
        (761, "#30a663", "3′ nucleotide"),
        (798, "#a1a6ad", "Membrane P atoms"),
    ]:
        draw.ellipse((987, y + 3, 1000, y + 16), fill=color)
        label(1011, y, text, 17)
    label(28, 903, "Actual saved frames | no interpolation | water and ions hidden | pore-aligned view", 16, "#566b7b")
    label(28, 930, "Nonequilibrium pathway: interaction counts are geometric, not binding energies or experimental kinetics.", 15, "#566b7b")
    end = float(summary["cumulative_simulated_path_ns"])
    draw.rectangle((0, 82, int(width * row["path_time_ns"] / end), 87), fill="#ee962c")
    return canvas


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--test-frame", type=int)
    parser.add_argument("--encode-only", action="store_true")
    args = parser.parse_args()
    folder = args.input.resolve()
    config = json.loads((folder / "render.json").read_text())
    summary = json.loads((folder / "summary.json").read_text())
    table = pd.read_csv(folder / "timeseries.csv")
    frames = folder / "raytrace"
    frames.mkdir(exist_ok=True)
    first = 0 if args.test_frame is None else args.test_frame
    last = len(table) - 1 if args.test_frame is None else args.test_frame

    if not args.encode_only:
        env = dict(
            os.environ,
            VMDFORCECPUCOUNT="2",
            VMDCUDA="0",
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
        )
        span = max(config["z_bounds"][1] - config["z_bounds"][0], 150.0)
        command = [
            os.environ.get("VMD_BIN", "/opt/bin/vmd"), "-dispdev", "text",
            "-e", str(ROOT / "render.tcl"), "-args", str(folder), str(frames),
            *map(str, config["center"]), str(span), str(first), str(last),
        ]
        log = folder / ("render_test.log" if args.test_frame is not None else "render.log")
        with log.open("w") as stream:
            process = subprocess.Popen(
                command, env=env, stdin=subprocess.PIPE, stdout=stream, stderr=subprocess.STDOUT
            )
            return_code = process.wait()
            process.stdin.close()
        text = log.read_text(errors="replace")
        if return_code != 0 or "RENDER_COMPLETE" not in text or "RENDER_FAILED" in text:
            raise RuntimeError(f"VMD rendering failed; inspect {log}")
        checks = pd.read_csv(folder / f"render_checks_{first}_{last}.csv")
        for _, check in checks.iterrows():
            source = table.iloc[int(check["frame"])]
            if not np.isclose(source["Rg_A"], check["rg_A"], atol=1e-3, rtol=1e-5):
                raise ValueError("Independent VMD Rg measurement disagrees")
            if not np.isclose(source["end_to_end_A"], check["end_to_end_A"], atol=1e-3, rtol=1e-5):
                raise ValueError("Independent VMD end-to-end measurement disagrees")

    if args.test_frame is not None:
        compose(Image.open(frames / f"frame_{first:04d}.tga"), table.iloc[first], summary).save(folder / "render_test.png")
        print(folder / "render_test.png")
        return

    movie = folder / "trajectory.mp4"
    writer = imageio_ffmpeg.write_frames(
        str(movie), (1280, 960), fps=10, codec="libx264", pix_fmt_out="yuv420p",
        quality=8, macro_block_size=16, output_params=["-movflags", "+faststart", "-threads", "2"],
    )
    writer.send(None)
    distinct = set()
    try:
        for index, row in table.iterrows():
            picture = Image.open(frames / f"frame_{index:04d}.tga").convert("RGB")
            pixels = np.asarray(picture)
            if np.mean(pixels.min(axis=2) < 220) < 0.005:
                raise ValueError(f"Blank molecular rendering at frame {index}")
            distinct.add(hashlib.sha256(pixels.tobytes()).hexdigest())
            canvas = compose(picture, row, summary)
            if index == len(table) - 1:
                canvas.save(folder / "poster.png")
            writer.send(np.asarray(canvas))
    finally:
        writer.close()
    if len(distinct) < len(table) * 0.90:
        raise ValueError("Frozen or excessively repeated molecular animation frames")

    reader = imageio_ffmpeg.read_frames(str(movie), pix_fmt="rgb24")
    next(reader)
    decoded = sum(1 for _ in reader)
    if decoded != len(table):
        raise ValueError(f"Encoded {decoded} frames, expected {len(table)}")
    validation = {
        "codec": "H.264",
        "container": "MP4",
        "dimensions": [1280, 960],
        "fps": 10,
        "frames": decoded,
        "duration_seconds": decoded / 10,
        "cumulative_simulated_path_ns": summary["cumulative_simulated_path_ns"],
        "interpolation": False,
        "decoded_all_frames": True,
        "distinct_molecular_frames": len(distinct),
        "water_and_ions_hidden": True,
        "protein_front_chains_transparent": True,
        "ffmpeg_version": imageio_ffmpeg.get_ffmpeg_version(),
    }
    (folder / "movie_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps(validation, indent=2), flush=True)


if __name__ == "__main__":
    main()
