# PyMOL display of real sampled representative frames from the matched 45.01–50.00 ns
# dA40/dT40 analysis.  PCA is used only for a rigid-body display alignment;
# internal coordinates and conformations are unchanged.
reinitialize
load ../dA40_representative.pdb, dA40
load ../dT40_representative.pdb, dT40

python
import numpy as np
from pymol import cmd

def principal_display_alignment(name, x_shift):
    coords = cmd.get_coords(name, 1)
    center = coords.mean(axis=0)
    centered = coords - center
    covariance = np.cov(centered.T)
    _, vectors = np.linalg.eigh(covariance)
    vectors = vectors[:, ::-1]
    if np.linalg.det(vectors) < 0:
        vectors[:, -1] *= -1
    # PC1 is displayed vertically, PC2 horizontally and PC3 in depth.
    transformed = centered @ vectors[:, [1, 0, 2]]
    transformed[:, 0] += x_shift
    cmd.load_coords(transformed, name, 1)

principal_display_alignment("dA40", -38.0)
principal_display_alignment("dT40", 38.0)
python end

bg_color white
set ray_opaque_background, on
set orthoscopic, on
set depth_cue, 0
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 0
set specular, 0.18
set ambient, 0.48
set direct, 0.55
set cartoon_smooth_loops, on
set cartoon_ladder_mode, 1
set cartoon_ring_mode, 3
set cartoon_nucleic_acid_mode, 4
set cartoon_tube_radius, 0.65
set stick_radius, 0.19
set sphere_scale, 0.65
set label_font_id, 7
set label_size, 24
set label_color, black
set label_outline_color, white

hide everything, all
show cartoon, dA40 or dT40
select dA_bases, dA40 and not name P+OP1+OP2+O5'+C5'+H5'+H5''+C4'+H4'+O4'+C3'+H3'+O3'+C2'+H2'+H2''+C1'+H1'+H5T+H3T
select dT_bases, dT40 and not name P+OP1+OP2+O5'+C5'+H5'+H5''+C4'+H4'+O4'+C3'+H3'+O3'+C2'+H2'+H2''+C1'+H1'+H5T+H3T
show sticks, dA_bases or dT_bases

set_color da_backbone, [0.10, 0.55, 0.78]
set_color da_base, [0.35, 0.78, 0.95]
set_color dt_backbone, [0.92, 0.40, 0.10]
set_color dt_base, [1.00, 0.68, 0.28]
color da_backbone, dA40
color da_base, dA_bases
color dt_backbone, dT40
color dt_base, dT_bases

select termini, (dA40 or dT40) and name C1' and (resi 1 or resi 40)
show spheres, termini
color forest, (dA40 or dT40) and resi 1 and name C1'
color red, (dA40 or dT40) and resi 40 and name C1'
label dA40 and resi 1 and name C1', "5'"
label dA40 and resi 40 and name C1', "3'"
label dT40 and resi 1 and name C1', "5'"
label dT40 and resi 40 and name C1', "3'"
pseudoatom dA_title, pos=[-38, 60, 0], label="poly(dA)40"
pseudoatom dT_title, pos=[38, 60, 0], label="poly(dT)40"
hide spheres, dA_title or dT_title
set label_size, 34, dA_title or dT_title
set label_position, [0, 0, 0], dA_title or dT_title

orient dA40 or dT40
zoom dA40 or dT40, 7
viewport 2400, 1400
png dA40_dT40_representatives_pymol.png, 2400, 1400, 300, 1
save dA40_dT40_representatives.pse

disable dT40
disable dA_title
disable dT_title
hide labels, dT40
orient dA40
zoom dA40, 7
viewport 1600, 1400
png dA40_representative_pymol.png, 1600, 1400, 300, 1

enable dT40
show labels, dT40 and name C1' and (resi 1 or resi 40)
disable dA40
hide labels, dA40
orient dT40
zoom dT40, 7
png dT40_representative_pymol.png, 1600, 1400, 300, 1

enable dA40
show labels, dA40 and name C1' and (resi 1 or resi 40)
enable dA_title
enable dT_title
orient dA40 or dT40
zoom dA40 or dT40, 7
save dA40_dT40_representatives.pse
quit
