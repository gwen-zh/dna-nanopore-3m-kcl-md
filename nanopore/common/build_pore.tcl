package require psfgen
package require solvate
package require autoionize

if {[llength $argv] != 13} {
    error "usage: pore_psf pore_pdb donor_pdb base remove_old_dna target_z minx miny minz maxx maxy maxz label"
}

lassign $argv pore_psf pore_pdb donor_pdb base remove_old_dna target_z \
    minx miny minz maxx maxy maxz label

if {$base ne "ADE" && $base ne "THY"} {
    error "base must be ADE or THY"
}

set toppar /opt/charmm/toppar
topology ${toppar}/top_all36_na.rtf

# Make a clean pore template. The legacy 7AHL source contains an old ssDNA
# segment and ions; the 3B07 template contains neither, but the catch calls are
# harmless and keep both build routes identical.
readpsf $pore_psf
coordpdb $pore_pdb
catch {delatom ION}
catch {delatom ION1}
if {$remove_old_dna} {
    delatom AN1
}
writepsf pore_clean.psf
writepdb pore_clean.pdb

# Orient the high-salt DNA with its 5' end toward the pore and its 1->40
# end-to-end vector along +Z. The external helper performs a deterministic
# azimuth/height search, protein-clash test, and periodic-box clearance test.
set placement [exec python3 ../../common/place_dna.py \
    pore_clean.pdb $donor_pdb placed_dna.pdb $target_z \
    $minx $miny $minz $maxx $maxy $maxz $label 10.0]
puts $placement

# Rebuild the selected 40-mer with the CHARMM36 deoxyribose patches, then add
# it to the clean pore/membrane/solvent template.
resetpsf
topology ${toppar}/top_all36_na.rtf
readpsf pore_clean.psf
coordpdb pore_clean.pdb
segment AN1 {
    first 5TER
    last 3TER
    for {set resid 1} {$resid <= 40} {incr resid} {
        residue $resid $base
    }
}
coordpdb placed_dna.pdb AN1
patch DEO5 AN1:1
for {set resid 2} {$resid <= 40} {incr resid} {
    patch DEOX AN1:$resid
}
regenerate angles dihedrals
guesscoord
writepsf combined_raw.psf
writepdb combined_raw.pdb

# Reject protein/lipid overlaps and remove complete pre-existing waters that
# are displaced by the newly inserted, folded DNA.
mol new combined_raw.psf type psf waitfor all
mol addfile combined_raw.pdb type pdb waitfor all
set dna [atomselect top "segname AN1"]
set protein [atomselect top "protein"]
set lipid [atomselect top "lipid"]
set protein_clashes [llength [lindex [measure contacts 1.2 $dna $protein] 0]]
set lipid_clashes [llength [lindex [measure contacts 1.2 $dna $lipid] 0]]
if {$protein_clashes > 0 || $lipid_clashes > 0} {
    error "invalid placement: protein_clashes=$protein_clashes lipid_clashes=$lipid_clashes"
}
set badwater [atomselect top "water and same residue as within 2.4 of segname AN1"]
set bad_pairs [lsort -unique [$badwater get {segname resid}]]
puts "REMOVING [llength $bad_pairs] displaced waters"
$badwater delete
$dna delete
$protein delete
$lipid delete
mol delete top

resetpsf
topology ${toppar}/top_all36_na.rtf
readpsf combined_raw.psf
coordpdb combined_raw.pdb
foreach pair $bad_pairs {
    delatom [lindex $pair 0] [lindex $pair 1]
}
writepsf pre_solv.psf
writepdb pre_solv.pdb

set minmax [list [list $minx $miny $minz] [list $maxx $maxy $maxz]]
solvate pre_solv.psf pre_solv.pdb -minmax $minmax -s NW -o solvated_noions
if {![file exists solvated_noions.psf] || ![file exists solvated_noions.pdb]} {
    error "solvate did not produce solvated_noions.psf/pdb"
}
source ../../common/fast_autoionize.tcl
autoionize -psf solvated_noions.psf -pdb solvated_noions.pdb \
    -sc 3.0 -cation POT -anion CLA -seg ION -from 5.0 -between 3.5 \
    -o system
if {![file exists system.psf] || ![file exists system.pdb]} {
    error "autoionize did not produce system.psf/pdb"
}

mol new system.psf type psf waitfor all
mol addfile system.pdb type pdb waitfor all
set nall [[atomselect top all] num]
set ndna [[atomselect top "segname AN1"] num]
set nk [[atomselect top "resname POT"] num]
set ncl [[atomselect top "resname CLA"] num]
set nwater [[atomselect top water] num]
if {$nk == 0 || $ncl == 0} {
    error "ionization validation failed: potassium=$nk chloride=$ncl"
}
puts "SYSTEM label=$label atoms=$nall dna_atoms=$ndna water_atoms=$nwater potassium=$nk chloride=$ncl box=$minmax"
quit
