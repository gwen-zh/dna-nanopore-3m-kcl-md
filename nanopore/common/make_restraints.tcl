mol new system.psf type psf waitfor all
mol addfile system.pdb type pdb waitfor all
set all [atomselect top all]
$all set beta 0.0
set ca [atomselect top "protein and name CA"]
set dna [atomselect top "segname AN1 and backbone"]
set lipidp [atomselect top "lipid and name P"]
$ca set beta 1.0
$dna set beta 1.0
$lipidp set beta 0.5
$all writepdb restraints_strong.pdb
$all set beta 0.0
$ca set beta 0.5
$lipidp set beta 0.1
$all writepdb restraints_prod.pdb
puts "RESTRAINTS protein_CA=[$ca num] DNA_backbone=[$dna num] lipid_P=[$lipidp num]"
quit
