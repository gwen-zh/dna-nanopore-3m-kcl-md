if {[catch {
    source ../../common/fast_autoionize.tcl
    autoionize -psf solvated_noions.psf -pdb solvated_noions.pdb \
        -sc 3.0 -cation POT -anion CLA -seg ION -from 5.0 -between 3.5 -o system
    set all [atomselect top all]
    set dna [atomselect top "segname AN1"]
    set k [atomselect top "resname POT"]
    set cl [atomselect top "resname CLA"]
    puts "SYSTEM atoms=[$all num] dna_atoms=[$dna num] potassium=[$k num] chloride=[$cl num]"
} message]} {
    puts stderr "FATAL ionization: $message"
    exit 1
}
quit
