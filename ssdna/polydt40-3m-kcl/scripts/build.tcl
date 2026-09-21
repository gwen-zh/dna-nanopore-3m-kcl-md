# Preserve the extracted DNA PSF; only solvent and ions are rebuilt.
proc build_box {} {
    global env
    cd $env(DT_REPAIR_BUILD)
    if {[file exists system.psf] || [file exists solvated.psf]} {
        error "Refusing to overwrite an existing solvent build"
    }
    package require psfgen
    package require solvate
    package require autoionize
    expr {srand(9152026)}
    solvate dna.psf dna.pdb -minmax {{-90 -90 -90} {90 90 90}} -o solvated
    # Retain the previous -sc 3.0 preparation convention. Actual NPT box
    # concentrations must be measured, not claimed to be exactly 3.000 M.
    autoionize -psf solvated.psf -pdb solvated.pdb \
        -sc 3.0 -cation POT -anion CLA -seg ION -from 5.0 -between 3.5 -o system
    mol delete all
    mol new system.psf type psf waitfor all
    mol addfile system.pdb type pdb waitfor all
    set all [atomselect top all]
    set dna [atomselect top "segname DNA"]
    if {[$dna num] != 1279} { error "DNA atom count changed" }
    set charge [vecsum [$all get charge]]
    if {abs($charge) > 0.001} { error "Non-neutral box: $charge" }
    $all set beta 0
    set bb [atomselect top "segname DNA and backbone"]
    $bb set beta 1.0
    $all writepdb restraints.pdb
    puts "BUILD_COMPLETE atoms=[$all num] DNA=[$dna num] restrained=[$bb num] charge=$charge"
}
if {[catch {build_box} message]} {
    puts stderr "BUILD_FAILED: $message"
    puts stderr $::errorInfo
    exit 1
}
quit
