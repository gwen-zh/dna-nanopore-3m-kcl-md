# Read one original frame, not an atom-wise time-unwrapped analysis export.
if {[catch {
    package require pbctools
    lassign $argv psf dcd frame output
    mol new $psf type psf waitfor all
    mol addfile $dcd type dcd first $frame last $frame waitfor all
    pbc join fragment -bondlist -sel nucleic -now
    set dna [atomselect top nucleic]
    set xyz [$dna get {x y z}]
    set indices [$dna get index]
    set coords [dict create]
    foreach i $indices r $xyz { dict set coords $i $r }
    set maximum 0.0
    set count 0
    foreach i $indices neighbors [$dna getbonds] {
        foreach j $neighbors {
            if {$i >= $j || ![dict exists $coords $j]} { continue }
            set d [veclength [vecsub [dict get $coords $i] [dict get $coords $j]]]
            if {$d < 0.5 || $d > 2.2} { error "Invalid donor bond $i $j length $d" }
            set maximum [expr {max($maximum,$d)}]
            incr count
        }
    }
    if {[$dna num] != 1279 || $count < 1278} { error "Incomplete DNA donor" }
    $dna writepdb $output
    puts "DONOR_VALIDATED raw_dcd_frame=$frame atoms=[$dna num] bonds=$count max_bond_A=$maximum output=$output"
} message]} {
    puts stderr "FATAL donor extraction: $message"
    exit 1
}
quit
