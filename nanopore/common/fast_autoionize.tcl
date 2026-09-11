# Keep the installed autoionize count, topology and water-replacement code.
# Replace only its quadratic site-search loop, with guarded source markers.
package require autoionize
set ::fast_ion_helper [file join [file dirname [info script]] fast_ion_positions.py]
proc ::fast_ion_sites {count from between} {
    set water [atomselect top "water and noh"]
    set fh [open ion_candidates.tsv w]
    foreach i [$water get index] r [$water get {x y z}] {
        puts $fh "$i [join $r { }]"
    }
    close $fh
    $water delete
    set solute [atomselect top "not water"]
    set fh [open ion_solute.tsv w]
    foreach r [$solute get {x y z}] { puts $fh [join $r { }] }
    close $fh
    $solute delete
    set sites [exec python3 $::fast_ion_helper $count $from $between $::env(SYSTEM_LABEL)]
    if {[llength $sites] != $count} { error "Wrong fast ion site count" }
    return $sites
}
set original_body [info body ::autoionize::autoionize_core]
set start [string first {    set nTries 0} $original_body]
set end [string first {    puts "Autoionize) Obtained positions for $nIons ions."} $original_body]
if {$start < 0 || $end <= $start} { error "Unsupported autoionize source version" }
set replacement {    set ionList [::fast_ion_sites $nIons $from $between]
}
proc ::autoionize::autoionize_core {args} [string replace $original_body $start [expr {$end-1}] $replacement]
puts "FAST_AUTOIONIZE enabled: original ion counts/topology, periodic distance checks"
