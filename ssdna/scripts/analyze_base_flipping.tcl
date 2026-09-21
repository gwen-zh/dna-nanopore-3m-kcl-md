set psf [lindex $argv 0]
set pdb [lindex $argv 1]
set dcd [lindex $argv 2]
set baseType [lindex $argv 3]
set prefix [lindex $argv 4]
set time_offset_ns [lindex $argv 5]
if {$time_offset_ns eq ""} { set time_offset_ns 0.0 }
set residue_start_ns [lindex $argv 6]
if {$residue_start_ns eq ""} { set residue_start_ns -1.0 }
set frame_stride [lindex $argv 7]
if {$frame_stride eq ""} { set frame_stride 50 }
set sampled_dt_ns [lindex $argv 8]
if {$sampled_dt_ns eq ""} { set sampled_dt_ns 0.5 }

mol new $psf type psf waitfor all
mol addfile $pdb type pdb waitfor all
# Production trajectory is saved every 10 ps. Analyze every 50th frame = 0.5 ns.
mol addfile $dcd type dcd first 0 last -1 step $frame_stride waitfor all

set dna [atomselect top "nucleic"]
set nframes [molinfo top get numframes]
set pi [expr {acos(-1.0)}]
set cos45 [expr {cos(45.0*$pi/180.0)}]

proc vsub {a b} {
    return [list [expr {[lindex $a 0]-[lindex $b 0]}] \
                 [expr {[lindex $a 1]-[lindex $b 1]}] \
                 [expr {[lindex $a 2]-[lindex $b 2]}]]
}
proc vdot {a b} {
    return [expr {[lindex $a 0]*[lindex $b 0]+[lindex $a 1]*[lindex $b 1]+[lindex $a 2]*[lindex $b 2]}]
}
proc vcross {a b} {
    return [list [expr {[lindex $a 1]*[lindex $b 2]-[lindex $a 2]*[lindex $b 1]}] \
                 [expr {[lindex $a 2]*[lindex $b 0]-[lindex $a 0]*[lindex $b 2]}] \
                 [expr {[lindex $a 0]*[lindex $b 1]-[lindex $a 1]*[lindex $b 0]}]]
}
proc vnorm {a} { return [expr {sqrt([vdot $a $a])}] }
proc distance {a b} { return [vnorm [vsub $a $b]] }

for {set r 1} {$r <= 40} {incr r} {
    set basesel($r) [atomselect top "nucleic and resid $r and not name P O1P O2P O5' C5' H5' H5'' C4' H4' O4' C1' H1' C2' H2' H2'' C3' H3' O3' H5T H3T"]
    if {$baseType eq "ADE"} {
        set planesel($r) [atomselect top "nucleic and resid $r and name N9 C4 C8"]
        set chisel($r) [atomselect top "nucleic and resid $r and name O4' C1' N9 C4"]
    } else {
        set planesel($r) [atomselect top "nucleic and resid $r and name N1 C2 C6"]
        set chisel($r) [atomselect top "nucleic and resid $r and name O4' C1' N1 C2"]
    }
    set sasa_sum($r) 0.0
    set chi_sum($r) 0.0
    set chi2_sum($r) 0.0
    set syn_count($r) 0
    set unstack_count($r) 0
}

set tf [open ${prefix}_base_flipping_time.csv w]
puts $tf "time_ns,mean_base_sasa_A2,stacked_adjacent_pairs,unstacked_base_fraction,syn_like_fraction"
set analyzed 0
for {set f 1} {$f < $nframes} {incr f} {
    $dna frame $f
    set time [expr {$time_offset_ns+($f-1)*$sampled_dt_ns}]
    set total_sasa 0.0
    set syn_total 0
    for {set r 1} {$r <= 40} {incr r} {
        $basesel($r) frame $f
        $planesel($r) frame $f
        $chisel($r) frame $f
        set sasa [measure sasa 1.4 $dna -restrict $basesel($r) -samples 200]
        set total_sasa [expr {$total_sasa+$sasa}]
        set xyz [$planesel($r) get {x y z}]
        set p1 [lindex $xyz 0]; set p2 [lindex $xyz 1]; set p3 [lindex $xyz 2]
        set normal($r) [vcross [vsub $p2 $p1] [vsub $p3 $p1]]
        set center($r) [measure center $basesel($r) weight mass]
        set ids [$chisel($r) get index]
        set chi [measure dihed $ids frame $f]
        if {$chi > 180.0} { set chi [expr {$chi-360.0}] }
        if {$time >= $residue_start_ns} {
            set sasa_sum($r) [expr {$sasa_sum($r)+$sasa}]
            set chi_sum($r) [expr {$chi_sum($r)+$chi}]
            set chi2_sum($r) [expr {$chi2_sum($r)+$chi*$chi}]
        }
        # Operational syn-like range: |chi| < 90 degrees.
        if {abs($chi) < 90.0} {
            if {$time >= $residue_start_ns} { incr syn_count($r) }
            incr syn_total
        }
    }
    set npairs 0
    for {set r 1} {$r < 40} {incr r} {
        set d [distance $center($r) $center([expr {$r+1}])]
        set nn [expr {[vnorm $normal($r)]*[vnorm $normal([expr {$r+1}])]}]
        set parallel 0.0
        if {$nn > 0} { set parallel [expr {abs([vdot $normal($r) $normal([expr {$r+1}])]/$nn)}] }
        set stacked($r) [expr {$d <= 5.5 && $parallel >= $cos45}]
        if {$stacked($r)} { incr npairs }
    }
    set nunstack 0
    for {set r 1} {$r <= 40} {incr r} {
        set left [expr {$r > 1 && $stacked([expr {$r-1}])}]
        set right [expr {$r < 40 && $stacked($r)}]
        if {!$left && !$right} {
            if {$time >= $residue_start_ns} { incr unstack_count($r) }
            incr nunstack
        }
    }
    puts $tf "$time,[expr {$total_sasa/40.0}],$npairs,[expr {$nunstack/40.0}],[expr {$syn_total/40.0}]"
    if {$time >= $residue_start_ns} { incr analyzed }
}
close $tf

set rf [open ${prefix}_base_flipping_residue.csv w]
puts $rf "resid,mean_sasa_A2,unstacked_fraction,mean_chi_deg,chi_sd_deg,syn_like_fraction"
for {set r 1} {$r <= 40} {incr r} {
    set mean [expr {$chi_sum($r)/$analyzed}]
    set var [expr {$chi2_sum($r)/$analyzed-$mean*$mean}]
    if {$var < 0} { set var 0 }
    puts $rf "$r,[expr {$sasa_sum($r)/$analyzed}],[expr {$unstack_count($r)/double($analyzed)}],$mean,[expr {sqrt($var)}],[expr {$syn_count($r)/double($analyzed)}]"
}
close $rf
quit
