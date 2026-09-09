package require pbctools

set psf [lindex $argv 0]
set pdb [lindex $argv 1]
set dcd [lindex $argv 2]
set outprefix [lindex $argv 3]
set selection_start_ns [lindex $argv 4]
if {$selection_start_ns eq ""} { set selection_start_ns 30.0 }
set time_offset_ns [lindex $argv 5]
if {$time_offset_ns eq ""} { set time_offset_ns 0.0 }

mol new $psf type psf waitfor all
mol addfile $pdb type pdb waitfor all
# 5000 stored frames -> analyze every tenth frame (100 ps spacing).
mol addfile $dcd type dcd first 0 last -1 step 10 waitfor all
pbc unwrap -all -sel "nucleic"

set dna [atomselect top "nucleic and noh"]
set first [atomselect top "nucleic and resid 1 and name C1'"]
set last [atomselect top "nucleic and resid 40 and name C1'"]
set n [molinfo top get numframes]
set fh [open ${outprefix}_metrics.csv w]
puts $fh "frame,time_ns,rgyr_A,end_to_end_A"
set records {}
for {set f 1} {$f < $n} {incr f} {
    $dna frame $f
    $first frame $f
    $last frame $f
    set rg [measure rgyr $dna weight mass]
    set a [lindex [$first get {x y z}] 0]
    set b [lindex [$last get {x y z}] 0]
    set dx [expr {[lindex $a 0]-[lindex $b 0]}]
    set dy [expr {[lindex $a 1]-[lindex $b 1]}]
    set dz [expr {[lindex $a 2]-[lindex $b 2]}]
    set ee [expr {sqrt($dx*$dx+$dy*$dy+$dz*$dz)}]
    set t [expr {$time_offset_ns+($f-1)*0.1}]
    puts $fh "$f,$t,$rg,$ee"
    if {$t >= $selection_start_ns} { lappend records [list $rg $f $t $ee] }
}
close $fh

# Pick the frame closest to the median Rg in the final 20 ns.
set sorted [lsort -real -index 0 $records]
set median [lindex $sorted [expr {[llength $sorted]/2}]]
set target [lindex $median 0]
set best {}
set bestdiff 1e30
foreach rec $records {
    set diff [expr {abs([lindex $rec 0]-$target)}]
    if {$diff < $bestdiff} { set bestdiff $diff; set best $rec }
}
set bf [lindex $best 1]
set all [atomselect top all frame $bf]
$all writepdb ${outprefix}_representative_full.pdb
set onlydna [atomselect top "nucleic" frame $bf]
$onlydna writepdb ${outprefix}_representative_dna.pdb
set sf [open ${outprefix}_selected.txt w]
puts $sf "selected_frame $bf"
puts $sf "time_ns [lindex $best 2]"
puts $sf "rgyr_A [lindex $best 0]"
puts $sf "end_to_end_A [lindex $best 3]"
close $sf
quit
