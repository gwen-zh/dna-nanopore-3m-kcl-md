# Headless, CPU-only rendering of actual extracted trajectory frames.
if {[catch {
    lassign $argv input output centerx centery centerz span first last
    file mkdir $output
    mol new [file join $input visualization.psf] type psf waitfor all
    mol addfile [file join $input visualization.dcd] type dcd waitfor all
    set mid [molinfo top]
    if {[molinfo top get numframes] <= $last} { error "Requested missing frame" }
    mol delrep 0 $mid
    display projection Orthographic
    display resize 800 800
    display depthcue off
    display shadows off
    display ambientocclusion off
    display aoambient 0.75
    display aodirect 0.25
    if {[catch {render aasamples TachyonInternal 4} dt_aa_message]} {
        puts "AA_DEFAULT: $dt_aa_message"
    } else {
        puts "AA_SAMPLES: [render aasamples TachyonInternal]"
    }
    axes location Off
    color Display Background white
    color change rgb 0 0.32 0.55 0.71
    color change rgb 1 0.84 0.18 0.16
    color change rgb 3 0.95 0.56 0.12
    color change rgb 4 0.19 0.65 0.39
    color change rgb 6 0.63 0.65 0.68
    material add PoreFront
    material change opacity PoreFront 0.22
    material change ambient PoreFront 0.30
    material change diffuse PoreFront 0.65
    material change shininess PoreFront 0.35
    set p [atomselect top protein frame 0]
    set front {}
    set back {}
    foreach segment [lsort -unique [$p get segname]] {
        set s [atomselect top "protein and segname $segment" frame 0]
        set y [lindex [measure center $s] 1]
        if {$y < $centery} { lappend front $segment } else { lappend back $segment }
        $s delete
    }
    $p delete
    foreach group [list $back $front] mat {Opaque PoreFront} {
        mol representation NewCartoon 0.35 6 4.1 0
        mol color ColorID 0
        mol selection "protein and segname [join $group { }]"
        mol material $mat
        mol addrep $mid
    }
    mol representation VDW 0.25 10
    mol color ColorID 6
    mol selection "resname POPC and name P"
    mol material Opaque
    mol addrep $mid
    mol representation Licorice 0.24 12 12
    mol color ColorID 3
    mol selection "segname AN1 and noh"
    mol material Opaque
    mol addrep $mid
    foreach resid {1 40} colorid {1 4} {
        mol representation VDW 0.7 16
        mol color ColorID $colorid
        mol selection "segname AN1 and resid $resid and name C1'"
        mol material Opaque
        mol addrep $mid
    }
    mol representation Licorice 0.26 10 10
    mol color ColorID 10
    mol selection "protein and within 3.5 of segname AN1"
    mol material Opaque
    mol addrep $mid
    mol selupdate 6 $mid on
    mol on $mid
    display update on
    display resetview
    molinfo $mid set center_matrix [list [transoffset [list [expr {-$centerx}] [expr {-$centery}] [expr {-$centerz}]]]]
    rotate x by -90
    scale to [expr {2.5/$span}]
    puts "CAMERA [molinfo $mid get {center_matrix rotate_matrix scale_matrix global_matrix}]"
    display update
    set qc [open [file join $input [format "render_checks_%d_%d.csv" $first $last]] w]
    puts $qc "frame,rg_A,end_to_end_A,contact_pairs_3p5A"
    foreach frame [lsort -unique [list $first $last]] {
        set d [atomselect top "segname AN1 and noh" frame $frame]
        set p [atomselect top "protein and noh" frame $frame]
        set a [atomselect top "segname AN1 and resid 1 and name C1'" frame $frame]
        set b [atomselect top "segname AN1 and resid 40 and name C1'" frame $frame]
        set ee [veclength [vecsub [lindex [$a get {x y z}] 0] [lindex [$b get {x y z}] 0]]]
        set pairs [llength [lindex [measure contacts 3.5 $d $p] 0]]
        puts $qc "$frame,[measure rgyr $d weight mass],$ee,$pairs"
        foreach s [list $d $p $a $b] { $s delete }
    }
    close $qc
    for {set frame $first} {$frame <= $last} {incr frame} {
        animate goto $frame
        display update
        set target [file join $output [format "frame_%04d.tga" $frame]]
        render TachyonInternal $target
        puts "RENDERED $frame $target"
    }
    puts "RENDER_COMPLETE"
} message]} {
    puts stderr "RENDER_FAILED: $message"
    puts stderr $::errorInfo
    exit 1
}
quit
