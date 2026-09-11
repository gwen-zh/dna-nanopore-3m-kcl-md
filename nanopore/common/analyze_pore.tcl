package require pbctools

set psf [lindex $argv 0]
set pdb [lindex $argv 1]
set dcd [lindex $argv 2]
set prefix [lindex $argv 3]

mol new $psf type psf waitfor all
mol addfile $pdb type pdb waitfor all
mol addfile $dcd type dcd first 0 last -1 step 1 waitfor all
pbc unwrap -all -sel "segname AN1"

set dna [atomselect top "segname AN1 and noh"]
set lead [atomselect top "segname AN1 and resid 1 and name C1'"]
set tail [atomselect top "segname AN1 and resid 40 and name C1'"]
set protein [atomselect top "protein and noh"]
set n [molinfo top get numframes]
set out [open ${prefix}_pore_metrics.csv w]
puts $out "frame,time_ns,rgyr_A,end_to_end_A,lead_z_A,dna_com_z_A,dna_min_z_A,dna_max_z_A,dna_protein_contacts_3p5A,dna_contact_residues"

for {set f 1} {$f < $n} {incr f} {
    foreach sel [list $dna $lead $tail $protein] { $sel frame $f }
    set rg [measure rgyr $dna weight mass]
    set a [lindex [$lead get {x y z}] 0]
    set b [lindex [$tail get {x y z}] 0]
    set ee [veclength [vecsub $a $b]]
    set com [measure center $dna weight mass]
    set bounds [measure minmax $dna]
    set pairs [measure contacts 3.5 $dna $protein]
    set dna_indices [lsort -unique [lindex $pairs 0]]
    set contact_resids {}
    foreach idx $dna_indices {
        set one [atomselect top "index $idx" frame $f]
        lappend contact_resids [lindex [$one get resid] 0]
        $one delete
    }
    set contact_resids [lsort -integer -unique $contact_resids]
    set t [expr {($f-1)*0.025}]
    puts $out "$f,$t,$rg,$ee,[lindex $a 2],[lindex $com 2],[lindex [lindex $bounds 0] 2],[lindex [lindex $bounds 1] 2],[llength [lindex $pairs 0]],[join $contact_resids {;}]"
}
close $out
quit
