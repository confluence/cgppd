# Add polyubiquitin directory locations

# diubiquitin

set ubq(diubq1) "polyubiquitin/diubiquitin_met_1_25542"
set ubq(diubq6) "polyubiquitin/diubiquitin_lys_6_13226"
set ubq(diubq11) "polyubiquitin/diubiquitin_lys_11_24746"
set ubq(diubq27) "polyubiquitin/diubiquitin_lys_27_7868"
set ubq(diubq29) "polyubiquitin/diubiquitin_lys_29_29440"
set ubq(diubq33) "polyubiquitin/diubiquitin_lys_33_16163"
set ubq(diubq48) "polyubiquitin/diubiquitin_lys_48_4168"
set ubq(diubq63) "polyubiquitin/diubiquitin_lys_63_31307"

# diubiquitin with flexible loop and longer tail

set ubq(diubq1ll) "polyubiquitin_ll/diubiquitin_met_1_longtail_loop_31998"
set ubq(diubq6ll) "polyubiquitin_ll/diubiquitin_lys_6_longtail_loop_1688"
set ubq(diubq11ll) "polyubiquitin_ll/diubiquitin_lys_11_longtail_loop_21139"
set ubq(diubq27ll) "polyubiquitin_ll/diubiquitin_lys_27_longtail_loop_13597"
set ubq(diubq29ll) "polyubiquitin_ll/diubiquitin_lys_29_longtail_loop_2190"
set ubq(diubq33ll) "polyubiquitin_ll/diubiquitin_lys_33_longtail_loop_9343"
set ubq(diubq48ll) "polyubiquitin_ll/diubiquitin_lys_48_longtail_loop_32531"
set ubq(diubq63ll) "polyubiquitin_ll/diubiquitin_lys_63_longtail_loop_32075"

# tetraubiquitin

set ubq(tetraubq1) "polyubiquitin/tetraubiquitin_met_1_14451"
set ubq(tetraubq6) "polyubiquitin/tetraubiquitin_lys_6_93879"
set ubq(tetraubq48) "polyubiquitin/tetraubiquitin_lys_48_98513"
set ubq(tetraubq63) "polyubiquitin/tetraubiquitin_lys_63_104466"

# tetraubiquitin with flexible loop and longer tail

set ubq(tetraubq1ll) "polyubiquitin_ll/tetraubiquitin_met_1_longtail_loop_31445"
set ubq(tetraubq6ll) "polyubiquitin_ll/tetraubiquitin_lys_6_longtail_loop_16256"
set ubq(tetraubq48ll) "polyubiquitin_ll/tetraubiquitin_lys_48_longtail_loop_23137"
set ubq(tetraubq63ll) "polyubiquitin_ll/tetraubiquitin_lys_63_longtail_loop_24659"

# octaubiquitin

set ubq(octaubq11) "polyubiquitin/octaubiquitin_lys_11_14908"

# octaubiquitin with flexible loop and longer tail

set ubq(octaubq11ll) "polyubiquitin_ll/octaubiquitin_lys_11_longtail_loop_27299"

# procedure for restoring default representations after clustering deletes them

proc restore_reps {m} {
	set r [molinfo $m get numreps]
    
    # molecule bodies
	mol addrep $m
	mol modselect $r $m all
	mol modstyle $r $m QuickSurf 1.0 0.5 1.0 1.0
	mol modcolor $r $m Chain
	mol modmaterial $r $m Translucent
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile44
	mol addrep $m
	mol modselect $r $m resid 8 or resid 44 or resid 68 or resid 70
	mol modstyle $r $m QuickSurf 1.5 0.5 1.0 1.0
	mol modcolor $r $m ColorID 0
	mol modmaterial $r $m Diffuse
    
    set r [expr $r + 1]

    # Ile44
	mol addrep $m
	mol modselect $r $m resid 44
	mol modstyle $r $m QuickSurf 1.5 0.5 1.0 1.0
	mol modcolor $r $m ColorID 0
	mol modmaterial $r $m Diffuse

    mol showrep $m $r 0
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile36
	mol addrep $m
	mol modselect $r $m resid 36 or resid 71 or resid 73
	mol modstyle $r $m QuickSurf 1.5 0.5 1.0 1.0
	mol modcolor $r $m ColorID 7
	mol modmaterial $r $m Diffuse
    
    set r [expr $r + 1]

    # Ile36
	mol addrep $m
	mol modselect $r $m resid 36
	mol modstyle $r $m QuickSurf 1.5 0.5 1.0 1.0
	mol modcolor $r $m ColorID 7
	mol modmaterial $r $m Diffuse

    mol showrep $m $r 0
    
    set r [expr $r + 1]

    # Leu8
	mol addrep $m
	mol modselect $r $m resid 8
	mol modstyle $r $m QuickSurf 1.5 0.5 1.0 1.0
	mol modcolor $r $m ColorID 10
	mol modmaterial $r $m Diffuse

    mol showrep $m $r 0
}

# Alignment

proc align {m seltext} {
    set ref [atomselect $m $seltext frame 0]
    set sel [atomselect $m $seltext]
    set all [atomselect $m all]
    set n [molinfo $m get numframes]

    for { set i 1 } { $i < $n } { incr i } {
        $sel frame $i
        $all frame $i
        $all move [measure fit $sel $ref]
    }
    $ref delete
    $all delete
    $sel delete
    return
}

# load a molecule

proc load_ubq {name} {
    set path "$::env(HOME)/repos/cgppd/output/$::ubq($name)/trajectory.pdb"
    mol new $path type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
    mol rename top $name
    mol delrep 0 top
    restore_reps top
    
    color Chain A yellow
    color Chain B orange
    
    if {[string match "diubq*" $name]} {
        align top "chain A"
    } elseif {[string match "tetraubq*" $name]} {
        align top "chain B or chain C"
        color Chain C orange2
        color Chain D red
    } elseif {[string match "octaubq*" $name]} {
        align top "chain C or chain D or chain E or chain F"
    }
    
    
    display resetview
}

# Set display

display projection Orthographic
color Display Background white
axes location Off

display ambientocclusion on
display aoambient 0.8
display aodirect 0.3
display shadows on

#                               

# comparison structures

set ubq_ref(diubq1) { "2w9n-AB.pdb" "3axc-AB.pdb" }
set ubq_ref(diubq6) { "2xk5.pdb" }
set ubq_ref(diubq11) { "2xew-AB.pdb" "3nob-AB.pdb" }
set ubq_ref(diubq27) { }
set ubq_ref(diubq29) { "4s22.pdb" "4s1z.pdb" }
set ubq_ref(diubq33) { "5af4.pdb" "4xyz.pdb" "5af6.pdb" }
set ubq_ref(diubq48) { "1aar.pdb" "2pe9.pdb" "3aul.pdb" }
set ubq_ref(diubq63) { "2jf5-AB.pdb" "3a1q-AB.pdb" "3a1q-DE.pdb" }

# frames for clusters

# 68%
set cluster_frame(diubq1) { 0 }
# 53%, 12%
set cluster_frame(diubq6) { 15 50 }
# 49%
set cluster_frame(diubq11) { 15 }
set cluster_frame(diubq27) {}
set cluster_frame(diubq29) {}
set cluster_frame(diubq33) {}
# 67%
set cluster_frame(diubq48) { 5 }
# 77%
set cluster_frame(diubq63) { 0 }

# 25%
set cluster_frame(diubq1ll) { 20 }
#42%, 11%, 10%
set cluster_frame(diubq6ll) { 10 25 5 }
# 18%, 10%, 10%
set cluster_frame(diubq11ll) { 10 16 54 }
set cluster_frame(diubq27ll) {}
set cluster_frame(diubq29ll) {}
set cluster_frame(diubq33ll) {}
# 40%, 28%
set cluster_frame(diubq48ll) { 5 0 }
# 19%, 10%
set cluster_frame(diubq63ll) { 0 20 }

# representation for comparison structures

proc comparison_reps {m} {
    mol delrep 0 $m

	set r [molinfo $m get numreps]
    
    # molecule bodies
	mol addrep $m
	mol modselect $r $m name CA
	mol modstyle $r $m QuickSurf 1.0 0.5 1.0 1.0
	mol modcolor $r $m Chain
	mol modmaterial $r $m Translucent
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile44
	mol addrep $m
	mol modselect $r $m (resid 8 or resid 44 or resid 68 or resid 70) and name CA
	mol modstyle $r $m QuickSurf 1.5 0.5 1.0 1.0
	mol modcolor $r $m ColorID 0
	mol modmaterial $r $m Diffuse
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile36
	mol addrep $m
	mol modselect $r $m (resid 36 or resid 71 or resid 73) and name CA
	mol modstyle $r $m QuickSurf 1.5 0.5 1.0 1.0
	mol modcolor $r $m ColorID 7
	mol modmaterial $r $m Diffuse

}

# align (WIP)

proc align_ref {mref f m seltext} {
    set ref [atomselect $mref $seltext frame $f]
    set sel [atomselect $m $seltext]
    set all [atomselect $m all]
    set n [molinfo $m get numframes]

    for { set i 0 } { $i < $n } { incr i } {
        $sel frame $i
        $all frame $i
        $all move [measure fit $sel $ref]
    }
    $ref delete
    $all delete
    $sel delete
    return
}

# load reference ubq

proc load_ubq_ref {name} {
    
    foreach pdbfile $::ubq_ref($name) {
        set path "$::env(HOME)/repos/cgppd/data/ubiquitin_comparisons/$pdbfile"
        mol new $path type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
        mol rename top "$name-$pdbfile"
        mol delrep 0 top
        comparison_reps top
    }
        
    color Chain A yellow
    color Chain B orange
        
    display resetview
}

# align these molecules to the given reference molecule

proc align_all {} {
    set all_mols [molinfo list]
    set mref [lindex $all_mols 0]
    set mols [lrange $all_mols 1 end]

    foreach m $mols {
        align_ref $mref 0 $m "resid 1 to 72 and name CA and chain A and (altloc \"\" or altloc \"A\")" 
    }
    
    display resetview
}

# Output files

proc load_diubq {linkage} {

    load_ubq "diubq${linkage}"
    load_ubq "diubq${linkage}ll"
    load_ubq_ref "diubq${linkage}"

    align_all
}

proc output_ubq {prefix} {
    
    set all_mols [molinfo list]
    
    foreach i $all_mols {
        foreach j $all_mols {
            mol off $j
        }
        
        mol on $i
        
        set name [molinfo $i get name]
        regexp "^diubq\[0-9\]+" $name dirname
        set path "$::env(HOME)/repos/cgppd/vmd_exports/$dirname"
        
        if { ![string match "*pdb" $name] } {
            set c 1
            foreach f $::cluster_frame($name) {
                animate goto $f
                render Tachyon $path/$prefix-$name-C$c-$f.dat /usr/local/lib/vmd/tachyon_LINUXAMD64 -aasamples 12 %s -format TARGA -o %s.tga
                set c [expr $c + 1]
            }
        } else {
            render Tachyon $path/$prefix-$name.dat /usr/local/lib/vmd/tachyon_LINUXAMD64 -aasamples 12 %s -format TARGA -o %s.tga
        }
    }

}

# TODO: We need to output 3 views of all molecules aligned *in the same way*. Use one dummy one-frame molecule to load always, and ignore it? Pick one which has a good view of the patches on chain A. Use chain A only?
# TODO: We should also save this cluster data so that we don't have to keep redoing it.

output_all_ubq {
# for each linkage
    set linkages { 1 6 11 27 29 33 48 63 }
    
    foreach l $linkages {
        # TODO load reference
        load_diubq $l
        scale by 0.8
        
        output_ubq ONE
        
        rotate x by 90
        output_ubq TWO
        
        rotate x by 270
        rotate y by 90
        output_ubq THREE
        
        # TODO delete all molecules
    }
}
