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

    mol showrep $m $r 0
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile44
	mol addrep $m
	mol modselect $r $m resid 44
	mol modstyle $r $m QuickSurf 1.5 0.5 1.0 1.0
	mol modcolor $r $m ColorID 0
	mol modmaterial $r $m Diffuse
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile36
	mol addrep $m
	mol modselect $r $m resid 36 or resid 71 or resid 73
	mol modstyle $r $m QuickSurf 1.5 0.5 1.0 1.0
	mol modcolor $r $m ColorID 7
	mol modmaterial $r $m Diffuse

    mol showrep $m $r 0
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile36
	mol addrep $m
	mol modselect $r $m resid 36
	mol modstyle $r $m QuickSurf 1.5 0.5 1.0 1.0
	mol modcolor $r $m ColorID 7
	mol modmaterial $r $m Diffuse
    
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
    if {[string match "diubq*" $name]} {
        align top "chain A"
    } elseif {[string match "tetraubq*" $name]} {
        align top "chain B or chain C"
    } elseif {[string match "octaubq*" $name]} {
        align top "chain C or chain D or chain E or chain F"
    }
}

# Set display

color Display Background white
axes location Off

# TODO turn on occlusion lighting and shadows
# TODO change chain colours to yellows and oranges
