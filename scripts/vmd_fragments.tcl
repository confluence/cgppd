# Add polyubiquitin directory locations

# diubiquitin

set ubq(diubqMet1) "best_diubiquitin/diubiquitin_met_1_8154"
set ubq(diubqLys6) "best_diubiquitin/diubiquitin_lys_6_97254"
set ubq(diubqLys11) "best_diubiquitin/diubiquitin_lys_11_97240"
set ubq(diubqLys27) "best_diubiquitin/diubiquitin_lys_27_1256"
set ubq(diubqLys29) "best_diubiquitin/diubiquitin_lys_29_27702"
set ubq(diubqLys33) "best_diubiquitin/diubiquitin_lys_33_27721"
set ubq(diubqLys48) "best_diubiquitin/diubiquitin_lys_48_6283"
set ubq(diubqLys63) "best_diubiquitin/diubiquitin_lys_63_11433"

# diubiquitin with flexible loop and longer tail

# tetraubiquitin

set ubq(tetraubqMet1) "tetraubiquitin/tetraubiquitin_met_1_14451"
set ubq(tetraubqLys6) "tetraubiquitin/tetraubiquitin_lys_6_93879"
set ubq(tetraubqLys48) "tetraubiquitin/tetraubiquitin_lys_48_98513"
set ubq(tetraubqLys63) "tetraubiquitin/tetraubiquitin_lys_63_104466"

# tetraubiquitin with flexible loop and longer tail

# octaubiquitin

# octaubiquitin with flexible loop and longer tail

# procedure for restoring default representations after clustering deletes them

proc restore_reps {m} {
	set r [molinfo $m get numreps]
    
    # molecule bodies
	mol addrep $m
	mol modselect $r $m all
	mol modstyle $r $m VDW 1.0 100.0
	mol modcolor $r $m Chain
	mol modmaterial $r $m Glass1
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile44
	mol addrep $m
	mol modselect $r $m resid 8 or resid 44 or resid 68 or resid 70
	mol modstyle $r $m VDW 1.0 100.0
	mol modcolor $r $m ColorID 0
	mol modmaterial $r $m Diffuse

    mol showrep $m $r 0
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile44
	mol addrep $m
	mol modselect $r $m resid 44
	mol modstyle $r $m VDW 1.1 100.0
	mol modcolor $r $m ColorID 0
	mol modmaterial $r $m Diffuse
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile36
	mol addrep $m
	mol modselect $r $m resid 36 or resid 71 or resid 73
	mol modstyle $r $m VDW 1.0 100.0
	mol modcolor $r $m ColorID 7
	mol modmaterial $r $m Diffuse

    mol showrep $m $r 0
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile36
	mol addrep $m
	mol modselect $r $m resid 36
	mol modstyle $r $m VDW 1.1 100.0
	mol modcolor $r $m ColorID 7
	mol modmaterial $r $m Diffuse
    
    set r [expr $r + 1]

    # Leu8
	mol addrep $m
	mol modselect $r $m resid 8
	mol modstyle $r $m VDW 1.1 100.0
	mol modcolor $r $m ColorID 10
	mol modmaterial $r $m Diffuse

    mol showrep $m $r 0
}

# load a molecule

proc load_ubq {name} {
    set path "$::env(HOME)/repos/cgppd/output/$::ubq($name)/trajectory.pdb"
    mol new $path type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
    mol rename top $name
    mol delrep 0 top
    restore_reps top
}
