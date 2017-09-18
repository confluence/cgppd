# Add polyubiquitin directory locations

# diubiquitin

set ubq(diubq1) "best_diubiquitin/diubiquitin_met_1_8154"
set ubq(diubq6) "best_diubiquitin/diubiquitin_lys_6_97254"
set ubq(diubq11) "best_diubiquitin/diubiquitin_lys_11_97240"
set ubq(diubq27) "best_diubiquitin/diubiquitin_lys_27_1256"
set ubq(diubq29) "best_diubiquitin/diubiquitin_lys_29_27702"
set ubq(diubq33) "best_diubiquitin/diubiquitin_lys_33_27721"
set ubq(diubq48) "best_diubiquitin/diubiquitin_lys_48_6283"
set ubq(diubq63) "best_diubiquitin/diubiquitin_lys_63_11433"

# diubiquitin with flexible loop and longer tail

#set ubq(diubq1ll) "diubiquitin_ll/"
#set ubq(diubq6ll) "diubiquitin_ll/"
set ubq(diubq11ll) "diubiquitin_ll/diubiquitin_lys_11_longtail_loop_21139"
set ubq(diubq27ll) "diubiquitin_ll/diubiquitin_lys_27_longtail_loop_13597"
#set ubq(diubq29ll) "diubiquitin_ll/"
#set ubq(diubq33ll) "diubiquitin_ll/"
#set ubq(diubq48ll) "diubiquitin_ll/"
#set ubq(diubq63ll) "diubiquitin_ll/"

# tetraubiquitin

set ubq(tetraubq1) "tetraubiquitin/tetraubiquitin_met_1_14451"
set ubq(tetraubq6) "tetraubiquitin/tetraubiquitin_lys_6_93879"
set ubq(tetraubq48) "tetraubiquitin/tetraubiquitin_lys_48_98513"
set ubq(tetraubq63) "tetraubiquitin/tetraubiquitin_lys_63_104466"

# tetraubiquitin with flexible loop and longer tail

#set ubq(tetraubq1ll) "tetraubiquitin_ll/"
#set ubq(tetraubq6ll) "tetraubiquitin_ll/"
set ubq(tetraubq48ll) "tetraubiquitin_ll/tetraubiquitin_lys_48_longtail_loop_23137"
#set ubq(tetraubq63ll) "tetraubiquitin_ll/"

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
