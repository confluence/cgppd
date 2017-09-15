# TODO TODO TODO redo all these trajectories with temperature 303.8K because I'm an idiot.
# It probably won't make much difference, but we should be rigorous.
# Open them, align them, and write them back.

# diubiquitin

mol new output/best_diubiquitin/diubiquitin_lys_48_6283/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top Lys48-diubq

mol new output/best_diubiquitin/diubiquitin_lys_63_11433/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top Lys63-diubq

mol new output/best_diubiquitin/diubiquitin_met_1_8154/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top Met1-diubq

mol new output/best_diubiquitin/diubiquitin_lys_11_97240/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top Lys11-diubq

mol new output/best_diubiquitin/diubiquitin_lys_27_1256/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top Lys27-diubq

mol new output/best_diubiquitin/diubiquitin_lys_29_27702/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top Lys29-diubq

mol new output/best_diubiquitin/diubiquitin_lys_33_27721/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top Lys33-diubq

mol new output/best_diubiquitin/diubiquitin_lys_6_97254/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top Lys6-diubq

# tetraubiquitin

# octaubiquitin

# procedure for restoring default representations after clustering deletes them

proc restore_reps {m} {
	set r [molinfo $m get numreps]
    
    # molecule bodies
	mol addrep $m
	mol modselect $r $m all
	mol modstyle $r $m VDW 0.5 100.0
	mol modcolor $r $m Chain
	mol modmaterial $r $m Glass1
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile44
	mol addrep $m
	mol modselect $r $m resid 8 or resid 44 or resid 68 or resid 70
	mol modstyle $r $m VDW 1.0 100.0
	mol modcolor $r $m ColorID 0
	mol modmaterial $r $m Opaque

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
	mol modmaterial $r $m Opaque

    mol showrep $m $r 0
    
    set r [expr $r + 1]

    # hydrophobic patch around Ile36
	mol addrep $m
	mol modselect $r $m resid 36
	mol modstyle $r $m VDW 1.1 100.0
	mol modcolor $r $m ColorID 7
	mol modmaterial $r $m Diffuse
}

