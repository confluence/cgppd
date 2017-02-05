# load the molecules

mol new output/diubiquitin/diubiquitin_2016-05-18_A/diubiquitin_lys_48_6283/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top LYS-48

mol new output/diubiquitin/diubiquitin_2016-05-18_A/diubiquitin_lys_63_11433/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top LYS-63

mol new output/diubiquitin/diubiquitin_2016-05-18_A/diubiquitin_met_1_8154/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top MET-1

# Add data from other runs

mol addfile output/diubiquitin/diubiquitin_2016-05-18_B/diubiquitin_lys_48_6297/trajectory.pdb type pdb first 0 last -1 step 1 waitfor 1 0
mol addfile output/diubiquitin/diubiquitin_2016-05-18_B/diubiquitin_lys_63_11908/trajectory.pdb type pdb first 0 last -1 step 1 waitfor 1 1
mol addfile output/diubiquitin/diubiquitin_2016-05-18_B/diubiquitin_met_1_8824/trajectory.pdb type pdb first 0 last -1 step 1 waitfor 1 2

mol addfile output/diubiquitin/diubiquitin_2015-12-07/diubiquitin_lys_48_21091/trajectory.pdb type pdb first 0 last -1 step 1 waitfor 1 0
mol addfile output/diubiquitin/diubiquitin_2015-12-07/diubiquitin_lys_63_16552/trajectory.pdb type pdb first 0 last -1 step 1 waitfor 1 1
mol addfile output/diubiquitin/diubiquitin_2015-12-07/diubiquitin_met_1_11116/trajectory.pdb type pdb first 0 last -1 step 1 waitfor 1 2

# procedure for restoring default representations after clustering deletes them

proc restore_reps {m} {
	set r [molinfo $m get numreps]

	mol addrep $m
	mol modselect $r $m all
	mol modstyle $r $m VDW 1.000000 12.000000
	mol modcolor $r $m Chain
	mol modmaterial $r $m Glass1

	mol addrep $m
	mol modselect [expr $r + 1] $m resid 74 to 76
	mol modstyle [expr $r + 1] $m VDW 1.000000 12.000000
	mol modcolor [expr $r + 1] $m ColorID 12
	mol modmaterial [expr $r + 1] $m Opaque

	mol addrep $m
	mol modselect [expr $r + 2] $m resid 8 or resid 44 or resid 68 or resid 70
	mol modstyle [expr $r + 2] $m VDW 1.000000 12.000000
	mol modcolor [expr $r + 2] $m ColorID 13
	mol modmaterial [expr $r + 2] $m Opaque
}

