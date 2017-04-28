# TODO TODO TODO redo all these trajectories with temperature 303.8K because I'm an idiot.
# It probably won't make much difference, but we should be rigorous.
# Open them, align them, and write them back.

# load the molecules

mol new output/diubiquitin/diubiquitin_2016-05-18_A/diubiquitin_lys_48_6283/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top LYS-48

mol new output/diubiquitin/diubiquitin_2016-05-18_A/diubiquitin_lys_63_11433/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top LYS-63

mol new output/diubiquitin/diubiquitin_2016-05-18_A/diubiquitin_met_1_8154/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top MET-1

# Other diubiquitin types:

mol new output/diubiquitin/diubiquitin_2017_03_08/diubiquitin_lys_11_97240/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top LYS-11

#mol new output/diubiquitin/diubiquitin_2017_03_08/diubiquitin_lys_27_27683/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
#mol new output/diubiquitin/diubiquitin_2017_04_10/diubiquitin_lys_27_1005/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
#mol new output/diubiquitin/diubiquitin_2017_04_10/diubiquitin_lys_27_25186/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
#mol rename top LYS-27
mol new output/diubiquitin/diubiquitin_2017_04_18/diubiquitin_lys_27_1256/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top LYS-27

mol new output/diubiquitin/diubiquitin_2017_03_08/diubiquitin_lys_29_27702/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top LYS-29

mol new output/diubiquitin/diubiquitin_2017_03_08/diubiquitin_lys_33_27721/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top LYS-33

mol new output/diubiquitin/diubiquitin_2017_03_08/diubiquitin_lys_6_97254/trajectory.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol rename top LYS-6

# Add data from other runs

#mol addfile output/diubiquitin/diubiquitin_2016-05-18_B/diubiquitin_lys_48_6297/trajectory.pdb type pdb first 0 last -1 step 1 waitfor all 0
#mol addfile output/diubiquitin/diubiquitin_2016-05-18_B/diubiquitin_lys_63_11908/trajectory.pdb type pdb first 0 last -1 step 1 waitfor all 1
#mol addfile output/diubiquitin/diubiquitin_2016-05-18_B/diubiquitin_met_1_8824/trajectory.pdb type pdb first 0 last -1 step 1 waitfor all 2

#mol addfile output/diubiquitin/diubiquitin_2016-02-25/diubiquitin_met_1_30787/trajectory.pdb type pdb first 0 last -1 step 1 waitfor all 2

#mol addfile output/diubiquitin/diubiquitin_2015-12-07/diubiquitin_lys_48_21091/trajectory.pdb type pdb first 0 last -1 step 1 waitfor all 0
#mol addfile output/diubiquitin/diubiquitin_2015-12-07/diubiquitin_lys_63_16552/trajectory.pdb type pdb first 0 last -1 step 1 waitfor all 1
#mol addfile output/diubiquitin/diubiquitin_2015-12-07/diubiquitin_met_1_11116/trajectory.pdb type pdb first 0 last -1 step 1 waitfor all 2

#mol addfile output/diubiquitin/diubiquitin_2017_02_03/diubiquitin_lys_48_6870/trajectory.pdb type pdb first 0 last -1 step 1 waitfor all 0

# procedure for restoring default representations after clustering deletes them

proc restore_reps {m} {
	set r [molinfo $m get numreps]
    
    # molecule bodies
	mol addrep $m
	mol modselect $r $m resid 1 to 73 or resid 77
	mol modstyle $r $m QuickSurf 1.000000 0.500000 1.000000 1.000000
	mol modcolor $r $m Chain
	mol modmaterial $r $m Opaque
    
    set r [expr $r + 1]

    # tails
	mol addrep $m
	mol modselect $r $m resid 73 to 77
	mol modstyle $r $m QuickSurf 1.000000 0.500000 1.000000 1.000000
	mol modcolor $r $m ColorID 12
	mol modmaterial $r $m Opaque
    
    set r [expr $r + 1]

    # hydrophobic patches
	mol addrep $m
	mol modselect $r $m resid 8 or resid 44 or resid 68 or resid 70
	mol modstyle $r $m QuickSurf 1.200000 0.500000 1.000000 1.000000
	mol modcolor $r $m ColorID 13
	mol modmaterial $r $m Opaque
    mol showrep $m $r 0
    
    set r [expr $r + 1]

    # only his-68
	mol addrep $m
	mol modselect $r $m resid 68
	mol modstyle $r $m QuickSurf 1.200000 0.500000 1.000000 1.000000
	mol modcolor $r $m ColorID 10
	mol modmaterial $r $m Opaque
    
    set r [expr $r + 1]

    # only val-70
	mol addrep $m
	mol modselect $r $m resid 70
	mol modstyle $r $m QuickSurf 1.200000 0.500000 1.000000 1.000000
	mol modcolor $r $m ColorID 3
	mol modmaterial $r $m Opaque
    
    set r [expr $r + 1]
    
    # N terminus
	mol addrep $m
	mol modselect $r $m index 0
	mol modstyle $r $m QuickSurf 1.200000 0.500000 1.000000 1.000000
	mol modcolor $r $m ColorID 19
	mol modmaterial $r $m Opaque
    
    set r [expr $r + 1]
    
    # C terminus
	mol addrep $m
	mol modselect $r $m index 152
	mol modstyle $r $m QuickSurf 1.200000 0.500000 1.000000 1.000000
	mol modcolor $r $m ColorID 29
	mol modmaterial $r $m Opaque
}

restore_reps 0
restore_reps 1
restore_reps 2
restore_reps 3
restore_reps 4
restore_reps 5
restore_reps 6
restore_reps 7

mol on 0
mol off 1
mol off 2
mol off 3
mol off 4
mol off 5
mol off 6
mol off 7
