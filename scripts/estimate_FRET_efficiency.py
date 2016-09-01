#!/usr/bin/env python
# -*- coding: utf-8 -*-

from plot_objects import DiubiquitinSimulationGroup
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate FRET efficiency from diubiquitin simulations produced by cgppd")
    parser.add_argument("dirs", help="Individual directories to process", nargs="+")
    parser.add_argument("-r", "--reference-length", help="Set R0 value", type=float, default=50.0)


    args = parser.parse_args()

    simulation_group = DiubiquitinSimulationGroup.from_dirs(args.dirs)
    R0 = args.reference_length
    
    for (name, sim) in simulation_group.sims:
        print "Simulation %s:" % name
        
        for i, cluster in enumerate(sim.clusters, 1):
            sample_E = []
            
            for sample in cluster.samples:
                R = sample.length
                E = 1.0 / (1.0 + (R / R0))**6
                sample_E.append(E)
                
            average_E = sum(sample_E) / len(sample_E)
            
            print "Cluster %d (%d samples): average E %.3f" % (i, len(sample_E), average_E)

