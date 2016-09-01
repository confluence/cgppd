#!/usr/bin/env python
# -*- coding: utf-8 -*-

from plot_objects import DiubiquitinSimulationGroup
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate FRET efficiency from diubiquitin simulations produced by cgppd")
    parser.add_argument("dirs", help="Individual directories to process", nargs="+")

    args = parser.parse_args()

    simulation_group = DiubiquitinSimulationGroup.from_dirs(args.dirs)
    
    for (name, sim) in simulation_group.sims:
        print name

