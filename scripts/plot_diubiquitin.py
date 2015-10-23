#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import argparse
import matplotlib.pyplot as plt
from plot_objects import Simulation

class DiubiquitinSimulationGroup(object):
    NAME = re.compile("diubiquitin_(lys|met)_(\d+)_.*")
    
    def __init__(self, sims):
        self.sims = sims
        
    @classmethod
    def from_dirs(cls, dirs, cutoff):
        sims = []
        
        for d in dirs:
            dirname = os.path.basename(d)
            name_match = cls.NAME.match(dirname)
            if name_match is None:
                sys.exit("'%s' does not look like a diubiquitin simulation." % dirname)
            res, index = name_match.groups()
            
            sims.append(("%s-%s" % (res.upper(), index), Simulation.from_dir(d, cutoff)))
            
        return cls(sims)

    def _plot_vs_time(self, measurement, args):
        rows = len(self.sims)

        for i, (name, sim) in enumerate(self.sims, 1):
            values = [getattr(s, measurement) for s in sim.samples]
            plt.subplot(rows,1,i)
            plt.plot(values)
            plt.title(name)
            plt.xlabel("Sample no.")
            plt.ylabel(u"%s (Å)" % measurement)        
            
    def plot_radius(self, args):
        self._plot_vs_time("radius", args)
            
    def plot_length(self, args):
        self._plot_vs_time("length", args)

    def _plot_histogram(self, measurement, args):
        rows = len(self.sims)
        
        for i, (name, sim) in enumerate(self.sims, 1):
            values = [getattr(s, measurement) for s in sim.samples]
            plt.subplot(rows,1,i)
            plt.hist(values)
            plt.title(name)
            plt.xlabel(u"%s (Å)" % measurement)
            plt.ylabel("No. of samples")      
            
    def plot_hist_radius(self, args):
        self._plot_histogram("radius", args)
            
    def plot_hist_length(self, args):
        self._plot_histogram("length", args)

    #some kind of meaningful cluster plot?
    #build the xvg plot into this?

PLOTS = tuple(n[5:] for n in DiubiquitinSimulationGroup.__dict__ if n.startswith("plot_"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("dirs", help="Individual directories to process", nargs="+")
    parser.add_argument("-c", "--cutoff", help="Use cluster with this cutoff")
    parser.add_argument("-p", "--plot", dest="plots", help="Type of plot", choices=PLOTS, action="append")

    args = parser.parse_args()

    simulation_group = DiubiquitinSimulationGroup.from_dirs(args.dirs, args.cutoff)

    for plot in args.plots:
        plt.figure()
        getattr(simulation_group, "plot_%s" % plot)(args)
    plt.show()
