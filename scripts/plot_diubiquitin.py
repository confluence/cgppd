#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import argparse
import matplotlib.pyplot as plt
from plot_objects import DiubiquitinSimulationGroup

class DiubiquitinPlots(DiubiquitinSimulationGroup):
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

    def _plot_histogram(self, measurement, args, xlim=100, ylim=250, units=u"Å"):
        rows = len(self.sims)

        for i, (name, sim) in enumerate(self.sims, 1):
            values = [getattr(s, measurement) for s in sim.samples]
            plt.subplot(rows,1,i)
            plt.hist(values, bins=100)
            plt.title(name)
            unit_str = " (%s)" % units if units else ""
            plt.xlabel(u"%s%s" % (measurement, unit_str))
            plt.ylabel("No. of samples")
            plt.xlim([0, xlim])
            plt.ylim([0, ylim])
            plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.97, wspace=0.2, hspace=0.7)      
            
    def plot_hist_radius(self, args):
        self._plot_histogram("radius", args)
            
    def plot_hist_length(self, args):
        self._plot_histogram("length", args)

    #some kind of meaningful cluster plot?
    
    def _plot_cluster_histogram(self, measurement, args, xlim=100, ylim=1000, units=u"Å"):
        rows = len(self.sims)
        cols = max([len(s.clusters) for (n, s) in self.sims])
        
        for i, (name, sim) in enumerate(self.sims):
            for j, cluster in enumerate(sim.clusters):
                values = [getattr(s, measurement) for s in cluster.samples]
                
                subplot_no = i * cols + j + 1
                plt.subplot(rows,cols,subplot_no)
                plt.hist(values, bins=100)
                plt.title(name)
                unit_str = " (%s)" % units if units else ""
                plt.xlabel(u"Cluster %d %s%s" % (j + 1, measurement, unit_str))
                plt.ylabel("No. of samples")
                plt.xlim([0, xlim])
                plt.ylim([0, ylim])
                plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.97, wspace=0.2, hspace=0.7)
                
    def plot_cluster_hist_radius(self, args):
        self._plot_cluster_histogram("radius", args)

    def plot_cluster_hist_length(self, args):
        self._plot_cluster_histogram("length", args)
        
    def _add_fret_efficiency(self, args):
        R0 = args.reference_length
        
        # It's hacktastic
        for (name, sim) in self.sims:
            for s in sim.samples:
                s.fret_efficiency = 1.0 / (1.0 + ((s.length + args.pad_length) / R0)**6)

    def plot_hist_fret_efficiency(self, args):
        self._add_fret_efficiency(args)
        
        # TODO: automatic limits
        self._plot_histogram("fret_efficiency", args, xlim=1, units=None)

    def plot_cluster_hist_fret_efficiency(self, args):
        self._add_fret_efficiency(args)
        
        # TODO: automatic limits
        self._plot_cluster_histogram("fret_efficiency", args, xlim=1, units=None)
        

        
        

PLOTS = tuple(n[5:] for n in DiubiquitinPlots.__dict__ if n.startswith("plot_"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("dirs", help="Individual directories to process", nargs="+")
    parser.add_argument("-p", "--plot", dest="plots", help="Type of plot", choices=PLOTS, action="append")
    parser.add_argument("-r", "--reference-length", help="Set R0 value", type=float, default=50.0)
    parser.add_argument("-l", "--pad-length", help="Add padding value to molecule length to simulate presence of a chromatophore pair", type=float, default=20.0)

    args = parser.parse_args()

    simulation_group = DiubiquitinPlots.from_dirs(args.dirs)

    for plot in args.plots:
        plt.figure()
        getattr(simulation_group, "plot_%s" % plot)(args)
    plt.show()
