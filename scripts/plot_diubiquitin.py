#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import argparse
import matplotlib.pyplot as plt
from collections import defaultdict
from plot_objects import DiubiquitinSimulationGroup

class DiubiquitinPlots(DiubiquitinSimulationGroup):
    def _ordered_sims(self, ordering):
        if ordering == "name":
            return sorted(self.sims, key=lambda s: s[0])
        else:
            return self.sims
    
    def _plot_vs_time(self, measurement, args):
        plt.figure()
        
        rows = len(self.sims)

        for i, (name, sim) in enumerate(self._ordered_sims(args.order_by), 1):
            values = [getattr(s, measurement) for s in sim.samples]
            plt.subplot(rows,1,i)
            plt.plot(values)
            plt.title(name)
            plt.xlabel("Sample no.")
            plt.ylabel(u"%s (Å)" % measurement)  
            plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.97, wspace=0.2, hspace=0.7)      
            
    def plot_radius(self, args):
        self._plot_vs_time("radius", args)
            
    def plot_length(self, args):
        self._plot_vs_time("length", args)
            
    def plot_fret_efficiency(self, args):
        self._add_fret_efficiency(args)
        self._plot_vs_time("fret_efficiency", args)
        
    def _plot_aggregate_histogram(self, measurement, args, units=u"Å"):
        plt.figure()
        
        aggregate_values = defaultdict(list)

        for name, sim in self.sims:
            aggregate_values[name].extend([getattr(s, measurement) for s in sim.samples])
        
        rows = len(aggregate_values)
        
        lastplot = None

        for i, name in enumerate(sorted(aggregate_values.keys()), 1):
            values = aggregate_values[name]
            if lastplot is not None:
                lastplot = plt.subplot(rows,1,i, sharex=lastplot, sharey=lastplot)
            else:
                lastplot = plt.subplot(rows,1,i)
            plt.hist(values, bins=100)
            plt.title(name)
            unit_str = " (%s)" % units if units else ""
            plt.xlabel(u"%s%s" % (measurement, unit_str))
            plt.ylabel("No. of samples")
            plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.97, wspace=0.2, hspace=0.7)

    def plot_agg_hist_length(self, args):        
        self._plot_aggregate_histogram("length", args, units=None)

    def plot_agg_hist_radius(self, args):        
        self._plot_aggregate_histogram("radius", args, units=None)

    def _plot_histogram(self, measurement, args, units=u"Å"):
        plt.figure()
        
        rows = len(self.sims)
        
        lastplot = None

        for i, (name, sim) in enumerate(self._ordered_sims(args.order_by), 1):
            values = [getattr(s, measurement) for s in sim.samples]
            if lastplot is not None:
                lastplot = plt.subplot(rows,1,i, sharex=lastplot, sharey=lastplot)
            else:
                lastplot = plt.subplot(rows,1,i)
            plt.hist(values, bins=100)
            plt.title(name)
            unit_str = " (%s)" % units if units else ""
            plt.xlabel(u"%s%s" % (measurement, unit_str))
            plt.ylabel("No. of samples")
            plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.97, wspace=0.2, hspace=0.7)      
            
    def plot_hist_radius(self, args):
        self._plot_histogram("radius", args)
            
    def plot_hist_length(self, args):
        self._plot_histogram("length", args)
    
    def _plot_cluster_histogram(self, measurement, args, units=u"Å"):
        for name, sim in self.sims:
            for description, clusters in sim.cluster_sets.items():
                plt.figure()
                
                rows = len(clusters)
                cols = 1
                        
                for i, cluster in enumerate(clusters):
                    values = [getattr(s, measurement) for s in cluster.samples]
                    plt.subplot(rows, cols, i + 1)
                    plt.hist(values, bins=100)
                    plt.title("%s: %s" % (name, description))
                    unit_str = " (%s)" % units if units else ""
                    plt.xlabel(u"Cluster %d (%d/%d samples) %s%s" % (i + 1, len(values), len(sim.samples), measurement, unit_str))
                    plt.ylabel("No. of samples")
                    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.97, wspace=0.2, hspace=0.7)
                
    def plot_cluster_hist_radius(self, args):
        self._plot_cluster_histogram("radius", args)

    def plot_cluster_hist_length(self, args):
        self._plot_cluster_histogram("length", args)
        
    def _add_fret_efficiency(self, args):
        R0 = args.reference_length
        
        if not "fret_efficiency" in self.extra_properties:
			# It's hacktastic
			for (name, sim) in self.sims:
				for s in sim.samples:
					s.fret_efficiency = 1.0 / (1.0 + ((s.length + args.pad_length) / R0)**6)
					
			self.extra_properties.add("fret_efficiency")

    def plot_hist_fret_efficiency(self, args):
        self._add_fret_efficiency(args)
        
        self._plot_histogram("fret_efficiency", args, units=None)
        
    def plot_agg_hist_fret_efficiency(self, args):
        self._add_fret_efficiency(args)
        
        self._plot_aggregate_histogram("fret_efficiency", args, units=None)

    def plot_cluster_hist_fret_efficiency(self, args):
        self._add_fret_efficiency(args)
        
        self._plot_cluster_histogram("fret_efficiency", args, units=None)
        

        
        

PLOTS = tuple(n[5:] for n in DiubiquitinPlots.__dict__ if n.startswith("plot_"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("dirs", help="Individual directories to process", nargs="+")
    parser.add_argument("-p", "--plot", dest="plots", help="Type of plot", choices=PLOTS, action="append")
    parser.add_argument("-r", "--reference-length", help="Set R0 value", type=float, default=50.0)
    parser.add_argument("-l", "--pad-length", help="Add padding value to molecule length to simulate presence of a chromatophore pair", type=float, default=20.0)
    parser.add_argument("-o", "--order-by", help="Order simulation subplots. If no ordering is specified, the order of the directory parameters will be preserved.", choices=(None, 'name'), default=None)

    args = parser.parse_args()

    simulation_group = DiubiquitinPlots.from_dirs(args.dirs)

    for plot in args.plots:
        getattr(simulation_group, "plot_%s" % plot)(args)
    plt.show()
