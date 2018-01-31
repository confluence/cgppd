#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import optimize
import sys

from plot_objects import PolyalanineSimulationSequence

class PolyalaninePlots(PolyalanineSimulationSequence):

    def _plot_vs_n(self, measurement, args):
        xvalues, sims = zip(*self.sims)
        values = [np.sqrt(np.mean([getattr(s, measurement)**2 for s in sim.samples])) for sim in sims]

        plt.plot(xvalues, values, 'bo')
        
        if args.lj == "off":
            power = 1.0/2.0
        elif args.lj == "repulsive":
            power = 2.0/3.0

        fitfunc = lambda p, N: p[0] * np.array(N) ** power # Target function
        errfunc = lambda p, N, y: fitfunc(p, N) - y # Distance to the target function
        p0 = [1.0] # Initial guess for the parameters
        p1, success = optimize.leastsq(errfunc, p0[:], args=(xvalues, values))

        #logging.debug("A = %g" % p1[0])

        plt.plot(xvalues, [fitfunc(p1, x) for x in xvalues], 'b-')
        
        #if args.lj == "off":
            #plt.title("Polyalanine with LJ potentials off")
        #elif args.lj == "repulsive":
            #plt.title("Polyalanine with repulsive LJ potentials")

        plt.xlabel("Number of residues")
        plt.ylabel(u"Mean %s (Å)" % measurement)

        plt.xscale('log', basex=2) # TODO: investigate using the loglog function instead; maybe add an option for it
        plt.yscale('log')
        
        if args.save_svg:
            fig = plt.gcf()
            fig.savefig("%s_vs_N_LJ_%s.svg" % (measurement, args.lj), format='svg')
            
    def plot_mean_radius(self, args):
        self._plot_vs_n("radius", args)

    def plot_mean_length(self, args):
        self._plot_vs_n("length", args)
        
    def _rows_and_cols(self):
        rows = np.floor(np.sqrt(len(self.sims)))
        cols = np.ceil(len(self.sims)/rows)
        return int(rows), int(cols)

    def _plot_vs_time(self, measurement, args):
        rows, cols = self._rows_and_cols()

        for i in range(rows):
            for j in range(cols):
                if i * cols + j >= len(self.sims):
                    break

                n, sim = self.sims[i * cols + j]
                values = [getattr(s, measurement) for s in sim.samples]
                plt.subplot(rows, cols, i * cols + j + 1)
                plt.plot(values)
                plt.title("%d residues" % n)
                plt.xlabel("Sample no.")
                plt.ylabel(u"%s (Å)" % measurement)
                plt.xticks(rotation='vertical')
        
        plt.subplots_adjust(left=0.06, bottom=0.08, right=0.99, top=0.94, wspace=0.3, hspace=0.4)
                
        if args.lj == "off":
            plt.suptitle("Polyalanine with LJ potentials off")
        elif args.lj == "repulsive":
            plt.suptitle("Polyalanine with repulsive LJ potentials")
            
        if args.save_svg:
            fig = plt.gcf()
            fig.savefig("%s_vs_time_LJ_%s.svg" % (measurement, args.lj), format='svg')
            
    def plot_radius(self, args):
        self._plot_vs_time("radius", args)
            
    def plot_length(self, args):
        self._plot_vs_time("length", args)

    def _plot_histogram(self, measurement, args):
        rows, cols = self._rows_and_cols()

        for i in range(rows):
            for j in range(cols):
                if i * cols + j >= len(self.sims):
                    break
                
                n, sim = self.sims[i * cols + j]
                values = [getattr(s, measurement) for s in sim.samples]
                plt.subplot(rows, cols, i * cols + j + 1)
                plt.hist(values, bins=100)
                plt.title("%d residues" % n)
                plt.xlabel(u"%s (Å)" % measurement)
                plt.ylabel("No. of samples") 
                plt.xticks(rotation='vertical')
        
        plt.subplots_adjust(left=0.06, bottom=0.08, right=0.99, top=0.94, wspace=0.3, hspace=0.4)
                
        if args.lj == "off":
            plt.suptitle("Polyalanine with LJ potentials off")
        elif args.lj == "repulsive":
            plt.suptitle("Polyalanine with repulsive LJ potentials")
        
        if args.save_svg:
            fig = plt.gcf()
            fig.savefig("%s_histogram_LJ_%s.svg" % (measurement, args.lj), format='svg')
            
    def plot_hist_radius(self, args):
        self._plot_histogram("radius", args)
            
    def plot_hist_length(self, args):
        self._plot_histogram("length", args)
        
    
PLOTS = tuple(n[5:] for n in PolyalaninePlots.__dict__ if n.startswith("plot_"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    
    parser.add_argument("lj", help="Type of Lennard-Jones interaction", choices=("off", "repulsive"))
    parser.add_argument("dirs", help="Individual directories to process", nargs="+")
    
    parser.add_argument("-p", "--plot", dest="plots", help="Type of plot", choices=PLOTS, action="append")
    parser.add_argument("-s", "--save-svg", help="Save plots as SVG", action="store_true", default=False)

    args = parser.parse_args()

    simulation_set = PolyalaninePlots.from_dirs(args)
    
    for plot in args.plots:
        plt.figure()
        getattr(simulation_set, "plot_%s" % plot)(args)
    plt.show()
