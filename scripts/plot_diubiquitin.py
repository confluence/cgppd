#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import math
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from plot_objects import DiubiquitinSimulationGroup

class DiubiquitinPlots(DiubiquitinSimulationGroup):
    def _ordered_sims(self, ordering):
        if ordering == "name":
            def ubq_sort(s):
                p = re.match("(Lys|Met).*?(\d+)", s[0]).groups()
                return int(p[1])
                
            return sorted(self.sims, key=ubq_sort)
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

        if args.save_svg:
            fig = plt.gcf()
            fig.savefig("%s_vs_time.svg" % measurement, format='svg')
            
    def plot_radius(self, args):
        self._plot_vs_time("radius", args)
            
    def plot_length(self, args):
        self._plot_vs_time("length", args)
            
    def plot_potential(self, args):
        self._plot_vs_time("potential", args)
            
    def plot_fret_efficiency(self, args):
        self._add_fret_efficiency(args)
        self._plot_vs_time("fret_efficiency", args)
        
    def _plot_aggregate_histogram(self, measurement, args, units=u"Å"):
        unit_str = " (%s)" % units if units else ""
        
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
            plt.ylabel("No. of samples")
            plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.97, wspace=0.2, hspace=0.7)
        
        plt.xlabel(u"%s%s" % (measurement, unit_str))

        if args.save_svg:
            fig = plt.gcf()
            fig.savefig("%s_aggregate_histogram.svg" % measurement, format='svg')

    def plot_agg_hist_length(self, args):        
        self._plot_aggregate_histogram("length", args, units=None)

    def plot_agg_hist_radius(self, args):        
        self._plot_aggregate_histogram("radius", args, units=None)

    def _plot_histogram(self, measurement, args, units=u"Å"):
        unit_str = " (%s)" % units if units else ""
        
        plt.figure()
        
        rows = len(self.sims)
        
        lastplot = None
        
        roman = ('i','ii','iii','iv','v','vi','vii','viii','ix','x','xi','xii','xiii','xiv','xv','xvi')

        for i, (name, sim) in enumerate(self._ordered_sims(args.order_by), 1):
            values = [getattr(s, measurement) for s in sim.samples]
            if lastplot is not None:
                lastplot = plt.subplot(rows,1,i, sharex=lastplot, sharey=lastplot)
            else:
                lastplot = plt.subplot(rows,1,i)
            plt.hist(values, bins=100)
            name = re.sub(r'(Lys|Met)(\d+) diUb.*', r'\1\2-linked diubiquitin', name)
            plt.title("(%s) %s" % (roman[i-1], name))
            plt.ylabel("No. of samples")
            #plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.97, wspace=0.2, hspace=0.7)
            plt.xlabel(u"%s%s" % (measurement.replace("_", " ").replace("fret", "FRET"), unit_str))

        if args.save_svg:
            fig = plt.gcf()
            fig.set_size_inches(5, 10)
            fig.savefig("%s_histogram.svg" % measurement, format='svg')
            
    def plot_hist_radius(self, args):
        self._plot_histogram("radius", args)
            
    def plot_hist_length(self, args):
        self._plot_histogram("length", args)
            
    def plot_hist_potential(self, args):
        self._plot_histogram("potential", args)
    
    def _plot_cluster_histogram(self, measurement, args, units=u"Å"):
        unit_str = " (%s)" % units if units else ""
        
        for name, sim in self.sims:
            for description, clusters in sim.cluster_sets.items():
                plt.figure()
                
                rows = len(clusters)
                cols = 1
                
                largest_clusters = [c for c in clusters if 100.0*len(c.samples)/len(sim.samples) >= args.cluster_cutoff]
                        
                for i, cluster in enumerate(largest_clusters):
                    values = [getattr(s, measurement) for s in cluster.samples]
                    plt.subplot(rows, cols, i + 1)
                    plt.hist(values, bins=100)
                    plt.ylabel("No. of samples")
                    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.97, wspace=0.2, hspace=0.7)
                    
                    print len(values)
                    print len(sim.samples)
                    print 100.0*len(values)/len(sim.samples)
                
                    plt.xlabel(u"Cluster %d (%d%%) %s%s" % (i, 100.0*len(values)/len(sim.samples), measurement, unit_str))
                
                #plt.title(name)

                if args.save_svg:
                    fig = plt.gcf()
                    fig.savefig("%s_cluster_histogram_%s_%s.svg" % (measurement, name, description), format='svg')
                
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
    
    def plot_average_contacts(self, args):
        plt.figure()
        
        contact_averages_per_sim = [sim.cached_contacts for (name, sim) in self._ordered_sims(args.order_by)]
            
        chain_names = {
            "A": "Distal",
            "B": "Proximal",
        }
                        
        rows = len(contact_averages_per_sim)
        cols = len(contact_averages_per_sim[0])

        lastplot = None
        
        roman = ('i','ii','iii','iv','v','vi','vii','viii','ix','x','xi','xii','xiii','xiv','xv','xvi')

        for i, (name, sim) in enumerate(self._ordered_sims(args.order_by), 1):
            contact_averages = sim.cached_contacts
            linkage, _ = name.split()[:2]
            
            for j, (chain, averages) in enumerate(contact_averages, 1):
                if args.collapse_sidechain and j == 1:
                    averages = averages[:-2] + [sum(averages[-2:])/2]
                
                outliers = []
                if args.truncate_outliers < 100:
                    for ai, a in enumerate(averages):
                        if a > args.truncate_outliers:
                            averages[ai] = args.truncate_outliers
                            outliers.append((ai, a))
                
                residues = range(1, len(averages) + 1)
                
                if lastplot is not None:
                    lastplot = plt.subplot(rows, cols, (i - 1) * cols + j, sharex=lastplot, sharey=lastplot)
                else:
                    lastplot = plt.subplot(rows, cols, (i - 1) * cols + j, xticks=range(0,80,10))
                
                if j==2:
                    pos = int(re.search(r'\d+', linkage).group(0))
                    plt.vlines(pos, 0, args.truncate_outliers, colors='y', linestyles='dotted', zorder=0)
                    
                plt.bar(residues, averages)
                
                if outliers:
                    xoffsets = [0] * len(outliers)
                    xoffsets[0] -= 7
                    xoffsets[-1] += 7
                
                    for (ai, a), xoff in zip(outliers, xoffsets):
                        lastplot.annotate('%.1f' % a, (ai+1, args.truncate_outliers), xytext=(0 + xoff, 4), textcoords='offset points', size=7, ha='center')
                
                #if j == 1:
                    #lastplot.annotate("(%s) %s" % (roman[i-1], linkage), (0, 0), xytext=(-80, 10), textcoords='offset points', xycoords='axes fraction')
                
                if i == 1:
                    plt.title(chain_names[chain])
                
                if i == rows:
                    plt.xlabel("Residue")
                else:
                    plt.setp(lastplot.get_xticklabels(), visible=False)
                
                if i == math.ceil(rows/2) and j == 1:
                    plt.ylabel(u"Mean no. of contacts (cutoff: %g Å)" % args.contact_cutoff)
        
        #plt.subplots_adjust(wspace=0.1)
        #plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.97, wspace=0.2, hspace=0.7)

        if args.save_svg:
            fig = plt.gcf()
            fig.set_size_inches(6, 10)
            fig.savefig("average_contacts_%g_Å.svg" % args.contact_cutoff, format='svg')


PLOTS = tuple(n[5:] for n in DiubiquitinPlots.__dict__ if n.startswith("plot_"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("dirs", help="Individual directories to process", nargs="+")
    parser.add_argument("-p", "--plot", dest="plots", help="Type of plot", choices=PLOTS, action="append")
    parser.add_argument("-r", "--reference-length", help="Set R0 value", type=float, default=50.0)
    parser.add_argument("-l", "--pad-length", help="Add padding value to molecule length to simulate presence of a chromatophore pair", type=float, default=20.0)
    parser.add_argument("-c", "--contact-cutoff", help="In the contact averages plot, cutoff for determining whether two residues are in contact", type=float, default=7.0)
    parser.add_argument("-u", "--cluster-cutoff", help="Cutoff cluster size", type=float, default=20.0)
    parser.add_argument("-a", "--collapse-sidechain", help="In the contact averages plot, collapse the fake sidechain residue into the end residue of the distal tail by averaging their values", action="store_true", default=False)
    parser.add_argument("-t", "--truncate-outliers", help="In the contact averages plot, truncate values to this maximum", type=float, default=100.0)
    parser.add_argument("-o", "--order-by", help="Order simulation subplots. If no ordering is specified, the order of the directory parameters will be preserved.", choices=(None, 'name'), default=None)
    parser.add_argument("-s", "--save-svg", help="Save plots as SVG", action="store_true", default=False)

    args = parser.parse_args()

    simulation_group = DiubiquitinPlots.from_dirs(args)

    for plot in args.plots:
        getattr(simulation_group, "plot_%s" % plot)(args)
    plt.show()
