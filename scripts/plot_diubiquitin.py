#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
from plot_objects import Simulation

# TODO do this properly; subgraphs and labels, etc. Commandline params for plots.
def plot_length(simulations):
    for simulation in simulations:
        plt.hist([s.length for s in simulation.samples])
        plt.show()

        plt.hist([s.cluster.length for s in simulation.samples])
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("dirs", help="Individual directories to process", nargs="+")
    parser.add_argument("-c", "--cutoff", help="Use cluster with this cutoff")

    args = parser.parse_args()

    simulations = [Simulation.from_dir(d, args.cutoff) for d in args.dirs]
    plot_length(simulations)

