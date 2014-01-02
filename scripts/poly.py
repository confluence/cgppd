#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import math
import subprocess

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate .rib and possibly .pdb files of polypeptides of various lengths")

    parser.add_argument("residue", help="The residue type")
    parser.add_argument("end", type=int, help="The end of the range")
    parser.add_argument("-s", "--start", type=int, help="The start of the range")

    parser.add_argument("-r", "--range", default="exp", help="The type of range: lin(ear) or exp(onential)")
    parser.add_argument("-b", "--base", default=2, type=int, help="A base for the exponential range")

    parser.add_argument("--phi", default=-120.0, type=float, help="The default phi value")
    parser.add_argument("--psi", default=120.0, type=float, help="The default psi value")

    parser.add_argument("-o", "--output-dir", default=".", help="The output directory")
    parser.add_argument("-p", "--path-to-ribosome-dir", help="Path to directory containing ribosome executable and res.zmat file")

    args = parser.parse_args()

    if args.range == 'lin':
        start = args.start or 4
        rng = range(start, args.end + 1)
    else:
        start = args.start or int(math.ceil(math.log(4, args.base)))
        rng = (args.base**i for i in range(start, args.end + 1))

    for length in rng:
        filename = "%s/%s_%s" % (args.output_dir, args.residue, length)
        ribfilename = "%s.rib" % filename

        with open(ribfilename, "w") as ribfile:
            ribfile.write("TITLE %s\n" % args.residue.lower())
            ribfile.write("DEFAULT PHI %.2f\nDEFAULT PSI %.2f\n" % (args.phi, args.psi))

            for res in range(length):
                ribfile.write("RES %s\n" % args.residue)

        if args.path_to_ribosome_dir:
            ribosome = "%s/ribosome" %  args.path_to_ribosome_dir
            res_zmat = "%s/res.zmat" %  args.path_to_ribosome_dir
            pdbfilename = "%s.pdb" % filename

            subprocess.call([ribosome, ribfilename, pdbfilename, res_zmat])
            # Set the occupancy to 1 to mark entire chain as flexible
            subprocess.call(['sed', '-ri', r"'s/(^ATOM.*)/\1 1.00/'", pdbfilename])

