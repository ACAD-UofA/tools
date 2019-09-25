#!/usr/bin/env python
"""
Import Fst data as calculated using vcftools with the --weir-fst-pop parameter.

Computes Z-scores (number of std.devs. from the mean).
"""

from __future__ import print_function
import sys
import math
import numpy as np
from scipy import stats

def parse_tsv(filename, column):
    with open(filename) as f:
        # skip header
        next(f)

        data = []

        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            fst = float(fields[column])
            if math.isnan(fst):
                continue

            data.append(fst)

    return data

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Obtain Fst from a Tab Separated Value file and append a zFst column")
    parser.add_argument("-c", "--column", type=int, help="column containing Fst values")
    parser.add_argument("tsv_file", help=".tsv file")

    args = parser.parse_args()
    args.column -= 1 # make it zero based
    return args

if __name__ == "__main__":

    args = parse_args()

    fst = parse_tsv(args.tsv_file, args.column)

    # fst normalised by subtracting the mean and dividing std.dev.
    zfst = stats.zscore(fst)
    ziter = iter(zfst)

    with open(args.tsv_file) as f:
        header = next(f)

        header = header.rstrip("\r\n")
        print("{}\tzFst".format(header))

        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")
            fst = float(fields[args.column])
            if math.isnan(fst):
                print("{}\t0.0".format(line))
            else:
                z = next(ziter)
                print("{}\t{}".format(line, z))
