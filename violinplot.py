#!/usr/bin/env python

from __future__ import print_function
import sys
import os.path
import collections
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np


def parse_input(args):
    ylists = collections.defaultdict(list)
    with open(args.infile) as f:
        while args.skip > 0:
            next(f)
            args.skip -= 1
        for line in f:
            fields = line.split()
            if len(fields) < max(args.horiz_column, args.vert_column):
                continue
            x = float(fields[args.horiz_column-1])
            y = float(fields[args.vert_column-1])
            if (np.isnan(x) or np.isinf(x)) or (np.isnan(y) or np.isinf(y)):
                continue
            ylists[x].append(y)

    return sorted(ylists.keys()), ylists

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="plot histogram for 1 column of input file")
    parser.add_argument("-x", "--horiz-column", type=int, required=True, help="1-based column to plot on horizontal axis")
    parser.add_argument("-y", "--vert-column", type=int, required=True, help="1-based column to plot on vertical axis")
    parser.add_argument("-s", "--skip", type=int, default=0, help="number of input lines to skip [%(default)s]")
    parser.add_argument("-w", "--wide", action="store_true", default=False, help="plot widescreen ratio (16x9) [%(default)s]")
    parser.add_argument("--logx", action="store_true", default=False, help="plot x-axis on log scale [%(default)s]")
    parser.add_argument("--logy", action="store_true", default=False, help="plot y-axis on log scale [%(default)s]")
    parser.add_argument("--title", type=str, help="plot title")
    parser.add_argument("--xlabel", type=str, help="x axis label")
    parser.add_argument("--ylabel", type=str, help="y axis label")
    parser.add_argument("-o", "--opdf", type=str, default="out.pdf", help="output filename [%(default)s]")
    parser.add_argument("infile", help="input file")
    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_args()
    x, y = parse_input(args)

    pdf = PdfPages(args.opdf)
    if args.wide:
        fig_w, fig_h = plt.figaspect(9.0/16.0)
    else:
        fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    gs1 = gridspec.GridSpec(1, 1)
    ax1 = fig1.add_subplot(gs1[0])

    if not args.title:
        args.title = os.path.basename(args.infile)

    if not args.ylabel:
        if args.logy:
            args.ylabel = "log(column {})".format(args.vert_column)
        else:
            args.ylabel = "column {}".format(args.vert_column)

    if not args.xlabel:
        if args.horiz_column > 0:
            if args.logx:
                args.xlabel = "log(column {})".format(args.horiz_column)
            else:
                args.xlabel = "column {}".format(args.horiz_column)
        else:
            if args.logx:
                args.xlabel = "log(line number)"
            else:
                args.xlabel = "line number"

    ax1.violinplot([y[i] for i in x], x, points=60, widths=0.7, showmedians=True, bw_method=0.5)

    ax1.set_title(args.title)
    ax1.set_xlabel(args.xlabel)
    ax1.set_ylabel(args.ylabel)

    if args.logx:
        plt.xscale('log', nonposy='clip')
    if args.logy:
        plt.yscale('log', nonposy='clip')

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
