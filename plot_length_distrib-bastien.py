#!/usr/bin/env python

import sys
import os.path
import collections
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_lengths(filename):
    counts = []
    lengths = []
    with open(filename) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if len(line) == 0:
                continue
            fields = line.split()
            count = int(fields[0])
            length = int(fields[1])
            counts.append(count)
            lengths.append(length)
    return np.array(counts, dtype=float), np.array(lengths)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Plot fragment length histogram")
    #parser.add_argument("-m", "--maxlen", type=int, default=0, help="maximum length to plot")
    #parser.add_argument("-t", "--title", help="title")
    parser.add_argument("-o", "--outpdf", default="plot_length_distrib.pdf", help="out.pdf")
    parser.add_argument("--hist", nargs='*', help="in.hist", required=True)
    #parser.add_argument("-l", "--label", nargs='*', help="hist labels")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    pdf = PdfPages(args.outpdf)
    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    # black, orange, sky blue, blueish green, yellow, dark blue, vermillion, reddish purple
#    pallete = ("#000000", "#906000", "#357090", "#006050", "#959025", "#004570", "#804000", "#806070")

    pallete = ("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
                "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
                "#ffff99", "#b15928")

# blue/green, orange, purple
#    pallete = ("#000000", "#1b9e77", "#d95f02", "#7570b3")
#    style = (":", "-", "--", ":")

    colours = iter(pallete)

    for fn in args.hist:
        c,l = parse_lengths(fn)
        c = c/sum(c) # normalise by total counts
        ax1.plot(l, c, color=next(colours), lw=2, label=fn)

    ax1.xaxis.set_tick_params(direction="inout")
    ax1.set_xlabel("Length (bp)")
    ax1.set_ylabel("Frequency")

    ax1.legend(numpoints=1, loc="upper left")

#    if args.title:
#        ax1.set_title(args.title)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
