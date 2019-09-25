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
    hist = []
    fits = collections.defaultdict(dict)
    with open(filename) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if len(line) == 0:
                continue
            if line[0] == "#":
                fields = line[1:].split(";")
                for ff in fields[1:]:
                    k, v = ff.split("=")
                    fits[fields[0]][k] = float(v)
                continue
            hist.append(int(line))
    return np.array(hist, dtype=float), fits

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Plot fragment length histogram")
    parser.add_argument("-m", "--maxlen", type=int, default=0, help="maximum length to plot")
    parser.add_argument("-t", "--title", help="title")
    parser.add_argument("-o", "--outpdf", default="plot_length_distrib.pdf", help="out.pdf")
    parser.add_argument("--hist", nargs='*', help="in.hist", required=True)
    parser.add_argument("-l", "--label", nargs='*', help="hist labels")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    if len(args.label) != 0 and len(args.hist) != len(args.label):
        print("--label list must match the --hist input files".format())
        exit(1)

    hist = []
    for fn in args.hist:
        h, _ = parse_lengths(fn)
        hist.append(h)

    for i,_ in enumerate(hist):
        area = np.sum(hist[i])
        #hist[i] /= area

    if args.maxlen:
        for i,_ in enumerate(hist):
            hist[i] = hist[i][:args.maxlen]

    pdf = PdfPages(args.outpdf)
    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    # black, orange, sky blue, blueish green, yellow, dark blue, vermillion, reddish purple
#    pallete = ("#000000", "#906000", "#357090", "#006050", "#959025", "#004570", "#804000", "#806070")
    # blue/green, orange, purple
    pallete = ("#000000", "#1b9e77", "#d95f02", "#7570b3")
    style = (":", "-", "--", ":")

    width=0.8
    for i, (h, col, ls) in enumerate(zip(hist, pallete, style)):
        bins =  np.arange(1, len(h)+1)
        if len(args.label) != 0:
            label = args.label[i]
        ax1.plot(bins, h, color=col, ls=ls, lw=2, label=label)

    ax1.xaxis.set_tick_params(direction="inout")
    ax1.set_xlabel("Length (bp)")
    ax1.set_ylabel("Frequency")

    if len(args.label) != 0:
        ax1.legend(numpoints=1)

    if args.title:
        ax1.set_title(args.title)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
