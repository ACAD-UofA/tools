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
    with open(filename) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if len(line) == 0:
                continue
            hist.append(int(line))
    return np.array(hist, dtype=float)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Plot coverage histogram")
    parser.add_argument("-w", "--winsz", type=int, default=0, help="window size")
    parser.add_argument("--title", help="title")
    parser.add_argument("in_hist", help="in.hist")
    parser.add_argument("out_pdf", default="plot_coverage.pdf", help="out.pdf")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    hist = parse_lengths(args.in_hist)

    area = np.sum(hist)
    hist /= area
    cdf = np.zeros(len(hist))
    pc = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
    pcx = np.zeros(len(pc))
    j = 0
    for i, x in enumerate(hist):
        if i == 0:
            cdf[0] = 0
        else:
            cdf[i] = cdf[i-1] + hist[i-1]
        if j<len(pc) and cdf[i] > pc[j]:
            pcx[j] = i
            j += 1

    cov = np.arange(1.0, len(hist)+1.0)
    if args.winsz:
         cov /= args.winsz
         pcx /= args.winsz

    pallete = ("#000000", "#1b9e77", "#d95f02", "#7570b3")
    pdf = PdfPages(args.out_pdf)
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)
    ax2 = ax1.twinx()

    l1, = ax1.plot(cov, hist*area, color=pallete[0], lw=2, label="Frequency")
    l2, = ax2.plot(cov, cdf, ls="--", color=pallete[1], lw=2, label="Cumulative")
    for tl in ax2.get_yticklabels():
        tl.set_color(pallete[1])

    l3 = ax2.vlines(pcx[:-1], 0, 1, colors=pallete[3], linestyles=":", label="Deciles")
    l4 = ax2.vlines(pcx[-1], 0, 1, colors=pallete[2], linestyles="-.", label="99th %ile")

    ax1.set_xlim(0, pcx[-1]*1.3)
    ax2.set_ylim(0, 1)
    ax1.xaxis.set_tick_params(direction="inout")

    ax1.set_xlabel("Coverage (averaged over the window)")
    ax1.set_ylabel("Frequency", color=pallete[0])
    ax2.set_ylabel("Cumulative", color=pallete[1])
    ax1.legend(handles=[l1, l2, l3, l4], numpoints=1, loc="center right")

    if args.title:
        ax1.set_title(args.title)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
