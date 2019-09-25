#!/usr/bin/env python

from __future__ import print_function
import sys
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np

from colours import pal16

def parse_Q(fn, skip=0):
    rows = []
    rlen = None
    with open(fn) as f:
        while skip:
            next(f)
            skip -= 1
        for lineno, line in enumerate(f, skip):
            row = list(map(float, line.split()))
            if rlen is None:
                rlen = len(row)
            elif rlen != len(row):
                print("{}:{}: row length ({}) doesn't match first row ({})".format(fn, lineno, len(row), rlen), file=sys.stderr)
                exit(1)
            rows.append(row)
    return np.array(rows)

def parse_ind(fn):
    indlist = []
    with open(fn) as f:
        for line in f:
            fields = line.split()
            ind = "/".join([fields[2], fields[0]])
            indlist.append(ind)
    return indlist

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="barplot of ADMIXTURE Q matrix")
    parser.add_argument("-w", "--wide", action="store_true", default=True, help="plot widescreen ratio (16x9) [%(default)s]")
    parser.add_argument("-o", "--opdf", type=str, default="out.pdf", help="output filename [%(default)s]")
    parser.add_argument("-i", "--indfile", help="file containing list of individuals in first column, matching Q matrix rows")
    parser.add_argument("Qfile", help="input file")
    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = parse_args()

    Q = parse_Q(args.Qfile)
    n_inds = Q.shape[0]

    if args.indfile:
        indnames = parse_ind(args.indfile)
        if len(indnames) != n_inds:
            print("{} has {} individuals and {} has {} rows".format(args.indfile, len(indnames), args.Qfile, n_inds))
            exit(1)

    pdf = PdfPages(args.opdf)
    if args.wide:
        aspect = 9.0/16.0
    else:
        aspect = 3.0/4.0
    fig_w, fig_h = plt.figaspect(aspect)
    #fig_w, fig_h = 120,5 #plt.figaspect(aspect)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    gs1 = gridspec.GridSpec(1, 1)
    ax1 = fig1.add_subplot(gs1[0])

    indnums = np.arange(n_inds)
    bottom = np.vstack((np.zeros((n_inds,), dtype=Q.dtype), np.cumsum(Q.T, axis=0)[:-1]))

    for dat, col, bot in zip(Q.T, reversed(pal16), bottom):
        ax1.bar(indnums+0.5, dat, color=col, bottom=bot, alpha=1.0, width=1,
                edgecolor="black", linewidth=1)

    ax1.set_xticks(np.arange(0.5,n_inds))
    ax1.tick_params(length=0)
    if args.indfile:
        ax1.set_xticklabels(indnames, rotation='vertical', fontsize=6)
    ax1.get_yaxis().set_visible(False)
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(0.0, n_inds)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
