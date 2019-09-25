#!/usr/bin/env python

from __future__ import print_function
import sys
import collections
import os.path

import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np

from colours import pal16

def parse_snp(fn):
    snp_dict = {}
    with open(fn) as f:
        for lineno, line in enumerate(f,1):
            fields = line.split()
            #id = fields[0]
            chrom = fields[1]
            #gpos = fields[2]
            pos = int(fields[3])
            #ref = fields[4]
            #alt = fields[5]
            snp_dict[lineno] = (chrom, pos)
    return snp_dict

def parse_fai(fn):
    chrom_dict = {}
    with open(fn) as f:
        for line in f:
            fields = line.split()
            chrom = fields[0]
            size = int(fields[1])

            if chrom.startswith("chr"):
                chrom = chrom[3:]

            chrom_dict[chrom] = size

    return chrom_dict

def parse_saguaro_out(fn):
    with open(fn) as f:
        for lineno, line in enumerate(f,1):
            if not line.startswith("cactus"):
                continue
            fields = line.split()
            cactus = int(fields[0][6:])
            chrom = fields[1][:-1]
            if chrom.startswith("chr"):
                chrom = chrom[len("chr"):]
            start = int(fields[2])
            end = int(fields[4])
            score = float(fields[-1][len("score="):])
            yield lineno, cactus, chrom, start, end, score

def breakpoints2recomb(breakpoints, chlen, step, window):
    n_steps = int(np.ceil(float(chlen)/step))
    steps = np.empty(n_steps)
    x = 0
    for i, j in enumerate(range(step,chlen,step)):
        last_x = x
        while x < len(breakpoints) and breakpoints[x] < j:
            x += 1
        steps[i] = x-last_x

    # Account for partial step at the end.
    if n_steps > i:
        steps[-1] = len(breakpoints)-x

    # Normalise such that each step contains a value corresponding to
    # Pr(recomb | recomb occured somewhere on the chromosome).
    steps /= n_steps
    steps /= np.sum(steps)

    spw = window/step # steps per window

    r = np.empty(n_steps-spw+1, dtype=float)
    r[0] = sum(steps[:spw])
    for j in range(n_steps-spw):
        r[j+1] = r[j] - steps[j] + steps[j+spw]

    # average the values of the steps for each window
    r /= spw

    # midpoint of the windows
    windowpos = np.arange(window/2, chlen-window/2+step+1, step)

    return r, windowpos

def parse_args():
    import argparse

    def km_int(s):
        if s[-1] in "kK":
            s2 = s[:-1]
            m = 1000
        elif s[-1] in "mM":
            s2 = s[:-1]
            m = 1000000
        else:
            s2 = s
            m = 1
        try:
            x = int(s2)
        except:
            raise argparse.ArgumentTypeError("unrecognised size `{}'".format(s))
        return x*m

    parser = argparse.ArgumentParser(description="plot Saguaro output")
    parser.add_argument("-w", "--window", type=km_int, default="5m", help="use windows of size WINDOW [%(default)s]")
    parser.add_argument("-s", "--step", type=km_int, default="200k", help="slide windows with step size STEP [%(default)s]")
    parser.add_argument("-o", "--opdf", type=str, default="out.pdf", help="output filename [%(default)s]")
    parser.add_argument("fai", type=str, help="ref.fa index")
    parser.add_argument("infiles", nargs="+", help="input files")
    args = parser.parse_args()

    args.snp_files = []
    args.saguaro_files = []
    for fn in args.infiles:
        if fn.endswith(".snp"):
            args.snp_files.append(fn)
        elif fn.endswith(".out"):
            args.saguaro_files.append(fn)

    if len(args.snp_files) != len(args.saguaro_files):
        print("number of .snp files doesn't match number of .out files", file=sys.stderr)
        exit(1)

    if len(args.snp_files) == 0 or len(args.saguaro_files) == 0:
        print("missing .snp and/or .out file", file=sys.stderr)
        exit(1)

    return args

if __name__ == "__main__":
    args = parse_args()

    chrom_size = parse_fai(args.fai)

    pdf = PdfPages(args.opdf)
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    gs1 = gridspec.GridSpec(1, 1)

    ax = fig1.add_subplot(gs1[0])

    linestyles = [':', '--', '-.', '-']
    colours = pal16

    for col, ls, snp_fn, saguaro_fn in zip(colours, linestyles, args.snp_files, args.saguaro_files):
        snps = parse_snp(snp_fn)
        loci = collections.defaultdict(list)
        for lineno, _, chrom, start, end, _ in parse_saguaro_out(saguaro_fn):
            ch, pos = snps[start]
            assert (chrom == ch), "chrom `{}' doesn't match `{}' from {}, at {}: line {}".format(ch, chrom, snp_fn, saguaro_fn, lineno)
            loci[ch].append(pos)

        # last one; use end
        ch, pos = snps[end]
        assert (chrom == ch), "chrom `{}' doesn't match `{}' from {}, at {}: line {}".format(ch, chrom, snp_fn, saguaro_fn, lineno)
        loci[ch].append(pos)

        assert len(loci.keys()) == 1, "multi-chromosome files unsupported"

        breakpoints = loci[ch]
        chlen = chrom_size[ch]
        r, x = breakpoints2recomb(breakpoints, chlen, args.step, args.window)

        label = os.path.basename(saguaro_fn)
        ax.plot(x, r, color=col, linestyle=ls, label=label)

    xt_interval = 2e7
    xticks = np.arange(xt_interval, chlen, xt_interval)
    xlabels = [str(int(x/1e6))+" mbp" for x in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, size="6")

    ax.legend()

    ax.set_title("chr" + ch)

    ax.set_xlim([-0.5, chlen+0.5])

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
