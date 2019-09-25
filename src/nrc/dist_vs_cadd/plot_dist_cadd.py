#!/usr/bin/env python
# Copyright (c) 2015 Graham Gower <graham.gower@gmail.com>
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

from __future__ import print_function
import sys
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_tsv(filename):
    dist = []
    cadd = []
    with open(filename) as f:
        next(f) # skip header
        for lineno, line in enumerate(f, 1):
#            if lineno % 100 != 0:
#                continue
            line = line.rstrip("\r\n")
            fields = line.split("\t")
            d = int(fields[2])
            n = int(fields[3])

            c_mean = sum(i*int(c) for i,c in enumerate(fields[8:],1)) / float(d)

            dist.append(float(d)/n) # mutations per base
            cadd.append(c_mean)
    return dist, cadd

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="plot impact vs. distance")
    parser.add_argument("in_tsv", help="in.tsv")
    parser.add_argument("out_pdf", default="plot.pdf", help="out.pdf")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    d, cadd = parse_tsv(args.in_tsv)

    pct = np.arange(0, 101, 1)
    dp = np.percentile(d, pct)
    dec_mean = []

    for i in range(len(pct)-2):
        a = [cc for dd, cc in zip(d, cadd) if dp[i] < dd < dp[i+1]]
        dec_mean.append(np.mean(a))
        #u, p = mannwhitneyu(a, b)
        #print(pct[i+1], u, p)
        #u, p = mannwhitneyu(b, a)
        #print(u, p)
    i+=1
    a = [cc for dd, cc in zip(d, cadd) if dp[i] < dd < dp[i+1]]
    dec_mean.append(np.mean(a))

    pdf = PdfPages(args.out_pdf)
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)

    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    colours = ("#ca0020", "#f4a582", "#92c5de", "#0571b0")

    d_sub = [dd for i, dd in enumerate(d) if i%100==0]
    c_sub = [cc for i, cc in enumerate(cadd) if i%100==0]

    ax1.plot(d_sub, c_sub, ls="none", c=colours[2], mfc=colours[2], mec=colours[2], marker="o", ms=2, zorder=1, label="mean C-score")

    ymin, ymax = ax1.get_ylim()
    #ax1.vlines(dp, ymin, ymax, color=colours[2], zorder=2, label="Deciles")
    ax1.hlines(dec_mean, dp[:-1], dp[1:], color=colours[0], zorder=2, label="Percentile mean")

    ax1.set_ylim(ymin, ymax)

    ax1.set_xlabel("Pairwise differences per base")
    ax1.set_ylabel("C-score")
    #ax1.legend(title="SnpEff Impact", numpoints=1, loc="center right")
    ax1.legend(numpoints=1)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
