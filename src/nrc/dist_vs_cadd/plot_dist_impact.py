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
    low = []
    modifier = []
    moderate = []
    high = []
#    sample1	sample2	distance	n_sites	impact_LOW	impact_MODIFIER	impact_MODERATE	impact_HIGH
    with open(filename) as f:
        next(f) # skip header
        for lineno, line in enumerate(f, 1):
            #if lineno % 1000 != 0:
            #    continue
            line = line.rstrip("\r\n")
            fields = line.split("\t")
            d = int(fields[2])
            n = int(fields[3])

            dist.append(float(d)/n) # mutations per base
            low.append(float(fields[4])/d)
            #modifier.append(float(fields[5])/d)
            #moderate.append(float(fields[6])/d)
            #high.append(float(fields[7])/d)
    return dist, low, modifier #, moderate, high

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="plot impact vs. distance")
    parser.add_argument("in_tsv", help="in.tsv")
    parser.add_argument("out_pdf", default="plot.pdf", help="out.pdf")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

#    d, low, modifier, moderate, high = parse_tsv(args.in_tsv)
    d, low, modifier = parse_tsv(args.in_tsv)
    low_mod = [l/m for l,m in zip(low,modifier)]

    pct = np.arange(0, 101, 1)
    dp = np.percentile(d, pct)
    dec_mean = []

    for i in range(len(pct)-2):
        a = [ll for dd, ll in zip(d, low) if dp[i] < dd < dp[i+1]]
        #b = [ll for dd, ll in zip(d, low) if dp[i+1] <= dd < dp[i+2]]
        dec_mean.append(np.mean(a))
        #u, p = mannwhitneyu(a, b)
        #print(pct[i+1], u, p)
        #u, p = mannwhitneyu(b, a)
        #print(u, p)
    i+=1
    a = [ll for dd, ll in zip(d, low) if dp[i] < dd < dp[i+1]]
    dec_mean.append(np.mean(a))

    pdf = PdfPages(args.out_pdf)
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)

    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    colours = ("#ca0020", "#f4a582", "#92c5de", "#0571b0")

#    for i, (impact, label) in enumerate(zip([low, modifier, low_mod], ["LOW", "MODIFIER", "LOW/MODIFIER"])):
#        if label in ["MODIFIER", "LOW/MODIFIER"]:
#            continue
#        ax1.plot(d, impact, ls="none", c=colours[i], mfc=colours[i], mec=colours[i], marker="o", ms=2, label=label)

    d_sub = [dd for i, dd in enumerate(d) if i%100==0]
    l_sub = [ll for i, ll in enumerate(low) if i%100==0]

    ax1.plot(d_sub, l_sub, ls="none", c=colours[2], mfc=colours[2], mec=colours[2], marker="o", ms=2, zorder=1, label="snpEff impact (LOW)")

    ymin, ymax = ax1.get_ylim()
    #ax1.vlines(dp, ymin, ymax, color=colours[2], zorder=2, label="Deciles")
    ax1.hlines(dec_mean, dp[:-1], dp[1:], color=colours[0], zorder=2, label="Percentile mean")

    ax1.set_ylim(ymin, ymax)

    ax1.set_xlabel("Pairwise differences per base")
    ax1.set_ylabel("Fraction of mutations")
    #ax1.legend(title="SnpEff Impact", numpoints=1, loc="center right")
    ax1.legend(numpoints=1)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
