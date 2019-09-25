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
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_tsv(filename):
    dist = []
    num = []
    low = []
    modifier = []
    moderate = []
    high = []
#    sample1	sample2	distance	n_sites	impact_LOW	impact_MODIFIER	impact_MODERATE	impact_HIGH
    with open(filename) as f:
        next(f) # skip header
        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")
            d = int(fields[2])
            n = int(fields[3])
            if float(d)/n > 0.5:
                # ignore outlier
                continue

            dist.append(d)
            num.append(n)
            low.append(float(fields[4])/d)
            modifier.append(float(fields[5])/d)
            moderate.append(float(fields[6])/d)
            high.append(float(fields[7])/d)
    return dist, num, low, modifier, moderate, high

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="plot impact vs. distance")
    parser.add_argument("in_tsv", help="in.tsv")
    parser.add_argument("out_pdf", default="plot.pdf", help="out.pdf")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    d, n, low, modifier, moderate, high = parse_tsv(args.in_tsv)

    pdf = PdfPages(args.out_pdf)
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)

    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    colours = ("#ca0020", "#f4a582", "#92c5de", "#0571b0")

    for i, (impact, label) in enumerate(zip([low, modifier, moderate, high], ["LOW", "MODIFIER", "MODERATE", "HIGH"])):
        ax1.plot(d, impact, c=colours[i], mfc=colours[i], mec=colours[i], marker="o", ms=2, label=label)

    ax1.set_xlabel("Pairwise distance")
    ax1.set_ylabel("Fraction of mutations")
    ax1.legend(title="SnpEff Impact")

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
