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
import parsers
import numpy as np
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

font = {
        'family' : 'sans-serif',
        'sans-serif': 'Computer Modern Sans serif',
        'size'   : 11}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Plot fragment length histograms of reads")
    parser.add_argument("--title", default=None, help="label for the title")
    parser.add_argument("--max_len", type=int, default=10000, help="reads longer than this are ignored")
    parser.add_argument("in_fa", help="in.{fa,fq}")
    parser.add_argument("out_pdf", default="plot_fraglen_hist.pdf", help="out.pdf")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    lengths = np.zeros(args.max_len)
    lmin = 1e6
    lmax = 0

    for line in parsers.parse_fa(args.in_fa):
        length = len(line[1])
        if length > args.max_len:
            continue

        lengths[length] += 1

        if length < lmin:
            lmin = length
        elif length > lmax:
            lmax = length

    pdf = PdfPages(args.out_pdf)
    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig_w = fig_h = 307.28987/72.27

    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    colour = "blue"
    bins =  np.arange(lmin, lmax+1) -0.5
    ll = lengths[lmin:lmax+1] / sum(lengths[lmin:lmax+1])
    ax1.bar(bins, ll, width=0.5, color=colour, edgecolor=colour)
    ax1.set_xlabel("Length (bp)")
    ax1.set_ylabel("Frequency")

    if args.title:
        title = "Fragment length histogram"
        title += " ({})".format(args.title)
        ax1.set_title(title)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
