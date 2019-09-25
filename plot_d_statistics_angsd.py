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

"""
Plot a series of 4 population D statistics.

E.g.
$ ./plot_d_statistics.py dstats.tsv
"""

from __future__ import print_function
import sys
import operator
import re
import numpy as np
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FixedLocator

def parse_dstats_tsv(filename, pop1, pop2, pop3):

    if pop1:
        re_p1 = [re.compile(p) for p in pop1]
    if pop2:
        re_p2 = [re.compile(p) for p in pop2]
    if pop3:
        re_p3 = [re.compile(p) for p in pop3]

    data = []
    with open(filename) as f:
        next(f) # skip header
        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            p1 = fields[0]
            p2 = fields[1]
            p3 = fields[2]
            abba = int(fields[3])
            baba = int(fields[4])
            d = float(fields[5])
            #jackest = float(fields[6])
            se = float(fields[7])
            z = float(fields[8])

            if pop1:
                for r in re_p1:
                    if r.search(p1):
                        break
                else:
                    continue

            if pop2:
                for r in re_p2:
                    if r.search(p2):
                        break
                else:
                    continue

            if pop3:
                for r in re_p3:
                    if r.search(p3):
                        break
                else:
                    continue

            data.append((p1, p2, p3, None, d, z, baba, abba, None))
    data.reverse()
    return data

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Plot D statistics")
    parser.add_argument("-s", "--figscale", default=1.5, type=float, help="[%(default)s]")
    parser.add_argument("-o", "--outpdf", default="plot_d_statistics.pdf", help="output file [%(default)s]")
    parser.add_argument("-1", "--pop1", action="append", help="filter p1")
    parser.add_argument("-2", "--pop2", action="append", help="filter p2")
    parser.add_argument("-3", "--pop3", action="append", help="filter p3")
    parser.add_argument("tsv", help="D_stats.tsv")
    #parser.add_argument("outgroup", help="Outgroup, to show in the title")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    dstats = parse_dstats_tsv(args.tsv, args.pop1, args.pop2, args.pop3)
    dstats.sort(key=operator.itemgetter(4))

    if len(dstats) == 0:
        raise Exception("empty dstats table, check population filters")

    p2set = set(p2 for _,p2,_,_,_,_,_,_,_ in dstats)
    p3set = set(p3 for _,_,p3,_,_,_,_,_,_ in dstats)

    title_p4 = "Turkey"

    if len(p3set) == 1:
        if len(p2set) == 1:
            p4_labels = ["{}".format(p1) for p1,p2,p3,p4,_,_,_,_,_ in dstats]
            ylabel = "$P_1$"
            title_p2 = p2set.pop()
            title_p3 = p3set.pop()
        else:
            p4_labels = ["{}, {}".format(p1,p2) for p1,p2,p3,p4,_,_,_,_,_ in dstats]
            title_p2 = "$P_2$"
            title_p3 = p3set.pop()
            ylabel = "$P_1$, $P_2$"
    else:
        p4_labels = ["({}, {}; {})".format(p1,p2,p3) for p1,p2,p3,p4,_,_,_,_,_ in dstats]
        title_p2 = "$P_2$"
        title_p3 = "$P_3$"
        ylabel = "($P_1$, $P_2$; $P_3$)"

    d = [dd for _,_,_,_,dd,_,_,_,_ in dstats]

    # 3 * standard error
    err = [3*abs(dd/zz) for _,_,_,_,dd,zz,_,_,_ in dstats]

    pdf = PdfPages(args.outpdf)
    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig_w, fig_h = plt.figaspect(3.0/4.0)

    fig1 = plt.figure(figsize=(args.figscale*fig_w, args.figscale*fig_h))
    ax1 = fig1.add_subplot(111)

    col0, col1 = ('k', 'k')

    colours = [col1 if abs(dd)-ee < 0 else col0 for dd, ee in zip(d, err)]
    edgecolour = colours

    y_pos = np.arange(len(dstats))
    ax1.scatter(d, y_pos, c=colours, edgecolors=edgecolour, lw=1.5, s=60)
    ax1.errorbar(d, y_pos, xerr=err, ecolor=col1, marker="none", fmt="none", capsize=0)

    ax1.set_ylim([-0.5, y_pos[-1]+0.5])
    xlim = ax1.get_xlim()
    if xlim[0] < 0 and xlim[1] < 0:
        xlim = (xlim[0], 0)
    elif xlim[0] > 0 and xlim[1] > 0:
        xlim = (0, xlim[1])
    ax1.set_xlim([max(-1, xlim[0]), min(1, xlim[1])])
    xlim = ax1.get_xlim()

    ax1.hlines(y_pos, xlim[0], xlim[1], colors='lightgrey', lw=0.1)
    ax1.axvline(c='k', ls=':')

    ax1.set_yticks(y_pos)
    ylabels = ax1.set_yticklabels(p4_labels)
    for l, ec in zip(ylabels, edgecolour):
        l.set_color(ec)
    ax1.set_xlabel("$D$", labelpad=20, size=args.figscale*10)
    ax1.set_ylabel(ylabel, labelpad=20, size=args.figscale*10)
    ax1.set_title("$D$($P_1$, {}; {}, {})".format(title_p2, title_p3, title_p4), y=1.01, size=args.figscale*10)

    # Do a bunch of coordinate transformations for some text...
    # I want the x coordinates to line up with the yaxis edges
    # and the y coordinate to be a small fraction of the figure size.
    x1_dcoords = ax1.transAxes.transform((0, 0))
    x2_dcoords = ax1.transAxes.transform((1, 0))
    xx1, _ = ax1.transData.inverted().transform(x1_dcoords)
    xx2, _ = ax1.transData.inverted().transform(x2_dcoords)
    y_dcoords = fig1.transFigure.transform((0, 0.01*(args.figscale**2)))
    _, yy = ax1.transData.inverted().transform(y_dcoords)

    # ABBA-BABA convention
    if xlim[0] < 0:
        text_left = "$P_1$ $\longleftrightarrow$ "+title_p3
    else:
        text_left = ""
    if xlim[1] > 0:
        text_right = title_p2+" $\longleftrightarrow$ "+title_p3
    else:
        text_right = ""
    
    ax1.annotate(text_left,
            xy=(0, 0), xytext=(xx1, yy),
            va='bottom', ha='left',
            annotation_clip=False)
    ax1.annotate(text_right,
            xy=(0, 0), xytext=(xx2, yy),
            va='bottom', ha='right',
            annotation_clip=False)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
