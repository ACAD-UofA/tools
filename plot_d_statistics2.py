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
Plot a series of 4 population D statistics.  The series constists of sets of
3 tests, which correspond to the 3 possible phylogenetic topologies.

E.g.
$ ./plot_d_statistics.py dstats.tsv
"""

from __future__ import print_function
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FixedLocator

def parse_dstats_tsv(filename):
    data = []
    with open(filename) as f:
        next(f) # skip header
        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            p1 = fields[0]
            p2 = fields[1]
            p3 = fields[2]
            p4 = fields[3]
            d = float(fields[4])
            z = float(fields[5])
            baba = int(fields[6])
            abba = int(fields[7])
            n_snps = int(fields[8])

            data.append((p1, p2, p3, p4, d, z, baba, abba, n_snps))
    data.reverse()
    return data

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Plot D statistics")
    parser.add_argument("-p", "--parsimonious", action="store_true", default=False, help="Plot only the most parsimonious from the block of 3")
    parser.add_argument("-u", "--ungrouped", action="store_true", default=False, help="Do not group into blocks of 3")
    parser.add_argument("-g", "--greyscale", action="store_true", default=False, help="Use grey/black instead of red/black")
    parser.add_argument("-s", "--figscale", default=2, type=float, help="[%(default)s]")
    parser.add_argument("-o", "--outpdf", default="plot_d_statistics.pdf", help="output file [%(default)s]")
    parser.add_argument("tsv", help="D_stats.tsv")
    parser.add_argument("outgroup", help="Outgroup, to show in the title")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    dstats_raw = parse_dstats_tsv(args.tsv)

    dstats = []
    if args.parsimonious:
        args.ungrouped = True
        for i in range(len(dstats_raw)/3):
            baba0, abba0 = dstats_raw[i*3][6:8]
            baba1, abba1 = dstats_raw[i*3+1][6:8]
            baba2, abba2 = dstats_raw[i*3+2][6:8]
            m = np.max((baba0,abba0,baba1,abba1,baba2,abba2))
            if m not in (baba0, abba0):
                dstats.append(dstats_raw[i*3])
            elif m not in (baba1, abba1):
                dstats.append(dstats_raw[i*3+1])
            else: #m not in (baba2, abba2):
                dstats.append(dstats_raw[i*3+2])

    tt = "$" + args.outgroup.replace(" ", "$ $") + "$"

    pdf = PdfPages(args.outpdf)
    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig_w, fig_h = plt.figaspect(3.0/4.0)

    fig1 = plt.figure(figsize=(args.figscale*fig_w, args.figscale*fig_h))
    ax1 = fig1.add_subplot(111)

    y_pos = np.arange(len(dstats))
    p4_labels = ["(({}, {}), {})".format(p1,p2,p3) for p1,p2,p3,p4,_,_,_,_,_ in dstats]
    d = [dd for _,_,_,_,dd,_,_,_,_ in dstats]

    # 3 * standard error
    err = [3*abs(dd/zz) for _,_,_,_,dd,zz,_,_,_ in dstats]

    if args.greyscale:
        col0, col1 = ('grey', 'k')
    else:
        col0, col1 = ('r', 'k')

    colours = [col1 if abs(dd)-ee < 0 else col0 for dd, ee in zip(d, err)]

    if args.ungrouped:
        edgecolour = colours
    else:
        # pick one of the group of tests to have a different edge colour
        edgecolour = []
        for i in range(len(dstats)/3):
            d0 = abs(dstats[i*3][4])
            d1 = abs(dstats[i*3+1][4])
            d2 = abs(dstats[i*3+2][4])
            if d0 < d1 and d0 < d2:
                edgecolour.extend([col1, col0, col0])
            elif d1 < d0 and d1 < d2:
                edgecolour.extend([col0, col1, col0])
            else:
                edgecolour.extend([col0, col0, col1])

    ax1.scatter(d, y_pos, c=colours, edgecolors=edgecolour, lw=1.5, s=60, zorder=3)
    ax1.errorbar(d, y_pos, xerr=err, ecolor=col1, marker="none", fmt="none", zorder=0, capsize=0)

    ax1.axvline(c='k', ls=':')

    if not args.ungrouped:
        # horizonal lines
        for i in range(1, len(dstats)/3):
            ax1.axhline(3*i-0.5, c='k', ls=':')

    ax1.set_ylim([-0.5, y_pos[-1]+0.5])
    xlim = ax1.get_xlim()
    ax1.set_xlim([max(-1, xlim[0]), min(1, xlim[1])])

    ax1.set_yticks(y_pos)
    ylabels = ax1.set_yticklabels(p4_labels)
    for l, ec in zip(ylabels, edgecolour):
        l.set_color(ec)
    ax1.set_xlabel("$D$", labelpad=20, size=20)
    ax1.set_ylabel("(($P_1$, $P_2$), $P_3$)", labelpad=20, size=20)
    ax1.set_title("$D$((($P_1$, $P_2$), $P_3$), "+tt+")", y=1.01, size=20)
    if not args.ungrouped:
        ax1.yaxis.set_minor_locator(FixedLocator([3*i-0.5 for i in range(1, len(dstats)/3)]))
        ax1.yaxis.set_tick_params('minor', length=4, direction='out')

    # Do a bunch of coordinate transformations for some text...
    # I want the x coordinates to line up with the yaxis edges
    # and the y coordinate to be a small fraction of the figure size.
    x1_dcoords = ax1.transAxes.transform((0, 0))
    x2_dcoords = ax1.transAxes.transform((1, 0))
    xx1, _ = ax1.transData.inverted().transform(x1_dcoords)
    xx2, _ = ax1.transData.inverted().transform(x2_dcoords)
    y_dcoords = fig1.transFigure.transform((0, 0.01*(args.figscale**2)))
    _, yy = ax1.transData.inverted().transform(y_dcoords)

    if args.figscale >= 2:
        text_left = "suggests $P_2 \longleftrightarrow P_3$ affinity"
        text_right = "suggests $P_1 \longleftrightarrow P_3$ affinity"
    else:
        text_left = "$P_2 \longleftrightarrow P_3$"
        text_right = "$P_1 \longleftrightarrow P_3$"
    
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
