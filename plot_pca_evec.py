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
import re
import collections
import numpy as np
import numpy.random as npr
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_evec(filename):
    re_s = re.compile(r"[ \t]+")
    evec = collections.defaultdict(dict)
    eigvals = None
    with open(filename) as f:
        for line in f:
            fields = re_s.split(line.strip())
            if fields[0][0] == "#":
                eigvals = [float(x) for x in fields[1:]]
                continue
            sample = fields[0]
            group = fields[-1]
            evec[group][sample] = [float(x) for x in fields[1:-1]]
    return eigvals, evec

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Plot PCA from eigensoft *.pca.evec file")
    parser.add_argument("--groups", required=True, help="groups")
    parser.add_argument("evec_file", help="file.pca.evec")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    grouplist = args.groups.split(",")
    eigvals, evec = parse_evec(args.evec_file)

    pdf = PdfPages("plot_pca_evec.pdf")
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)

    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    pc1=0
    pc2=1

    pal = ("#a6cee3", "#b15928", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
            "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#1f78b4")

    for group, colour in zip(grouplist, pal):
        xx = [vec[pc1] for vec in evec[group].values()]
        yy = [vec[pc2] for vec in evec[group].values()]
        ax1.scatter(xx, yy, c=colour, edgecolors=colour, label=group)

        x_jitter = npr.random(len(evec[group])) / 100
        y_jitter = npr.random(len(evec[group])) / 100
        for (sample, vec), xx, yy in zip(evec[group].iteritems(), x_jitter, y_jitter):
            #http://stackoverflow.com/questions/8850142/matplotlib-overlapping-annotations
            ax1.annotate(sample, xy=(vec[pc1], vec[pc2]))#, xytext=(vec[pc1]+xx, vec[pc2]+yy))

    ax1.set_xlabel("PC{} ({:.2f}%)".format(pc1+1, eigvals[pc1]))
    ax1.set_ylabel("PC{} ({:.2f}%)".format(pc2+1, eigvals[pc2]))
    ax1.legend(numpoints=1, scatterpoints=1)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
