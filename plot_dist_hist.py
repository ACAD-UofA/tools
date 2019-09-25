#!/usr/bin/env python

from __future__ import print_function
import sys
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import parsers

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} in.{{fa,fq}}".format(sys.argv[0]), file=sys.stderr)

    dist = []

    for seq in parsers.parse_fa(sys.argv[1]):
        label = seq[0]
        extra = label.split(" ")[1]
        for ef in extra.split(";"):
            k, v = ef.split("=")
            if k == "DIST":
                dist.append(int(v))

    pdf = PdfPages("plot_dist_hist.pdf")
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    ax1.hist(dist, bins=range(0,max(dist)+1))
    ax1.set_xlabel("Edit Distance")
    ax1.set_ylabel("Frequency")
    ax1.set_title("Histogram of Edit Distances")

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
