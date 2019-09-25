#!/usr/bin/env python

from __future__ import print_function
import sys
import collections
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np

def parse_log(fn, reads):
    rset = set(reads)
    data = collections.defaultdict(lambda: [[], []])
    state = 0
    sumll = {}
    with open(fn) as f:
        for line in f:
            if line[0] == '@':
                read = line.strip()[1:]
                if read in rset:
                    state = 1
                else:
                    state = 0
                continue

            if state == 0:
                continue

            fields = line.split()
            ll = float(fields[2])

            if fields[1] == "sum":
                sumll[fields[0]] = ll
                continue

            if fields[0] == "HP":
                data[read][0].append(ll+sumll["HP"])
            elif fields[0] == "YY":
                data[read][1].append(ll+sumll["YY"])
    return data

if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print("usage: {} ll.log read1 [... readN]".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    reads = sys.argv[2:]
    data = parse_log(sys.argv[1], reads)

    with PdfPages("ll.pdf") as pdf:
        #fig_w, fig_h = plt.figaspect(9.0/16.0)
        fig_w, fig_h = plt.figaspect(3.0/4.0)

        for read in reads:
            hp, yy = data[read]

            fig1 = plt.figure(figsize=(fig_w, fig_h))
            gs1 = gridspec.GridSpec(1, 1)
            ax1 = fig1.add_subplot(gs1[0])

            lblue="#9dafff"
            cyan="#29d0d0"

            ax1.scatter(np.arange(1,len(yy)+1), yy, color=lblue, marker="o", label="m=YY")
            ax1.hlines(yy[-1], 1, len(yy)+1, color=lblue, linestyle=":")

            ax1.scatter(np.arange(1,len(hp)+1), hp, color="k", marker="x", label="m=HP")
            ax1.hlines(hp[-1], 1, len(yy)+1, color="k", linestyle=":")

            ax1.set_title(read)
            ax1.set_xlabel("molecule length (bp)")
            ax1.set_ylabel("log[P(data|LEN=i,MOL=m)]")
            ax1.legend(loc="lower right", scatterpoints=1)

            ax1.set_xlim([-3, len(yy)+3])
            ax1.set_ylim([1.05*min(min(hp),min(yy)), 20])

            plt.tight_layout()
            pdf.savefig()
            plt.close()
