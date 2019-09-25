#!/usr/bin/env python
# Copyright (c) 2016 Graham Gower <graham.gower@gmail.com>
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

import sys
import collections
import numpy as np
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

if __name__ == "__main__":

    header = next(sys.stdin).rstrip("\r\n").split("\t")
    lengths = [int(l) for l in header[1:]]

    lengths = lengths[3:]

    counts = collections.defaultdict(list)
    nreps = 100

    for line in sys.stdin:
        line = line.rstrip("\r\n")
        fields = line.split("\t")
        for l, c_str in zip(lengths, fields[1:]):
            counts[l].append(float(c_str)/nreps)

    pdf = PdfPages("min_length_plot.pdf")
    fig_w, fig_h = plt.figaspect(9.0/(2*16.0))
    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    #fig1 = plt.figure(figsize=(4, 2))
    ax1 = fig1.add_subplot(111)

    ax1.violinplot([counts[l] for l in lengths], lengths, widths=1, showmedians=True, showextrema=False)

    ax1.set_xlabel("Length (bp)")
    ax1.set_ylabel("Proportion of unique mappings")
    #ax1.legend(numpoints=1)

    ax1.invert_xaxis()

    ax1.set_xlim((43,12))
    #ax1.set_xscale('log')
    ax1.set_ylim((-0.1,1.1))

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
